#include "csr_graph.h"

#include <stdlib.h>
#include <sys/mman.h>

static inline void parse_id(char *data, size_t *p, long long *v)
{
    while (data[*p] < '0' || data[*p] > '9')
        (*p)++;

    *v = 0;
    while (data[*p] >= '0' && data[*p] <= '9')
        *v = (*v) * 10 + data[(*p)++] - '0';
}

csr_graph graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    NodeID N, t;
    EdgeID M;
    parse_id(data, &p, &N);
    parse_id(data, &p, &M);
    parse_id(data, &p, &t);

    EdgeID *V = malloc(sizeof(EdgeID) * (N + 1));
    NodeID *E = malloc(sizeof(NodeID) * (M * 2));

    NodeWeight *W = malloc(sizeof(NodeWeight) * N);

    EdgeID ei = 0;
    for (NodeID u = 0; u < N; u++)
    {
        parse_id(data, &p, W + u);
        V[u] = ei;
        while (ei < M * 2)
        {
            while (data[p] == ' ')
                p++;
            if (data[p] == '\n')
                break;
            
            EdgeID e;
            parse_id(data, &p, &e);
            E[ei++] = e - 1;;
        }
        p++;
    }
    V[N] = ei;

    munmap(data, size);

    return (csr_graph){.N = N, .V = V, .E = E, .W = W};
}

void graph_free(csr_graph g)
{
    free(g.V);
    free(g.E);
    free(g.W);
}

static inline int compare(const void *a, const void *b)
{
    return (*(NodeID *)a - *(NodeID *)b);
}

int graph_validate(NodeID N, const EdgeID *V, const NodeID *E)
{
    EdgeID M = 0;
    for (NodeID u = 0; u < N; u++)
    {
        if (V[u + 1] - V[u] < 0)
            return 0;

        M += V[u + 1] - V[u];

        for (EdgeID i = V[u]; i < V[u + 1]; i++)
        {
            if (i < 0 || i >= V[N])
                return 0;

            NodeID v = E[i];
            if (v < 0 || v >= N || v == u || (i > V[u] && v <= E[i - 1]))
                return 0;

            if (bsearch(&u, E + V[v], V[v + 1] - V[v], sizeof(int), compare) == NULL)
                return 0;
        }
    }

    if (M != V[N])
        return 0;

    return 1;
}

csr_graph graph_subgraph(csr_graph g, NodeID *mask, NodeID *reverse_map)
{
    int *forward_map = malloc(sizeof(int) * g.N);
    NodeID N = 0; 
    EdgeID M = 0;
    for (NodeID u = 0; u < g.N; u++)
    {
        if (!mask[u])
            continue;

        forward_map[u] = N;
        reverse_map[N] = u;
        N++;

        for (EdgeID i = g.V[u]; i < g.V[u + 1]; i++)
            if (mask[g.E[i]])
                M++;
    }

    csr_graph sg = {.N = N};
    sg.V = malloc(sizeof(EdgeID) * (g.N + 1));
    sg.E = malloc(sizeof(NodeID) * M);
    sg.W = malloc(sizeof(NodeWeight) * g.N);

    M = 0;
    for (NodeID u = 0; u < g.N; u++)
    {
        if (!mask[u])
            continue;

        sg.W[forward_map[u]] = g.W[u];
        sg.V[forward_map[u]] = M;

        for (EdgeID i = g.V[u]; i < g.V[u + 1]; i++)
        {
            NodeID v = g.E[i];
            if (!mask[v])
                continue;

            sg.E[M] = forward_map[v];
            M++;
        }
    }
    sg.V[sg.N] = M;
    free(forward_map);
    return sg;
}