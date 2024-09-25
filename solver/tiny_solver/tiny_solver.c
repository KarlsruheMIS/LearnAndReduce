#include "tiny_solver.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <limits.h>

#define MAX_NODES 128

tiny_solver *tiny_solver_init(NodeID N)
{
    tiny_solver *solver = malloc(sizeof(tiny_solver));

    solver->subgraph = malloc(sizeof(char *) * MAX_NODES);
    char *data = aligned_alloc(32, sizeof(char) * MAX_NODES * MAX_NODES);
    for (NodeID i = 0; i < MAX_NODES; i++)
        solver->subgraph[i] = data + i * MAX_NODES;

    solver->subgraph_W = aligned_alloc(32, sizeof(NodeWeight) * MAX_NODES);

    solver->forward_map = malloc(sizeof(NodeID) * N);
    for (NodeID i = 0; i < N; i++)
        solver->forward_map[i] = 0;
    solver->reverse_map = malloc(sizeof(NodeID) * MAX_NODES);
    for (NodeID i = 0; i < MAX_NODES; i++)
        solver->reverse_map[i] = 0;

    solver->independent_set = malloc(sizeof(NodeID) * MAX_NODES);

    char *solution_data = aligned_alloc(32, sizeof(char) * MAX_NODES * MAX_NODES);
    solver->per_layer_solution = malloc(sizeof(char *) * MAX_NODES);
    for (NodeID i = 0; i < MAX_NODES; i++)
        solver->per_layer_solution[i] = solution_data + i * MAX_NODES;
    solver->per_layer_weight = aligned_alloc(32, sizeof(NodeWeight) * MAX_NODES);
    solver->branch = malloc(sizeof(NodeID) * MAX_NODES);

    return solver;
}

void tiny_solver_free(tiny_solver *solver)
{
    free(*solver->subgraph);
    free(solver->subgraph);
    free(solver->subgraph_W);
    free(solver->forward_map);
    free(solver->reverse_map);
    free(solver->independent_set);
    free(*solver->per_layer_solution);
    free(solver->per_layer_solution);
    free(solver->per_layer_weight);
    free(solver->branch);
    free(solver);
}

void tiny_solver_clear_solution(tiny_solver *solver)
{
    solver->independent_set_weight = 0;
    for (int i = 0; i < MAX_NODES; i++)
        solver->independent_set[i] = 0;

    solver->time_limit_exceeded = 0;
    solver->weight_limit_exceeded = 0;
    solver->node_limit_exceeded = 0;
}
void tiny_solver_clear(tiny_solver *solver)
{
    solver->subgraph_N = 0;
    memset(*solver->subgraph, 0, sizeof(char) * MAX_NODES * MAX_NODES);
    tiny_solver_clear_solution(solver);
}

static inline void tiny_solver_reduce(tiny_solver *solver, int layer)
{
    char *S = solver->per_layer_solution[layer];
    char **g = solver->subgraph;
    NodeWeight *W = solver->subgraph_W;
    NodeID N = solver->subgraph_N;
    int imp = 1;

    while (imp)
    {
        imp = 0;
        // Trivial rule
        for (NodeID u = 0; u < N; u++)
        {
            if (S[u] == 0 && W[u] == 0)
            {
                imp = 1;
                S[u] = -1;
            }
        }

        // Neighbourhood rule
        for (NodeID u = 0; u < N; u++)
        {
            if (S[u] != 0)
                continue;

            NodeWeight nw = 0;
            for (NodeID v = 0; v < N; v++)
                if (v != u && g[u][v] && S[v] == 0)
                    nw += W[v];

            if (nw <= W[u])
            {
                imp = 1;
                solver->per_layer_weight[layer] += W[u];
                S[u] = 1;
                for (NodeID v = 0; v < N; v++)
                    if (u != v && g[u][v])
                        S[v] = -1;
            }
        }
        if (imp)
            continue;

        // Domination rule
        for (NodeID u = 0; u < N; u++)
        {
            if (S[u] != 0)
                continue;
            for (NodeID v = 0; v < N; v++)
            {
                if (S[v] != 0 || v == u || W[u] > W[v] || !g[u][v])
                    continue;

                int dom = 1;
                for (NodeID i = 0; i < N; i++)
                {
                    if (i == u || i == v)
                        continue;
                    if (S[i] == 0 && g[v][i] && !g[u][i])
                    {
                        dom = 0;
                        break;
                    }
                }

                if (dom)
                {
                    imp = 1;
                    S[u] = -1;
                    break;
                }
            }
        }
    }
}

void tiny_solver_solve(tiny_solver *solver, double tl, NodeWeight Wl)
{
    double start = omp_get_wtime();
    int layer = 0, vertex = 0, branch_count = 0;

    solver->time_limit_exceeded = 0;
    solver->weight_limit_exceeded = 0;

    for (int i = 0; i < solver->subgraph_N; i++)
        solver->per_layer_solution[layer][i] = 0;

    solver->per_layer_weight[layer] = 0;
    solver->independent_set_weight = 0;

    solver->branch[layer] = 1;

    tiny_solver_reduce(solver, layer);

    while (layer >= 0)
    {
        if ((branch_count++ & 1023) == 0)
        {
            branch_count = 0;
            double elapsed = omp_get_wtime() - start;
            if (elapsed > tl)
            {
                solver->time_limit_exceeded = 1;
                return;
            }
        }

        // Find max degree vertex
        vertex = solver->subgraph_N;
        NodeID max_degree = 0;
        for (NodeID i = 0; i < solver->subgraph_N; i++)
        {
            if (solver->per_layer_solution[layer][i] != 0)
                continue;
            NodeID degree = 0;
            for (NodeID j = 0; j < solver->subgraph_N; j++)
                if (i != j && solver->subgraph[i][j] &&
                    solver->per_layer_solution[layer][j] == 0)
                    degree++;
            if (degree > max_degree)
            {
                max_degree = degree;
                vertex = i;
            }
        }

        // Fast upperbound
        NodeWeight ub = solver->per_layer_weight[layer];
        for (NodeID i = 0; i < solver->subgraph_N; i++)
            if (solver->per_layer_solution[layer][i] == 0)
                ub += solver->subgraph_W[i];

        // Base case
        if (vertex == solver->subgraph_N || ub <= solver->independent_set_weight)
        {
            // New best solution
            if (solver->per_layer_weight[layer] > solver->independent_set_weight)
            {
                solver->independent_set_weight = solver->per_layer_weight[layer];
                for (NodeID i = 0; i < solver->subgraph_N; i++)
                    solver->independent_set[i] = solver->per_layer_solution[layer][i];

                if (solver->independent_set_weight > Wl)
                {
                    solver->weight_limit_exceeded = 1;
                    return;
                }
            }

            // Backtrack to last branch
            layer--;
            while (layer >= 0 && solver->branch[layer] == 0)
                layer--;

            if (layer >= 0)
                solver->branch[layer] = 0;
        }
        // Include branch
        else if (solver->branch[layer] == 1)
        {
            layer++;
            solver->branch[layer] = 1;
            for (int i = 0; i < solver->subgraph_N; i++)
                solver->per_layer_solution[layer][i] = solver->per_layer_solution[layer - 1][i];

            solver->per_layer_solution[layer][vertex] = 1;
            solver->per_layer_weight[layer] = solver->per_layer_weight[layer - 1] + solver->subgraph_W[vertex];

            for (NodeID u = 0; u < solver->subgraph_N; u++)
                if (u != vertex && solver->subgraph[vertex][u])
                    solver->per_layer_solution[layer][u] = -1;

            tiny_solver_reduce(solver, layer);
        }
        // Exclude branch
        else if (solver->branch[layer] == 0)
        {
            layer++;
            solver->branch[layer] = 1;
            for (NodeID i = 0; i < solver->subgraph_N; i++)
                solver->per_layer_solution[layer][i] = solver->per_layer_solution[layer - 1][i];

            solver->per_layer_solution[layer][vertex] = -1;
            solver->per_layer_weight[layer] = solver->per_layer_weight[layer - 1];

            tiny_solver_reduce(solver, layer);
        }
    }
}
