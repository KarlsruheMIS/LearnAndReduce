#include "include/reductions.h"
#include "definitions.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

typedef struct
{
    NodeID *S, *NS, *T;
    int *S_B, *NSI_B, *IS_B;
} reduction_data;

void *reduction_init(NodeID N, EdgeID M)
{
    reduction_data *rp = malloc(sizeof(reduction_data));

    rp->S = malloc(sizeof(NodeID) * N);
    rp->NS = malloc(sizeof(NodeID) * N);
    rp->T = malloc(sizeof(NodeID) * N);
    rp->S_B = malloc(sizeof(int) * N);
    rp->NSI_B = malloc(sizeof(int) * N);
    rp->IS_B = malloc(sizeof(int) * N);

    for (NodeID i = 0; i < N; i++)
    {
        rp->S_B[i] = 0;
        rp->NSI_B[i] = 0;
        rp->IS_B[i] = 0;
    }

    return rp;
}

void reduction_free(void *R)
{
    reduction_data *rp = (reduction_data *)R;

    free(rp->S);
    free(rp->NS);
    free(rp->T);
    free(rp->S_B);
    free(rp->NSI_B);
    free(rp->IS_B);
    free(rp);
}

int reduction_neighborhood_csr(void *R, NodeID N, const EdgeID *V, const NodeID *E,
                               const NodeWeight *W, const int *A, NodeID u)
{
    if (!A[u])
        return 0;

    NodeWeight nw = 0;
    for (EdgeID i = V[u]; i < V[u + 1]; i++)
        if (A[E[i]])
            nw += W[E[i]];

    return nw <= W[u];
}

int reduction_unconfined_csr(void *R, NodeID N, const EdgeID *V, const NodeID *E,
                             const NodeWeight *W, const int *A, NodeID u)
{
    if (!A[u])
        return 0;

    reduction_data *rp = (reduction_data *)R;

    NodeID n = 0, m = 0;

    rp->S[n++] = u;
    rp->S_B[u] = 1;
    rp->NSI_B[u] = 1;
    for (EdgeID i = V[u]; i < V[u + 1]; i++)
    {
        if (!A[E[i]])
            continue;
        rp->NS[m++] = E[i];
        rp->NSI_B[E[i]] = 1;
    }

    int res = 0, first = 1;
    while (m > 0)
    {
        NodeID v = rp->NS[--m];
        if (rp->S_B[v])
            continue;

        NodeWeight sw = 0;
        if (first)
        {
            sw = W[u];
        }
        else
        {
            for (EdgeID i = V[v]; i < V[v + 1] && sw <= W[v]; i++)
                if (A[E[i]] && rp->S_B[E[i]])
                    sw += W[E[i]];
        }

        if (sw > W[v])
            continue;

        NodeID p = 0;
        NodeWeight is = 0, min = LLONG_MAX;

        for (EdgeID i = V[v]; i < V[v + 1] && sw + is <= W[v] + min; i++)
        {
            NodeID w = E[i];
            if (A[w] && !rp->NSI_B[w])
            {
                rp->T[p++] = w;
                is += W[w];
                rp->IS_B[w] = 1;
                if (W[w] < min)
                    min = W[w];
            }
        }

        int ind = sw + is <= W[v] + min;
        for (NodeID i = 0; i < p && ind; i++)
        {
            NodeID w = rp->T[i];
            for (EdgeID j = V[w]; j < V[w + 1] && ind; j++)
                if (A[E[j]] && rp->IS_B[E[j]])
                    ind = 0;
        }

        for (NodeID i = 0; i < p; i++)
            rp->IS_B[rp->T[i]] = 0;

        if (ind && sw + is <= W[v]) // Can reduce u
        {
            res = 1;
            break;
        }
        else if (ind && sw + is > W[v] && sw + is  <= W[v] + min) // Extend S
        {
            first = 0;

            for (NodeID i = 0; i < p; i++)
            {
                NodeID w = rp->T[i];

                rp->S[n++] = w;
                rp->S_B[w] = 1;
                rp->NSI_B[w] = 1;
                for (EdgeID j = V[w]; j < V[w + 1]; j++)
                {
                    if (!A[E[j]])
                        continue;
                    rp->NS[m++] = E[j];
                    rp->NSI_B[E[j]] = 1;
                }
            }
        }
    }

    for (NodeID i = 0; i < n; i++)
    {
        NodeID v = rp->S[i];
        for (EdgeID j = V[v]; j < V[v + 1]; j++)
            rp->NSI_B[E[j]] = 0;
        rp->S_B[v] = 0;
        rp->NSI_B[v] = 0;
    }

    return res;
}