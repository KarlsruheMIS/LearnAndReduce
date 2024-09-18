
#ifndef TINY_SOLVER_H
#define TINY_SOLVER_H

#include "definitions.h"

#ifdef __cplusplus
extern "C"
{
#endif
    typedef struct
    {
        // Subgraph
        NodeID subgraph_N;
        char **subgraph;
        NodeWeight *subgraph_W;

        // Vertex mapping
        NodeID *forward_map, *reverse_map;

        // Solution
        NodeID *independent_set;
        NodeWeight independent_set_weight;

        // Flags
        int time_limit_exceeded,
            weight_limit_exceeded,
            node_limit_exceeded;

        // Internal structures
        char **per_layer_solution;
        NodeWeight *per_layer_weight;
        int *branch;
    } tiny_solver;

    tiny_solver *tiny_solver_init(NodeID N);

    void tiny_solver_free(tiny_solver *solver);

    void tiny_solver_clear(tiny_solver *solver);
    void tiny_solver_clear_solution(tiny_solver *solver);

    void tiny_solver_solve(tiny_solver *solver, double tl, NodeWeight Wl);

#ifdef __cplusplus
}
#endif

#endif // TINY_SOLVER_H
