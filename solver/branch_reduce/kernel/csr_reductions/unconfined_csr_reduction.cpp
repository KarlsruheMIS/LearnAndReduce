#include "unconfined_csr_reduction.h"
#include "general_reduction.h"
#include "branch_and_reduce_algorithm.h"

#include <limits>

#ifdef __cplusplus
extern "C" {
#endif
    #include "include/graph.h"
    #include "include/reductions.h"
#ifdef __cplusplus
}
#endif

typedef branch_and_reduce_algorithm::IS_status IS_status;

graph unconfined_csr_reduction::build_graph(branch_and_reduce_algorithm *br_alg) {

    auto &status = br_alg->status;
    auto &reverse_map = br_alg->buffers[0]; 
    auto &map = br_alg->buffers[1]; 
    auto &V = br_alg->buffers[2]; 
    auto &E = br_alg->buffers[3];
    auto &W = br_alg->weight_buffer;

    reverse_map.resize(status.remaining_nodes, status.n);
    map.resize(status.n, status.n);

    NodeID N = status.remaining_nodes;
    EdgeID M = 0;
    NodeID u = 0;
    for (NodeID v = 0; v < status.n; v++)
    {
        if (! status.node_status[v] == IS_status::not_set) continue;
        map[v] = u;
        reverse_map[u] = v;
        u++;
        M += status.graph[v].size();
    }
    V.resize(N + 1);
    E.resize(M);
    W.resize(N);

    // build CSR
    std::vector<NodeID> tmp;
    EdgeID ei = 0;
    for (NodeID u : reverse_map)
    {
        assert(map[u] < status.n && "map is not set");
        V[map[u]] = ei;
        W[map[u]] = status.weights[u];

        // get the sorted, unreduced neighborhood
        tmp.clear();
        for (NodeID neighbor : status.graph[u])
        {
            if (status.node_status[neighbor] == IS_status::not_set)
                tmp.push_back(map[neighbor]);
        }

        std::sort(tmp.begin(), tmp.end());
        for (NodeID v : tmp)
            E[ei++] = v;
    }
    V[N] = ei;
    return (graph){N, V.data(), E.data(), W.data()};
}


bool unconfined_csr_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_unconfined) return false;
    if (br_alg->blowing_up) return false;
    auto &status = br_alg->status;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    size_t oldn = status.remaining_nodes;
    graph g = build_graph(br_alg);
    void *R = reduction_init(g.N, g.V[g.N]);
    int *A = malloc(sizeof(int) * g.N);
    for (NodeID i = 0; i < g.N; i++)
        A[i] = 1;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v, g, A, R); });

    free(A);
    reduction_free(R);
#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool unconfined_csr_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v, graph &g, int *A, void *R)
{
    auto &status = br_alg->status;
    auto &map = br_alg->buffers[1]; 
    size_t oldn = status.remaining_nodes;
    assert(status.node_status[v] == IS_status::not_set);
    if (reduction_unconfined_csr(R, g.N, g.V, g.E, g.W, A, map[v]))
    {
        br_alg->set(v, IS_status::excluded);
        A[map[v]] = 0;
        return true;
    }
    else 
        return false;
}

















// inline bool unconfined_csr_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID global_u, graph &g)
// {
//     auto &status = br_alg->status;
//     auto &reverse_map = br_alg->buffers[0]; 
//     auto &map = br_alg->buffers[1]; 
//     size_t oldn = status.remaining_nodes;
//     assert(status.node_status[global_u] == IS_status::not_set);
//     NodeID u = map[global_u];
//     assert(u < status.n && "Node not in graph");
//      if (!A[u])
//         return false;
//     NodeID n = 0, m = 0;
//     S.resize(g.N); 
//     NS.resize(g.N);
//     T.resize(g.N);
//     S_B.assign(g.N, 0);
//     NSI_B.assign(g.N, 0);
//     IS_B.assign(g.N, 0);
//     S[n++] = u;
//     S_B[u] = 1;
//     NSI_B[u] = 1;
//     for (EdgeID i = g.V[u]; i < g.V[u + 1]; i++)
//     {
//         if (!A[g.E[i]])
//             continue;
//         NS[m++] = g.E[i];
//         NSI_B[g.E[i]] = 1;
//     }
//     bool first = true;
//     while (m > 0)
//     {
//         NodeID v = NS[--m];
//         if (S_B[v])
//             continue;
//         NodeWeight sw = 0;
//         if (first)
//         {
//             sw = g.W[u];
//         }
//         else
//         {
//             for (EdgeID i = g.V[v]; i < g.V[v + 1] && sw <= g.W[v]; i++)
//                 if (A[g.E[i]] && S_B[g.E[i]])
//                     sw += g.W[g.E[i]];
//         }
//         if (sw > g.W[v])
//             continue;
//         NodeID p = 0;
//         NodeWeight is = 0;
//         NodeWeight min = std::numeric_limits<NodeWeight>::max();
//         for (EdgeID i = g.V[v]; i < g.V[v + 1] && sw + is - min <= g.W[v]; i++)
//         {
//             NodeID w = g.E[i];
//             if (A[w] && !NSI_B[w])
//             {
//                 T[p++] = w;
//                 is += g.W[w];
//                 IS_B[w] = 1;
//                 if (g.W[w] < min)
//                     min = g.W[w];
//             }
//         }
//         bool ind = sw + is - min <= g.W[v];
//         for (NodeID i = 0; i < p && ind; i++)
//         {
//             NodeID w = T[i];
//             for (EdgeID j = g.V[w]; j < g.V[w + 1] && ind; j++)
//                 if (A[g.E[j]] && IS_B[g.E[j]])
//                     ind = 0;
//         }
//         for (NodeID i = 0; i < p; i++)
//             IS_B[T[i]] = 0;
//         if (ind && sw + is <= g.W[v]) // Can reduce u
//         {
//             br_alg->set(global_u, IS_status::excluded, true);
//             A[u] = false;   
//             return true;
//         }
//         else if (ind && sw + is > g.W[v] && sw + is - min <= g.W[v]) // Extend S
//         {
//             first = false;
//             for (NodeID i = 0; i < p; i++)
//             {
//                 NodeID w = T[i];
//                 S[n++] = w;
//                 S_B[w] = 1;
//                 NSI_B[w] = 1;
//                 for (EdgeID j = g.V[w]; j < g.V[w + 1]; j++)
//                 {
//                     if (!A[g.E[j]])
//                         continue;
//                     NS[m++] = g.E[j];
//                     NSI_B[g.E[j]] = 1;
//                 }
//             }
//         }
//     }
//     for (NodeID i = 0; i < n; i++)
//     {
//         NodeID v = S[i];
//         for (EdgeID j = g.V[v]; j < g.V[v + 1]; j++)
//             NSI_B[g.E[j]] = 0;
//         S_B[v] = 0;
//         NSI_B[v] = 0;
//     }
//     return false;
// }
