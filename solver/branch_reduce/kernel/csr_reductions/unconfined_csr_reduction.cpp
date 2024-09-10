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
    tmp.reserve(N);
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
    int *A = static_cast<int*>(malloc(sizeof(int) * g.N));
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

bool unconfined_csr_reduction::generate_global_data(branch_and_reduce_algorithm *br_alg, std::vector<NodeID> &label)
{
    auto &status = br_alg->status;
    auto &map = br_alg->buffers[1]; 
    size_t oldn = status.remaining_nodes;
    graph g = build_graph(br_alg);
    void *R = reduction_init(g.N, g.V[g.N]);
    int *A = static_cast<int*>(malloc(sizeof(int) * g.N));
    for (NodeID i = 0; i < g.N; i++)
        A[i] = 1;

    for (NodeID v = 0; v < g.N; v++)
    {
        if (reduction_unconfined_csr(R, g.N, g.V, g.E, g.W, A, map[v]))
            label.push_back(v);
    }

    free(A);
    reduction_free(R);
    return label.size() > 0;
}

