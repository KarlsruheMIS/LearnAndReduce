#include "unconfined_csr_reduction.h"
#include "general_reduction.h"
#include "reduce_algorithm.h"

#include <limits>

#ifdef __cplusplus
extern "C"
{
#endif
#include "csr_graph.h"
#include "include/reductions.h"
#ifdef __cplusplus
}
#endif

typedef reduce_algorithm::IS_status IS_status;

bool unconfined_csr_reduction::reduce(reduce_algorithm *r_alg)
{
    if (r_alg->config.disable_unconfined)
        return false;
    if (r_alg->blowing_up)
        return false;
    auto &status = r_alg->status;
#ifdef REDUCTION_INFO
    r_alg->reduction_timer.restart();
#endif
    size_t oldn = status.remaining_nodes;
    // map and reverse map are stored in
    // reverse_map = buffers[0];
    // map = buffers[1];
    csr_graph g = r_alg->build_global_graph_csr();
    void *R = reduction_init(g.N, g.V[g.N]);
    int *A = static_cast<int *>(malloc(sizeof(int) * g.N));
    for (NodeID i = 0; i < g.N; i++)
        A[i] = 1;

    for_each_changed_vertex(r_alg, [&](NodeID v)
                            { reduce_vertex(r_alg, v, g, A, R); });

    gnn_filter_marker(r_alg->config, marker.current.back());
    free(A);
    reduction_free(R);
#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += r_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool unconfined_csr_reduction::reduce_vertex(reduce_algorithm *r_alg, NodeID v, csr_graph &g, int *A, void *R)
{
    auto &status = r_alg->status;
    auto &map = r_alg->buffers[1];
    size_t oldn = status.remaining_nodes;
    assert(status.node_status[v] == IS_status::not_set);
    if (reduction_unconfined_csr(R, g.N, g.V, g.E, g.W, A, map[v]))
    {
        r_alg->set(v, IS_status::excluded);
        A[map[v]] = 0;
        return true;
    }
    else
        return false;
}

void unconfined_csr_reduction::generate_global_data(reduce_algorithm *r_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index)
{
    auto &status = r_alg->status;
    auto &map = r_alg->buffers[1];
    size_t oldn = status.remaining_nodes;
    csr_graph g = r_alg->build_global_graph_csr();
    void *R = reduction_init(g.N, g.V[g.N]);
    int *A = static_cast<int *>(malloc(sizeof(int) * g.N));
    for (NodeID i = 0; i < g.N; i++)
        A[i] = 1;

    for (NodeID v = 0; v < g.N; v++)
    {
        if (reduction_unconfined_csr(R, g.N, g.V, g.E, g.W, A, map[v]))
            reduction_data[reduction_index][v] = 1;
    }
    free(A);
    reduction_free(R);
}