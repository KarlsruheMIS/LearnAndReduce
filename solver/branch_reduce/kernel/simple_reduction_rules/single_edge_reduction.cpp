#include "single_edge_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;
bool single_edge_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_basic_se) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_basic_se) return false;
    if (br_alg->deg(v) == 0)
        return false;

    auto &status = br_alg->status;
    auto &graph = br_alg->status.graph;
    auto &weights = br_alg->status.weights;
    auto &neighbors = br_alg->set_1;
    size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, neighbors);

    for (NodeID neighbor : graph[v])
    {
        if (is_reduced(neighbor, br_alg))
            continue;
        if (weights[v] <= weights[neighbor])
        { // otherwise not applicable to this edge
            NodeWeight partial_neighbor_sum = 0;
            for (NodeID second_neighbor : status.graph[neighbor])
            {
                if (!is_reduced(second_neighbor, br_alg) && !neighbors.get(second_neighbor))
                {
                    partial_neighbor_sum += status.weights[second_neighbor];
                    if (partial_neighbor_sum > weights[neighbor])
                        break;
                }
            }

// note: weight of v is in partial_neighbor_sum included
// if N(neighbor) \subset N(v) partial_neighbor_sum = weights[v]
#ifdef gen_training_data
            if (partial_neighbor_sum == weights[v])
                break; // should be domination
#endif
            if (partial_neighbor_sum <= weights[neighbor])
            {
                br_alg->set(v, IS_status::excluded);
                break;
            }
        }
    }

    return oldn != status.remaining_nodes;
}
