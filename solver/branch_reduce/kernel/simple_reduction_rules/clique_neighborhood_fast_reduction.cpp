#include "clique_neighborhood_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_clique_neighborhood_fast) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

	// if (oldn != status.remaining_nodes) std::cout << "clique_neighborhood fast redu -> " << (oldn - status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes;
}
inline bool clique_neighborhood_reduction_fast::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
// if (br_alg->config.disable_clique_neighborhood_fast) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    auto &neighbors = br_alg->buffers[0];
    auto &neighborhood = br_alg->set_1;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (this->try_neighborhood_reduction(v, br_alg, neighbor_weights))
    {
        return oldn != br_alg->status.remaining_nodes;
    }
    this->get_neighborhood_vector(v, br_alg, neighbors);
    this->get_neighborhood_set(v, br_alg, neighborhood);

    std::sort(neighbors.begin(), neighbors.end(), [&status](const NodeID lhs, const NodeID rhs)
              { return status.weights[lhs] > status.weights[rhs]; });

    bool is_reducible = false;

    for (size_t i = 0; i < neighbors.size() && !is_reducible; i++)
    {
        NodeID neighbor1 = neighbors[i];

        if (!neighborhood.get(neighbor1))
            continue;

        for (NodeID neighbor2 : status.graph[neighbor1])
        {
            if (neighbor2 != neighbor1 && neighborhood.get(neighbor2))
            {
                // triangle [v, neighbor1, neighbor2] found
                neighbor_weights -= std::min(status.weights[neighbor1], status.weights[neighbor2]);
                neighborhood.remove(neighbor1);
                neighborhood.remove(neighbor2);

                if (status.weights[v] >= neighbor_weights)
                {
                    is_reducible = true;
                }

                break;
            }
        }
    }

    if (is_reducible)
    {
        br_alg->set(v, IS_status::included);
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
