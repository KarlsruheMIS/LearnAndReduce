#include "clique_neighborhood_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"
#include "clique_neighborhood_fast_reduction.h"

typedef reduce_algorithm::IS_status IS_status;

bool clique_neighborhood_reduction_fast::reduce(reduce_algorithm *r_alg)
{
    // if (r_alg->config.disable_clique_neighborhood_fast) return false;
    auto &status = r_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(r_alg, [&](NodeID v)
                            { reduce_vertex(r_alg, v); });

	// if (oldn != status.remaining_nodes) std::cout << "clique_neighborhood fast redu -> " << (oldn - status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes;
}
inline bool clique_neighborhood_reduction_fast::reduce_vertex(reduce_algorithm *r_alg, NodeID v)
{
// if (r_alg->config.disable_clique_neighborhood_fast) return false;
#ifdef REDUCTION_INFO
    r_alg->reduction_timer.restart();
#endif
    auto &status = r_alg->status;
    size_t oldn = status.remaining_nodes;
    auto &neighbors = r_alg->buffers[0];
    auto &neighborhood = r_alg->set_1;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, r_alg);
    if (this->try_neighborhood_reduction(v, r_alg, neighbor_weights))
    {
        return oldn != r_alg->status.remaining_nodes;
    }
    this->get_neighborhood_vector(v, r_alg, neighbors);
    this->get_neighborhood_set(v, r_alg, neighborhood);

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
        r_alg->set(v, IS_status::included);
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += r_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline int clique_neighborhood_reduction_fast::generate_data(reduce_algorithm *r_alg, NodeID v, std::vector<NodeID>& label)
{
    auto &status = r_alg->status;
    size_t oldn = status.remaining_nodes;
    auto &neighbors = r_alg->buffers[0];
    auto &neighborhood = r_alg->set_1;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, r_alg);
    this->get_neighborhood_vector(v, r_alg, neighbors);
    this->get_neighborhood_set(v, r_alg, neighborhood);

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
        label.push_back(v);
    }
    return is_reducible;
}
