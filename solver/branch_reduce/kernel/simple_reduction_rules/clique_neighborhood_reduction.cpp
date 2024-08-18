#include "clique_neighborhood_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;
bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_clique_neighborhood) return false;
    auto &status = br_alg->status;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool clique_neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_clique_neighborhood) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    if (partition_into_cliques(v, br_alg))
    {
        br_alg->set(v, IS_status::included);
    }

    return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction::partition_into_cliques(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &weights = br_alg->status.weights;
    auto &neighbors_vec = br_alg->buffers[0];
    auto &clique_neighbors_set = br_alg->set_1;
    target_weight = weights[v];

    neighbor_weights = get_neighborhood_weight(v, br_alg);
#ifdef gen_training_data
    if (neighbor_weights <= target_weight)
        return false;
#endif
#ifndef gen_training_data
    if (neighbor_weights <= target_weight)
        return true;
#endif
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    // partition neigbors of v into cliques
    while (neighbors_vec.size() >= 2 && br_alg->config.reduction_time_limit > br_alg->t.elapsed())
    {
        NodeID max_neighbor = *std::max_element(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b)
                                                { return weights[a] < weights[b]; });
        clique_neighbors_set.clear();
        for (auto neighbor : neighbors_vec)
        {
            clique_neighbors_set.add(neighbor);
        }
        clique_neighbors_set.remove(max_neighbor);

        if (expand_clique(max_neighbor, neighbors_vec, clique_neighbors_set, br_alg))
            return true;
    }

    return false;
}
bool clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, std::vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set, branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& temp_set = br_alg->set_2;
	auto& clique_set = br_alg->double_set;
	clique_set.clear();

    size_t local_max = status.n;
    NodeWeight local_max_weight = 0;
    bool intersection_empty = true;

    while (true)
    {
        // intersect neighbors of clique with neighbors of max_neighbor
        intersection_empty = true;
        local_max_weight = 0;
        clique_set.add(max_neighbor);
        temp_set.clear();

        for (auto neighbor : status.graph[max_neighbor])
        {
            if (clique_neighbors_set.get(neighbor))
            {
                temp_set.add(neighbor);
                intersection_empty = false;

                if (status.weights[neighbor] > local_max_weight)
                {
                    local_max = neighbor;
                    local_max_weight = status.weights[neighbor];
                }
            }
        }
        if (intersection_empty)
            break;
        // if (local_max == status.n) break;

        // add local_max to current clique
        neighbor_weights -= local_max_weight;
        if (neighbor_weights <= target_weight)
            return true;

        std::swap(clique_neighbors_set, temp_set);
        clique_neighbors_set.remove(local_max);
        max_neighbor = local_max;
    }

    auto &reamining_neighbors = br_alg->buffers[1];
    reamining_neighbors.clear();

    // adjust neigbors_vec
    for (auto neighbor : neighbors_vec)
    {
        if (!clique_set.get(neighbor))
            reamining_neighbors.push_back(neighbor);
    }

    std::swap(reamining_neighbors, neighbors_vec);
    return false;
}
