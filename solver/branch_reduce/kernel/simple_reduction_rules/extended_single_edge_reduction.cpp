#include "extended_single_edge_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;
bool extended_single_edge_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_extended_se) return false;
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
inline bool extended_single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_extended_se) return false;
    if (br_alg->deg(v) == 0)
        return false;

	auto& status = br_alg->status;
    auto& weights = status.weights;
	size_t oldn = status.remaining_nodes;
	auto& neighbors = br_alg->set_1;
	auto& neighbors_vec = br_alg->buffers[0];


    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    NodeID max_weight_neighbor = get_max_weight_neighbor(v, br_alg);
    assert(max_weight_neighbor != v && "ERROR: max_weight_neighbor == v");
    if (weights[v] < neighbors_weight - weights[max_weight_neighbor])
        return false;
    if (try_neighborhood_reduction(v, br_alg, neighbors_weight))
        return oldn != br_alg->status.remaining_nodes;
    
    get_neighborhood_set(v, br_alg, neighbors);
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    std::sort(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b)
              { return weights[a] > weights[b]; });
    for (NodeID max_neighbor : neighbors_vec)
    {
        if (v > max_neighbor)
            continue;
        if (weights[v] < neighbors_weight - weights[max_neighbor])
            return oldn != status.remaining_nodes;

        bool progress = false;
        for (NodeID neighbor : status.graph[max_neighbor])
        {
            if (neighbor == v)
                continue;
            if (is_reduced(neighbor, br_alg))
                continue;
            // exclude neighborhood intersection and update neighborhood
            if (neighbors.get(neighbor))
            {
                br_alg->set(neighbor, IS_status::excluded);
                progress = true;
            }
        }

        // if (progress)
        //     break;
    }

    return oldn != status.remaining_nodes;
}

inline int extended_single_edge_reduction::generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    if (br_alg->deg(v) == 0)
        return 0;
	auto& status = br_alg->status;
    auto& weights = status.weights;
	size_t oldn = status.remaining_nodes;
	auto& neighbors = br_alg->set_1;
	auto& neighbors_vec = br_alg->buffers[0];


    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    NodeID max_weight_neighbor = get_max_weight_neighbor(v, br_alg);
    assert(max_weight_neighbor != v && "ERROR: max_weight_neighbor == v");
    assert(label.size() == 0 && "ERROR: labeld_vertices not empty");
    if (weights[v] < neighbors_weight - weights[max_weight_neighbor])
        return 0;
    
    get_neighborhood_set(v, br_alg, neighbors);
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    std::sort(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b)
              { return weights[a] > weights[b]; });
    for (NodeID max_neighbor : neighbors_vec)
    {
        if (weights[v] < neighbors_weight - weights[max_neighbor])
            return 0;

        for (NodeID neighbor : status.graph[max_neighbor])
        {
            if (neighbor == v)
                continue;
            // exclude neighborhood intersection and update neighborhood
            if (neighbors.get(neighbor))
                label.push_back(neighbor);
        }

    }

    return label.size() > 0;
}
