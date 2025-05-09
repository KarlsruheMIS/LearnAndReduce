#include "domination_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool domination_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_domination)
        return false;
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
	// if (oldn != br_alg->status.remaining_nodes) std::cout << "domination redu -> " << (oldn - br_alg->status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes;
}
inline bool domination_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
{
    if (br_alg->config.disable_domination)
        return false;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->set_1;
    size_t oldn = status.remaining_nodes;
#ifndef gen_training_data
    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (try_neighborhood_reduction(v, br_alg, neighbor_weights))
    {
        return oldn != status.remaining_nodes;
    }
#endif

    size_t neighbors_count = 0;
    bool is_subset;

    neighbors.clear();
    neighbors.add(v);
    for (NodeID neighbor : status.graph[v])
    {
        neighbors.add(neighbor);
        neighbors_count++;
    }

    for (NodeID neighbor : status.graph[v])
    {
        if (br_alg->deg(neighbor) > neighbors_count)
            continue;
        if (status.weights[neighbor] < status.weights[v])
            continue;

        is_subset = true;

        for (NodeID neighbor2 : status.graph[neighbor])
        {
            if (!neighbors.get(neighbor2))
            {
                is_subset = false;
                break;
            }
        }

        if (is_subset)
        {
            br_alg->set(v, IS_status::excluded);
            break;
        }
    }

	return oldn != status.remaining_nodes;
}
inline int domination_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)  
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->set_1;
    size_t oldn = status.remaining_nodes;
    size_t neighbors_count = 0;
    bool is_subset;

    neighbors.clear();
    neighbors.add(v);
    for (NodeID neighbor : status.graph[v])
    {
        neighbors.add(neighbor);
        neighbors_count++;
    }

    for (NodeID neighbor : status.graph[v])
    {
        if (br_alg->deg(neighbor) > neighbors_count)
            continue;
        if (status.weights[neighbor] < status.weights[v])
            continue;

        is_subset = true;

        for (NodeID neighbor2 : status.graph[neighbor])
        {
            if (!neighbors.get(neighbor2))
            {
                is_subset = false;
                break;
            }
        }

        if (is_subset)
        {
            label.push_back(v);
            return true;
        }
    }

	return false;
}
