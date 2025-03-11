#include "twin_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;
bool twin_reduction::reduce(reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_twin) return false;
    size_t oldn = br_alg->status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
        reduced_nodes += (oldn - br_alg->status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != br_alg->status.remaining_nodes) std::cout << "twin redu -> " << (oldn - br_alg->status.remaining_nodes) << std::endl;
	return oldn != br_alg->status.remaining_nodes;
}
inline bool twin_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_twin) return false;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &twin_candidates_set = br_alg->set_1;
    auto &tmp_set = br_alg->set_2;

    size_t oldn = status.remaining_nodes;
    NodeID twin;
    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    if (try_neighborhood_reduction(v, br_alg, neighbors_weight))
    {
        return oldn != status.remaining_nodes;
    }
    get_neighborhood_vector(v, br_alg, neighbors);

    twin_candidates_set.clear();
    bool candidates_empty = true;

    for (NodeID neighbor : status.graph[neighbors[0]])
    {
        if (is_reduced(neighbor, br_alg))
            continue;
        if (neighbor != v && br_alg->deg(neighbor) == neighbors.size())
        {
            twin_candidates_set.add(neighbor);
            candidates_empty = false;
            twin = neighbor;
        }
    }

    for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++)
    {
        NodeID neighbor = neighbors[i];
        tmp_set.clear();
        candidates_empty = true;

        for (NodeID candidate : status.graph[neighbor])
        {
            if (is_reduced(candidate, br_alg))
                continue;
            if (twin_candidates_set.get(candidate))
            {
                tmp_set.add(candidate);
                candidates_empty = false;
                twin = candidate;
            }
        }

        std::swap(twin_candidates_set, tmp_set);
    }

    if (candidates_empty)
        return false;

    if (status.weights[v] + status.weights[twin] >= neighbors_weight)
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(twin, IS_status::included);
    }
    else
    {
        fold(br_alg, v, twin);
    }

    return oldn != status.remaining_nodes;
}
void twin_reduction::fold(reduce_algorithm *br_alg, NodeID main, NodeID twin)
{
    auto &status = br_alg->status;

    restore_vec.push_back({main, twin});

    br_alg->set(twin, IS_status::folded, true);
    status.weights[main] += status.weights[twin];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(main);
    br_alg->add_next_level_neighborhood(main);
}
void twin_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.twin);
    status.weights[data.main] -= status.weights[data.twin];

    restore_vec.pop_back();
}
void twin_reduction::apply(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto main = restore_vec.back().main;
    auto twin = restore_vec.back().twin;

    restore(br_alg);

    if (status.node_status[main] == IS_status::included)
    {
        status.node_status[twin] = IS_status::included;
    }
    else
    {
        status.node_status[twin] = IS_status::excluded;
    }
}

inline int twin_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &twin_candidates_set = br_alg->set_1;
    auto &tmp_set = br_alg->set_2;

    size_t oldn = status.remaining_nodes;
    NodeID twin;
    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    get_neighborhood_vector(v, br_alg, neighbors);

    twin_candidates_set.clear();
    bool candidates_empty = true;

    for (NodeID neighbor : status.graph[neighbors[0]])
    {
        if (neighbor != v && br_alg->deg(neighbor) == neighbors.size())
        {
            twin_candidates_set.add(neighbor);
            candidates_empty = false;
            twin = neighbor;
        }
    }

    for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++)
    {
        NodeID neighbor = neighbors[i];
        tmp_set.clear();
        candidates_empty = true;

        for (NodeID candidate : status.graph[neighbor])
        {
            if (twin_candidates_set.get(candidate))
            {
                tmp_set.add(candidate);
                candidates_empty = false;
                twin = candidate;
            }
        }

        std::swap(twin_candidates_set, tmp_set);
    }

    if (candidates_empty)
        return false;

    label.push_back(v);
    label.push_back(twin);
    return true;
}