#include "extended_twin_reduction.h"

#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool extended_twin_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_extended_twin)
        return false;
    if (br_alg->blowing_up)
        return false;
    auto &status = br_alg->status;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    size_t oldn = status.remaining_nodes;
    bool progress = false;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { progress = reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return progress;
}
inline bool extended_twin_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &candidate_neighbors = br_alg->buffers[1];
    auto &neighbors_set = br_alg->set_2;

    size_t oldn = status.remaining_nodes;
    get_neighborhood_vector(v, br_alg, neighbors);
    get_neighborhood_set(v, br_alg, neighbors_set);

    std::sort(neighbors.begin(), neighbors.end(), [&](NodeID a, NodeID b)
              { return br_alg->deg(a) < br_alg->deg(b); });

    for (NodeID cand : status.graph[neighbors[0]])
    {
        if (is_reduced(cand, br_alg))
            continue;
        if (cand != v && br_alg->deg(cand) >= neighbors.size())
        {
            NodeID common_node_count = 0;
            NodeWeight cand_neighbors_weight = get_neighborhood_weight(cand, br_alg);
            if (status.weights[v] + status.weights[cand] < cand_neighbors_weight)
                continue;
            // check if N(v) \subset N(cand)
            for (NodeID n : status.graph[cand])
            {
                if (!neighbors_set.get(n))
                    continue;
                common_node_count++;
                if (common_node_count == neighbors.size())
                {
                    if (br_alg->deg(v) == br_alg->deg(cand)) // in this case v and cand are actual twins
                        br_alg->set(cand, IS_status::included);
                    br_alg->add_next_level_node(cand);
                    br_alg->set(v, IS_status::included);
                    return true;
                }
            }
        }
    }

    return false;
}

int extended_twin_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label)
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &candidate_neighbors = br_alg->buffers[1];
    auto &neighbors_set = br_alg->set_2;

    size_t oldn = status.remaining_nodes;
    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    get_neighborhood_vector(v, br_alg, neighbors);
    get_neighborhood_set(v, br_alg, neighbors_set);

std:
    sort(neighbors.begin(), neighbors.end(), [&](NodeID a, NodeID b)
         { return br_alg->deg(a) < br_alg->deg(b); });

    for (NodeID cand : status.graph[neighbors[0]])
    {
        if (cand != v && br_alg->deg(cand) >= neighbors.size())
        {
            NodeID common_node_count = 0;
            for (NodeID n : status.graph[cand])
            {
                if (neighbors_set.get(n))
                    common_node_count++;
            }
            // if N[cand] is a superset of N[v]
            if (common_node_count == neighbors.size())
            {
                return true;
            }
        }
    }

    return false;
}