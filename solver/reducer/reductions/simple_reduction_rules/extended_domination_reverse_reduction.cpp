#include "extended_domination_reverse_reduction.h"

#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool extended_domination_reverse_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_extended_domination_reverse)
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
inline bool extended_domination_reverse_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
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
            // check if N(v) \subset N(cand)
            for (NodeID n : status.graph[cand])
            {
                if (!neighbors_set.get(n))
                    continue;
                common_node_count++;
                if (common_node_count == neighbors.size())
                {
                    if (status.weights[v] + status.weights[cand] >= cand_neighbors_weight)
                    {
                        br_alg->set(v, IS_status::included);
                        br_alg->add_next_level_node(cand);
                        // break;
                    }
                    else
                    {
                        fold(br_alg, v, cand);
                    }
                    return true;
                }
            }
        }
    }
    br_alg->config.disable_extended_domination_reverse = true;
    return false;
}
void extended_domination_reverse_reduction::fold(reduce_algorithm *br_alg, NodeID v, NodeID cand)
{

    auto &status = br_alg->status;

    status.modified_stack.push_back(br_alg->MODIFIED_TOKEN);
    restore_vec.push_back({v, cand, status.weights[cand]});

    status.weights[cand] += status.weights[v];
    status.graph.add_edge_undirected(v, cand);

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(v);
    br_alg->add_next_level_neighborhood(v);
    br_alg->add_next_level_neighborhood(cand);
}
void extended_domination_reverse_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    status.weights[data.cand] = data.candWeight;
    status.graph.hide_edge_undirected(data.v, data.cand);
    restore_vec.pop_back();
}
void extended_domination_reverse_reduction::apply(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto v = restore_vec.back().v;
    auto cand = restore_vec.back().cand;
    auto cand_status = status.node_status[cand];
    auto v_status = status.node_status[v];
    assert(status.graph.adjacent(v, cand) && "edge not present");
    assert(!(v_status == IS_status::included && cand_status == IS_status::included) && "neighbor is included");
    restore(br_alg);
    if (cand_status == IS_status::included)
        status.node_status[v] = IS_status::included;
}
int extended_domination_reverse_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label)
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