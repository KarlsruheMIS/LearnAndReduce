#include "fold1_reduction.h"
#include "general_reduction.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;
bool fold1_reduction::reduce(reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_fold1) return false;
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
inline bool fold1_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_fold1) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    assert(status.node_status[v] == IS_status::not_set);

    if (status.weights[v] == 0)
    {
        br_alg->set(v, IS_status::excluded);
    }
    else if (br_alg->deg(v) == 0)
    {
        br_alg->set(v, IS_status::included);
    }
    else if (br_alg->deg(v) == 1)
    {
        NodeID neighbor = status.graph[v][0];
        if (status.weights[neighbor] > status.weights[v])
            fold(br_alg, {v, neighbor});
        else
            br_alg->set(v, IS_status::included, true);
    }

    return oldn != status.remaining_nodes;
}
void fold1_reduction::fold(reduce_algorithm *br_alg, const fold_nodes &nodes)
{

    auto &status = br_alg->status;

    restore_vec.push_back({nodes, status.weights[nodes.deg1_node]});
    br_alg->set(nodes.deg1_node, IS_status::folded);

    status.reduction_offset += status.weights[nodes.deg1_node];
    status.weights[nodes.fold_node] -= status.weights[nodes.deg1_node];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.fold_node);
    br_alg->add_next_level_neighborhood(nodes.fold_node);
}
void fold1_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.nodes.deg1_node);

    status.weights[data.nodes.fold_node] += data.deg1_weight;
    status.reduction_offset -= data.deg1_weight;

    restore_vec.pop_back();
}
void fold1_reduction::apply(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto nodes = restore_vec.back().nodes;
    auto main_status = status.node_status[nodes.fold_node];
    restore(br_alg);

    if (main_status == IS_status::included)
    {
        status.node_status[nodes.fold_node] = IS_status::included;
        status.node_status[nodes.deg1_node] = IS_status::excluded;

	} 
    else if (main_status == IS_status::excluded) {
		status.node_status[nodes.fold_node] = IS_status::excluded;
		status.node_status[nodes.deg1_node] = IS_status::included;
	}
	status.is_weight += status.weights[nodes.deg1_node];
}
bool fold1_reduction::is_suited(NodeID v, reduce_algorithm *br_alg)
{
    return br_alg->deg(v) <= 1;
}
inline int fold1_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    // if (br_alg->config.disable_fold1) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    assert(status.node_status[v] == IS_status::not_set);

    if (status.weights[v] == 0 || br_alg->deg(v) == 0 || br_alg->deg(v) == 1)
    {
        label.push_back(v);
    }
    return label.size() > 0;
}