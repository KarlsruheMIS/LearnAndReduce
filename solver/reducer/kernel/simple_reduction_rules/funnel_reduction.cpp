#include "funnel_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "fast_set.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool funnel_reduction::reduce(reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_funnel) return false;
    size_t oldn = br_alg->status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&](NodeID node)
                            { reduce_vertex(br_alg, node); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != br_alg->status.remaining_nodes;
}
inline bool funnel_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_funnel) return false;
    if (br_alg->deg(v) <= 2) return false; // fold1 or 2
	auto& weights = br_alg->status.weights;
	auto& remaining_n = br_alg->status.remaining_nodes;
	auto& funnel_set = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	size_t oldn = remaining_n;
    NodeID funnel_neighbor = br_alg->status.n;

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    // need one vertex of weight >= v
    if (std::any_of(neighbors.begin(), neighbors.end(), [&](NodeID neighbor)
                    { return weights[neighbor] >= weights[v]; }))
    {
        get_neighborhood_set(v, br_alg, funnel_set);
        funnel_set.add(v);

        if (is_funnel(v, funnel_neighbor, br_alg, funnel_set, neighbors))
        {
            fold({v, funnel_neighbor}, funnel_set, br_alg);
        }
    }

    return oldn != remaining_n;
}
bool funnel_reduction::is_funnel(NodeID node, NodeID& funnel_neighbor, reduce_algorithm* br_alg, fast_set& funnel_set, std::vector<NodeID>& funnel_nodes) 
{
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    if (is_clique(br_alg, funnel_set, funnel_nodes))
        return false;
    funnel_neighbor = status.n;
    for (NodeID v : funnel_nodes)
    {
        assert(v != node && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(v, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");
        if (weights[v] >= weights[node])
        {

            bool skip = false;
            for (auto neighbor : graph[node])
            {
                if (is_reduced(neighbor, br_alg))
                    continue;
                if (neighbor == v)
                    continue;
                if (weights[neighbor] > weights[node])
                {
                    skip = true;
                    break;
                }
            }
            if (skip)
                continue;

            // remove v from funnel_nodes
            auto v_iter = std::find(funnel_nodes.begin(), funnel_nodes.end(), v);
            assert(v_iter != funnel_nodes.end());
            std::swap(*v_iter,funnel_nodes.back());
            funnel_nodes.pop_back();
            funnel_set.remove(v);
            if (is_clique(br_alg, funnel_set, funnel_nodes))
            {
                funnel_neighbor = v;
                funnel_nodes.push_back(funnel_neighbor);
                funnel_set.add(funnel_neighbor);
                return true;
            }
            else
            {
                funnel_nodes.push_back(v);
                funnel_set.add(v);
            }
        }
    }
    return false;
}
bool funnel_reduction::is_clique(reduce_algorithm* br_alg, fast_set& clique_set, std::vector<NodeID>& clique_nodes) {
    auto& status = br_alg->status;
    auto& graph = status.graph;

    for (auto neighbor : clique_nodes)
    {
        size_t count = 0;

        for (NodeID neighbor_2nd : graph[neighbor])
        {
            if (status.node_status[neighbor_2nd] != IS_status::not_set)
                continue;
            if (clique_set.get(neighbor_2nd))
                count++;
        }

        if (count != clique_nodes.size())
            return false; // not all neighbors are in the clique
    }
    return true;
}
void funnel_reduction::fold(const fold_data &data, fast_set &funnel_set, reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &weights = status.weights;
    auto &graph = status.graph;
    auto &f_neighbors = br_alg->set_2;
    NodeID node = data.node;
    NodeID funnel_neighbor = data.funnel_neighbor;

    f_neighbors.clear();
    std::vector<NodeID> outside_funnel_neighbors;
    outside_funnel_neighbors.reserve(status.graph[funnel_neighbor].size());
    std::vector<NodeID> remaining_neighbors;
    remaining_neighbors.reserve(graph[node].size());

    // collect neighbors of funnel neighbor and those outside the funnel set
    for (auto neighbor : graph[funnel_neighbor])
    {
        if (!is_reduced(neighbor, br_alg))
        {
            f_neighbors.add(neighbor);
            if (!funnel_set.get(neighbor))
                outside_funnel_neighbors.push_back(neighbor);
        }
    }

    // common neighbors are not part of the IS in either case
    for (auto neighbor : graph[node])
    {
        if (is_reduced(neighbor, br_alg))
            continue;
        if (neighbor == funnel_neighbor)
            continue;
        if (f_neighbors.get(neighbor))
        {
            br_alg->set(neighbor, IS_status::excluded);
        }
        else
        {
            remaining_neighbors.push_back(neighbor);
        }
    }

    if (remaining_neighbors.size() == 0)
    {
        assert(br_alg->deg(node) == 1 && "ERROR: funnel_reduction::fold: remaining neighbors must be non-empty");
        if (weights[node] >= weights[funnel_neighbor])
        {
            br_alg->set(node, IS_status::included);
        }
        return;
    }

    // remaining_neighbors.resize(remaining_neighbors.size());
    status.folded_stack.push_back(get_reduction_type());
    status.reduction_offset += weights[node];

    // add additional edges connecting the remaining neighbors with funnel neighors outside the funnel set
    restore_vec.push_back({data.node, data.funnel_neighbor, remaining_neighbors, {}});
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++)
    {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];
        assert(!is_reduced(neighbor, br_alg) && "ERROR: funnel_reduction::fold: remaining neighbors must be unset");
        assert(node != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        assert(funnel_neighbor != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        for (auto neighbor_2nd : outside_funnel_neighbors)
        {
            if (!graph.adjacent(neighbor, neighbor_2nd))
            {
                graph.add_edge_undirected(neighbor, neighbor_2nd);
                restore_vec.back().node_vecs[idx].push_back(neighbor_2nd);
            }
        }
        assert(weights[neighbor] + weights[funnel_neighbor] >= weights[node] && "ERROR: funnel_reduction::fold: weight must be larger than node weight");
        br_alg->add_next_level_node(neighbor);
    }
    weights[funnel_neighbor] -= weights[node];
    br_alg->set(node, IS_status::folded, true);
    br_alg->add_next_level_neighborhood(node);
}
void funnel_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.node);
    status.weights[data.funnel_neighbor] += status.weights[data.node];

    // remove added edges:
    for (size_t idx = 0; idx < data.remaining_neighbors.size(); idx++)
    {
        NodeID neighbor = data.remaining_neighbors[idx];
        for (auto neighbor_2nd : restore_vec.back().node_vecs[idx])
        {
            status.graph.hide_edge_undirected(neighbor, neighbor_2nd);
        }
    }
    status.reduction_offset -= status.weights[data.node];

    restore_vec.pop_back();
}
void funnel_reduction::apply(reduce_algorithm *br_alg)
{
    auto &node_status = br_alg->status.node_status;
    auto &is_weight = br_alg->status.is_weight;
    auto &weights = br_alg->status.weights;
    auto data = restore_vec.back();
    auto &remaining = data.remaining_neighbors;

    bool include_node = node_status[data.funnel_neighbor] == IS_status::excluded;
    if (include_node)
        include_node = std::none_of(remaining.begin(), remaining.end(), [&](NodeID neighbor)
                                    { return node_status[neighbor] == IS_status::included; });

    restore(br_alg);

    if (include_node)
    {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
    } else {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
    }
    is_weight += weights[data.node];
}
inline int funnel_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    if (br_alg->deg(v) <= 2) return false; // fold1 or 2
	auto& weights = br_alg->status.weights;
	auto& funnel_set = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
    NodeID funnel_neighbor = br_alg->status.n;

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    // need one vertex of weight >= v
    get_neighborhood_set(v, br_alg, funnel_set);
    funnel_set.add(v);

    if (is_funnel(v, funnel_neighbor, br_alg, funnel_set, neighbors))
    {
        label.push_back(v);
        return true;
    }

    return false;
}