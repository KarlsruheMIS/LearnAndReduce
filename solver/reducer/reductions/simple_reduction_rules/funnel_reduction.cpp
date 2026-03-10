#include "funnel_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "fast_set.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool funnel_reduction::reduce(reduce_algorithm *br_alg)
{
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
    if (br_alg->deg(v) <= 2)
        return false; // fold1 or 2
    auto &weights = br_alg->status.weights;

    NodeID u = br_alg->status.n;
    int higher_weight_neighbor_count = 0;
    for (NodeID neighbor : br_alg->status.graph[v])
    {
        if (weights[neighbor] > weights[v])
        {
            if (higher_weight_neighbor_count > 0)
                return false;
            higher_weight_neighbor_count++;
            u = neighbor;
        }
    }

    if (is_funnel(v, u, br_alg))
    {
        fold({v, u}, br_alg);
        return true;
    }

    return false;
}

bool funnel_reduction::is_funnel(NodeID v, NodeID &u, reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &weights = status.weights;
    auto &funnel_set = br_alg->set_1;
    auto &funnel_nodes = br_alg->buffers[0];

    funnel_set.clear();
    get_neighborhood_vector(v, br_alg, funnel_nodes);
    get_neighborhood_set(v, br_alg, funnel_set);
    funnel_set.add(v);

    if (is_clique(br_alg, funnel_set, funnel_nodes))
        return false; // TODO reduce this here

    // if weight[u] > weight[v], there is only one possible funnel:
    if (u < status.n)
    {
        auto iter = std::find(funnel_nodes.begin(), funnel_nodes.end(), u);
        assert(iter != funnel_nodes.end());
        std::swap(*iter, funnel_nodes.back());
        funnel_nodes.pop_back();
        funnel_set.remove(u);
        if (is_clique(br_alg, funnel_set, funnel_nodes))
            return true;
        return false;
    }

    // Here, v has max weight and we dont know u yet.
    for (NodeID x : funnel_nodes)
    {
        assert(x != v && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(x, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");

        // remove x from funnel_nodes
        auto iter = std::find(funnel_nodes.begin(), funnel_nodes.end(), x);
        assert(iter != funnel_nodes.end());
        std::swap(*iter, funnel_nodes.back());
        funnel_nodes.pop_back();
        funnel_set.remove(x);

        if (is_clique(br_alg, funnel_set, funnel_nodes))
        {
            u = x;
            return true;
        }
        funnel_nodes.push_back(x);
        funnel_set.add(x);
    }
    return false;
}

bool funnel_reduction::is_clique(reduce_algorithm *br_alg, fast_set &clique_set, std::vector<NodeID> &clique_nodes)
{
    auto &status = br_alg->status;
    auto &graph = status.graph;

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

void funnel_reduction::fold(const fold_data &data, reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &weights = status.weights;
    auto &graph = status.graph;
    auto &funnel_set = br_alg->set_1;
    auto &u_neighbors = br_alg->set_2;
    NodeID node = data.node;
    NodeID u = data.funnel_neighbor;
    auto &outside_funnel_neighbors = br_alg->buffers[0];
    auto &remaining_neighbors = br_alg->buffers[1];

    funnel_set.add(data.funnel_neighbor);
    u_neighbors.clear();
    outside_funnel_neighbors.clear();
    remaining_neighbors.clear();

    // collect neighbors of u and those outside the funnel set
    for (auto neighbor : graph[u])
    {
        if (!is_reduced(neighbor, br_alg))
        {
            u_neighbors.add(neighbor);
            if (!funnel_set.get(neighbor))
                outside_funnel_neighbors.push_back(neighbor);
        }
    }

    // common neighbors are not part of the IS in either case
    for (auto neighbor : graph[node])
    {
        if (is_reduced(neighbor, br_alg) || neighbor == u)
            continue;

        if (u_neighbors.get(neighbor))
        {
            br_alg->set(neighbor, IS_status::excluded);
        }
        else
        {
            if (weights[neighbor] + weights[u] <= weights[node])
                br_alg->set(neighbor, IS_status::excluded);
            else
                remaining_neighbors.push_back(neighbor);
        }
    }

    if (remaining_neighbors.size() == 0)
    {
        assert(br_alg->deg(node) == 1 && "ERROR: funnel_reduction::fold: remaining neighbors must be non-empty");
        if (weights[node] >= weights[u])
        {
            br_alg->set(node, IS_status::included);
        }
        return;
    }

    status.folded_stack.push_back(get_reduction_type());
    status.reduction_offset += weights[node];

    // add additional edges connecting the remaining neighbors with funnel neighors outside the funnel set
    restore_vec.push_back({data.node, u, remaining_neighbors, {}, weights[node] >= weights[u]});
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++)
    {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];

        assert(!is_reduced(neighbor, br_alg) && "ERROR: funnel_reduction::fold: remaining neighbors must be unset");
        assert(node != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        assert(u != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        for (auto neighbor_2nd : outside_funnel_neighbors)
        {
            if (!graph.adjacent(neighbor, neighbor_2nd))
            {
                graph.add_edge_undirected(neighbor, neighbor_2nd);
                restore_vec.back().node_vecs[idx].push_back(neighbor_2nd);
            }
        }
        br_alg->add_next_level_node(neighbor);
    }

    if (weights[u] > weights[node])
    {
        weights[u] -= weights[node];
    }
    else
    {
        for (NodeID n : remaining_neighbors)
        {
            assert(weights[n] >= weights[node] - weights[u] && "ERROR: funnel_reduction::fold: weight must be larger than funnel weight");
            weights[n] += weights[u];
            weights[n] -= weights[node];
        }
        br_alg->set(u, IS_status::folded, false);
    }

    br_alg->set(node, IS_status::folded, true);
}

void funnel_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();
    auto &weights = status.weights;

    br_alg->unset(data.node);
    if (data.node_is_max_weight)
    {
        br_alg->unset(data.funnel_neighbor);
        for (NodeID neighbor : data.remaining_neighbors)
        {
            weights[neighbor] += weights[data.node];
            weights[neighbor] -= weights[data.funnel_neighbor];
        }
    }
    else
    {
        status.weights[data.funnel_neighbor] += status.weights[data.node];
    }

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

    bool include_node = false; 
    if (node_status[data.funnel_neighbor] == IS_status::excluded || 
        node_status[data.funnel_neighbor] == IS_status::folded )
        include_node = std::none_of(remaining.begin(), remaining.end(), [&](NodeID neighbor)
                                    { return node_status[neighbor] == IS_status::included; });

    restore(br_alg);

    if (include_node)
    {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
    }
    else
    {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
    }
    is_weight += weights[data.node];
}

inline int funnel_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label)
{
    if (br_alg->deg(v) <= 2)
        return false; // fold1 or 2
    auto &weights = br_alg->status.weights;
    auto &funnel_set = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    NodeID funnel_neighbor = br_alg->status.n;

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    // need one vertex of weight >= v
    get_neighborhood_set(v, br_alg, funnel_set);
    funnel_set.add(v);

    if (is_funnel(v, funnel_neighbor, br_alg))
    {
        label.push_back(v);
        return true;
    }

    return false;
}