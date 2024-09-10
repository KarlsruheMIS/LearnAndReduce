#include "generalized_fold_reduction.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    // if (br_alg->config.disable_generalized_fold) return false;
    if (br_alg->heuristically_reducing) return false;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            {
        NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
        if(try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;
        reduce_vertex(br_alg, v); });

    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
    // has_run = false; // check marker with gnn in next round again
    // has_filtered_marker = true;
	return oldn != status.remaining_nodes;
}
inline bool generalized_fold_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_generalized_fold) return false;
    if (br_alg->deg(v) <= 1)
        return false;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &neighbors_set = br_alg->set_1;
    auto &MWIS_set = br_alg->set_2;
    auto &reverse_mapping = br_alg->buffers[1];
    size_t oldn = status.remaining_nodes;

    graph_access neighborhood_graph;

    get_neighborhood_set(v, br_alg, neighbors_set);
    get_neighborhood_vector(v, br_alg, neighbors);

    /* NodeID max_neighbor = get_max_weight_neighbor(v, br_alg); */
    /* NodeWeight max_neighbor_weight = status.weights[max_neighbor]; */

    if (status.graph[v].size() > br_alg->config.subgraph_node_limit)
    {
        return false;
    }

    // compute MWIS in N(v)

    NodeWeight MWIS_weight = 0;
    NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();
    bool solved_exact = solve_induced_subgraph_from_set(min_MWIS_neighbor_weight, MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping, true);
    if (!solved_exact)
    {
#ifdef REDUCTION_INFO
        reduction_time += br_alg->reduction_timer.elapsed();
#endif
        return false;
    }

    if (status.weights[v] >= MWIS_weight)
    {
        // }
        // same as in generalized_neighborhood_reduction
        br_alg->set(v, IS_status::included);
#ifdef REDUCTION_INFO
        reduction_time += br_alg->reduction_timer.elapsed();
#endif
        return oldn != status.remaining_nodes;
    }

    MWIS_set.clear();

    forall_nodes(neighborhood_graph, node)
    {
        if (neighborhood_graph.getPartitionIndex(node) == 1)
        {
            const NodeID neighbor = neighbors[node];
            MWIS_set.add(neighbor);

            if (status.weights[neighbor] < min_MWIS_neighbor_weight)
                min_MWIS_neighbor_weight = status.weights[neighbor];
        }
    }
    endfor

        if (status.weights[v] < MWIS_weight - min_MWIS_neighbor_weight)
    {
        // multiple IS exist that have bigger weight than v
        return false;
    }

    bool check_failed = false;

    // check that no other IS in N(v) exists with weight greater than v
    for (const NodeID neighbor : status.graph[v])
    {
        if (!MWIS_set.get(neighbor))
            continue;

        auto iter = std::find(neighbors.begin(), neighbors.end(), neighbor);
        assert(iter != neighbors.end());
        std::swap(*iter, neighbors.back());
        neighbors.pop_back();
        neighbors_set.remove(neighbor);

        NodeWeight MWIS_weight = 0;
        bool solved_exact = solve_induced_subgraph_from_set(status.weights[v], MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping);
        if (!solved_exact)
        {
            check_failed = true;
        }
        else if (MWIS_weight >= status.weights[v])
        {
            check_failed = true;
        }

        neighbors.push_back(neighbor);
        neighbors_set.add(neighbor);

        if (check_failed)
            break;
    }

    if (!check_failed)
    {
        fold(br_alg, v, MWIS_set, MWIS_weight);
        return oldn != status.remaining_nodes;
    }

    auto &neighborhood_intersection_set = MWIS_set;
    bool remove_node;

    // we can't fold but we can possibly remove some neighbors of v
    do
    {
        for (const NodeID node : status.graph[v])
        {
            neighborhood_intersection_set.clear();

            for (const NodeID neighbor : status.graph[node])
            {
                if (neighbors_set.get(neighbor))
                {
                    neighborhood_intersection_set.add(neighbor);
                }
            }

            // "force" node into an IS (= remove node and its neighbors from N(v) and compute MWIS in remaining N(v))
            auto iter = std::find(neighbors.begin(), neighbors.end(), node);
            assert(iter != neighbors.end());
            std::swap(*iter, neighbors.back());
            neighbors.pop_back();
            neighbors_set.remove(node);

            for (const NodeID neighbor : status.graph[node]) {
                if (neighborhood_intersection_set.get(neighbor)) {
                    auto iter = std::find(neighbors.begin(), neighbors.end(), neighbor);
                    assert(iter != neighbors.end());
                    std::swap(*iter, neighbors.back());
                    neighbors.pop_back();
                    neighbors_set.remove(neighbor);
                }
            }

            if (status.weights[v] < status.weights[node])
            {
                remove_node = false;
            }
            else
            {
                NodeWeight MWIS_weight = 0;
                NodeWeight bound = status.weights[v] - status.weights[node];
                bool solved_exact = solve_induced_subgraph_from_set(bound, MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping);
                if (!solved_exact)
                {
                    remove_node = false;
                }
                else
                {
                    // if the weight of every MWIS in N(v) which contains "node" is smaller than w(v) then we can remove "node"
                    remove_node = MWIS_weight <= bound;
                }
            }

            for (const NodeID neighbor : status.graph[node])
            {
                if (neighborhood_intersection_set.get(neighbor))
                {
                    neighbors.push_back(neighbor);
                    neighbors_set.add(neighbor);
                }
            }

            if (remove_node)
            {
                br_alg->set(node, IS_status::excluded);
                break; // break and restart loop because set(..) modifies the range which we currently iterate
            }

            neighbors.push_back(node);
            neighbors_set.add(node);
        }
    } while (remove_node);

    return oldn != status.remaining_nodes;
}
void generalized_fold_reduction::fold(branch_and_reduce_algorithm *br_alg, NodeID main_node, fast_set &MWIS_set, NodeWeight MWIS_weight)
{
    auto &status = br_alg->status;

    restore_vec.emplace_back();
    restore_data &data = restore_vec.back();
    data.main_weight = status.weights[main_node];
    data.MWIS_weight = MWIS_weight;

    auto &nodes = data.nodes;
    nodes.main = main_node;

    // temporary copy for iteration
    data.main_neighbor_list = status.graph[main_node];

    for (auto neighbor : data.main_neighbor_list)
    {
        if (MWIS_set.get(neighbor))
            nodes.MWIS.push_back(neighbor);
        else
            br_alg->set(neighbor, IS_status::excluded);
    }

    // reverse order because of later "restore_edge_and_replace"
    for (int i = nodes.MWIS.size() - 1; i >= 1; i--)
    {
        br_alg->set(nodes.MWIS[i], IS_status::folded, false);
    }

    br_alg->set(nodes.MWIS[0], IS_status::folded, true);

    data.main_neighbor_list = status.graph[main_node];

    // "move" weight into redu offset
    status.reduction_offset += data.main_weight;

    status.weights[nodes.main] = MWIS_weight - data.main_weight;

    std::vector<NodeID> new_neighbors;
    auto &neighbors = MWIS_set;
    neighbors.clear();
    neighbors.add(main_node);

    for (NodeID MWIS_node : nodes.MWIS)
    {
        std::vector<NodeID> node_vec;

        for (auto neighbor : status.graph[MWIS_node])
        {
            if (neighbors.add(neighbor))
            {
                new_neighbors.push_back(neighbor);
                status.graph.add_edge_directed(neighbor, nodes.main);
                node_vec.push_back(neighbor);
            }
        }

        data.MWIS_node_vecs.push_back(std::move(node_vec));
    }

    status.graph[nodes.main] = dynamic_graph::neighbor_list(std::move(new_neighbors));
    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.main);
    br_alg->add_next_level_neighborhood(nodes.main);
}
void generalized_fold_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    // is "restored" in following loops
    status.graph.hide_node(data.nodes.main);
    status.graph[data.nodes.main] = std::move(data.main_neighbor_list);

    for (size_t i = 0; i < data.nodes.MWIS.size(); i++)
    {
        br_alg->unset(data.nodes.MWIS[i]);

        for (auto neighbor : data.MWIS_node_vecs[i])
        {
            status.graph.replace_last_restored_edge(neighbor, data.nodes.MWIS[i]);
        }
    }

    status.weights[data.nodes.main] = data.main_weight;
    status.reduction_offset -= data.main_weight;

    restore_vec.pop_back();
}
void generalized_fold_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto nodes = restore_vec.back().nodes;
    auto MWIS_weight = restore_vec.back().MWIS_weight;
    auto main_status = status.node_status[nodes.main];
    restore(br_alg);

    if (main_status == IS_status::included)
    {
        status.node_status[nodes.main] = IS_status::excluded;

        for (auto node : nodes.MWIS)
        {
            status.node_status[node] = IS_status::included;
        }

        status.is_weight += MWIS_weight;
    }
    else
    {
        status.node_status[nodes.main] = IS_status::included;

        for (auto node : nodes.MWIS)
        {
            status.node_status[node] = IS_status::excluded;
        }

        status.is_weight += status.weights[nodes.main];
    }
}
inline int generalized_fold_reduction::generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    if (br_alg->deg(v) <= 1)
        return 0;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &neighbors_set = br_alg->set_1;
    auto &MWIS_set = br_alg->set_2;
    auto &reverse_mapping = br_alg->buffers[1];
    size_t oldn = status.remaining_nodes;

    graph_access neighborhood_graph;

    get_neighborhood_set(v, br_alg, neighbors_set);
    get_neighborhood_vector(v, br_alg, neighbors);

    if (status.graph[v].size() > br_alg->config.subgraph_node_limit)
    {
        label.push_back(v);
        return 2;
    }

    // compute MWIS in N(v)

    NodeWeight MWIS_weight = 0;
    NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();
    bool solved_exact = solve_induced_subgraph_from_set(min_MWIS_neighbor_weight, MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping, true);
    if (!solved_exact)
    {
        label.push_back(v);
        return 2;
    }

    if (status.weights[v] >= MWIS_weight)
    {
        label.push_back(v);
        return 1;
    }

    MWIS_set.clear();

    forall_nodes(neighborhood_graph, node)
    {
        if (neighborhood_graph.getPartitionIndex(node) == 1)
        {
            const NodeID neighbor = neighbors[node];
            MWIS_set.add(neighbor);

            if (status.weights[neighbor] < min_MWIS_neighbor_weight)
                min_MWIS_neighbor_weight = status.weights[neighbor];
        }
    }
    endfor

    if (status.weights[v] < MWIS_weight - min_MWIS_neighbor_weight)
    {
        // multiple IS exist that have bigger weight than v
        label.push_back(v);
        return 0;
    }

    bool check_failed = false;
    int l = 0;

    // check that no other IS in N(v) exists with weight greater than v
    for (const NodeID neighbor : status.graph[v])
    {
        if (!MWIS_set.get(neighbor))
            continue;

        auto iter = std::find(neighbors.begin(), neighbors.end(), neighbor);
        assert(iter != neighbors.end());
        std::swap(*iter, neighbors.back());
        neighbors.pop_back();
        neighbors_set.remove(neighbor);

        NodeWeight MWIS_weight = 0;
        bool solved_exact = solve_induced_subgraph_from_set(status.weights[v], MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping,l);
        if (!solved_exact)
        {
            check_failed = true;
        }
        else if (MWIS_weight >= status.weights[v])
        {
            l = 0;
            check_failed = true;
        }

        neighbors.push_back(neighbor);
        neighbors_set.add(neighbor);

        if (check_failed)
            break;
    }

    if (!check_failed)
    {
        for (const NodeID neighbor : status.graph[v])
        {
            if (MWIS_set.get(neighbor))
                label.push_back(neighbor);
        }
        return 1;
    }

    return l;
}