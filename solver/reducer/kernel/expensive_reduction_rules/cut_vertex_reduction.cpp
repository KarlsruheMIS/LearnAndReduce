#include "cut_vertex_reduction.h"

#include "reductions.h"
#include "general_reduction.h"
#include "reduce_algorithm.h"
#include "../tiny_solver/tiny_solver.h"

#include <stack>

typedef reduce_algorithm::IS_status IS_status;

bool cut_vertex_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->blowing_up ||  br_alg->config.disable_cut_vertex)
        return false;

#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    auto &cut_component = br_alg->buffers[0];
    auto &cut_component_set = br_alg->set_1;
    auto &cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto &reverse_mapping = br_alg->buffers[1];
    auto config = br_alg->config;
    size_t oldn = status.remaining_nodes;
    auto &map = br_alg->buffers[2];
    auto &reverse_map = br_alg->buffers[3];
    auto solver = br_alg->subgraph_solver;
    get_mappings_to_remaining_graph(br_alg, map, reverse_map);
    if (articulation_points.size() == 0)
    {
        get_articulation_points(br_alg, articulation_points, reverse_map, map);
        std::sort(articulation_points.begin(), articulation_points.end(), [&](NodeID a, NodeID b)
              { return status.graph[a].size() < status.graph[b].size(); });
    }

    bool progress = false;
    while (articulation_points.size() > 0 && !progress)
    {
        NodeID cut_v = articulation_points.back();
        articulation_points.pop_back();
        if (status.node_status[cut_v] != IS_status::not_set)
            continue;
        if (br_alg->deg(cut_v) < 3)
            break;
        if (!check_components(br_alg, cut_v, cut_component))
            continue;
        if (cut_component.size() <= 2)
            continue;
        assert(std::all_of(cut_component.begin(), cut_component.end(), [&](NodeID v)
                           { return status.node_status[v] == IS_status::not_set; }) &&
               "all nodes need to be unset");
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component)
            cut_component_set.add(n);

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut = false;
        for (NodeID neighbor : status.graph[cut_v])
        {
            if (!cut_component_set.get(neighbor))
            {
                real_cut = true;
                break;
            }
        }

        if (!real_cut)
        { // directly solve the component without fold
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
            if (!solve_induced_subgraph_from_set(no_limit, MWIS_weight, br_alg, cut_component, cut_component_set))
                return false;

            for (NodeID node : cut_component)
            {
                NodeID subgraph_node = solver->forward_map[node];
                if (solver->independent_set[subgraph_node] == 1)
                    br_alg->set(node, IS_status::included);
            }
        }
        else
        {
            NodeWeight small_cutMWIS_weight = 0;
            NodeWeight large_cutMWIS_weight = 0;
            std::vector<NodeID> cut_v_included_i(config.subgraph_node_limit);
            std::vector<NodeID> cut_v_included_e(config.subgraph_node_limit);
            std::vector<NodeID> cut_v_excluded_i(config.subgraph_node_limit);
            std::vector<NodeID> cut_v_excluded_e(config.subgraph_node_limit);
            if (get_fold_data(br_alg, cut_v, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e, large_cutMWIS_weight, small_cutMWIS_weight))
            {
                fold_data data = {cut_v, status.weights[cut_v], large_cutMWIS_weight, small_cutMWIS_weight, cut_component};
                fold(br_alg, data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e);
            }
            progress = true;
        }
    }

    if (articulation_points.size() == 0)
        br_alg->config.disable_cut_vertex = true;

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}

void cut_vertex_reduction::get_mappings_to_remaining_graph(reduce_algorithm *br_alg, std::vector<NodeID> &map, std::vector<NodeID> &reverse_map)
{
    auto &status = br_alg->status;
    NodeID u = 0;
    map.resize(status.n, status.n);
    reverse_map.resize(status.remaining_nodes, status.n);
    for (NodeID n = 0; n < status.n; ++n)
    {
        if (status.node_status[n] == IS_status::not_set)
        {
            map[n] = u;
            reverse_map[u] = n;
            u++;
        }
    }
}

bool cut_vertex_reduction::check_components(reduce_algorithm *br_alg, NodeID u, std::vector<NodeID> &smallComponent)
{
    auto &status = br_alg->status;
    auto &config = br_alg->config;

    // Mark all nodes that are not to be considered as already visited
    std::vector<bool> visited_neighbors(status.n, true);
    for (NodeID n = 0; n < status.n; ++n)
    {
        if (status.node_status[n] == IS_status::not_set)
            visited_neighbors[n] = false;
    }

    visited_neighbors[u] = true;

    for (NodeID n : status.graph[u])
    {
        if (!visited_neighbors[n])
        {
            if (br_alg->config.reduction_time_limit < br_alg->t.elapsed())
                return false;
            smallComponent.clear();

            if (!build_small_component(n, br_alg, smallComponent, visited_neighbors))
            {
                // If component is too large, early terminate and set the visited status for this path
                dfs_fill_visited(n, br_alg, visited_neighbors);
            }
            else if (smallComponent.size() <= config.subgraph_node_limit)
            {
                // A small enough component was found, no need to check further
                return true;
            }
        }
    }
    // After exploring all components and no valid small component was found
    return false;
}
bool cut_vertex_reduction::build_small_component(NodeID u, reduce_algorithm *br_alg, std::vector<NodeID> &component, std::vector<bool> &visited)
{
    auto &status = br_alg->status;
    std::vector<NodeID> stack;
    stack.reserve(status.remaining_nodes);
    stack.push_back(u);

    while (!stack.empty())
    {
        NodeID c = stack.back();
        stack.pop_back();

        if (visited[c])
            continue;

        visited[c] = true;
        assert(status.node_status[c] == IS_status::not_set && "ERROR: cut_vertex_reduction::build_small_component: node already set");
        component.push_back(c);

        if (component.size() > br_alg->config.subgraph_node_limit)
        {
            for (NodeID n : component)
                visited[n] = false;
            component.clear();
            return false;
        }

        for (NodeID neighbor : status.graph[c])
        {
            if (!visited[neighbor])
                stack.push_back(neighbor);
        }
    }
    return true;
}
void cut_vertex_reduction::dfs_fill_visited(NodeID u, reduce_algorithm *br_alg, std::vector<bool> &visited)
{
    auto &status = br_alg->status;
    auto &visited_dfs = br_alg->bool_buffer;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed())
    {
        NodeID c = stack.back();
        stack.pop_back();

        if (visited[c])
        {
            continue;
        }

        // Mark the current node as visited
        visited[c] = true;
        visited_dfs[c] = true;

        for (NodeID neighbor : br_alg->status.graph[c])
        {
            if (!visited[neighbor])
            {
                stack.push_back(neighbor);
            }
        }
    }
}
void cut_vertex_reduction::fold(reduce_algorithm *br_alg, fold_data &data, std::vector<NodeID> &cut_v_included_i, std::vector<NodeID> &cut_v_included_e, std::vector<NodeID> &cut_v_excluded_i, std::vector<NodeID> &cut_v_excluded_e)
{
    auto &status = br_alg->status;
    restore_vec.push_back({data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e});
    for (int i = data.cut_component.size() - 1; i >= 1; i--)
    {
        br_alg->set(data.cut_component[i], IS_status::folded, false);
    }
    br_alg->set(data.cut_component[0], IS_status::folded);

    status.reduction_offset += data.large_cutMWIS_weight;
    status.weights[data.cut_vertex] += data.small_cutMWIS_weight;
    status.weights[data.cut_vertex] -= data.large_cutMWIS_weight;

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(data.cut_vertex);
    br_alg->add_next_level_neighborhood(data.cut_vertex);
}
void cut_vertex_reduction::restore(reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back().data;

    status.reduction_offset -= data.large_cutMWIS_weight;
    status.weights[data.cut_vertex] = data.cut_vertex_weight;

    for (size_t i = 0; i < data.cut_component.size(); i++)
    {
        br_alg->unset(data.cut_component[i]);
    }

    restore_vec.pop_back();
}
void cut_vertex_reduction::apply(reduce_algorithm *br_alg)
{
    auto &node_status = br_alg->status.node_status;
    NodeWeight &is_weight = br_alg->status.is_weight;
    auto restore_info = restore_vec.back();
    NodeID cut_v = restore_info.data.cut_vertex;
    assert(restore_info.data.large_cutMWIS_weight < restore_info.data.cut_vertex_weight + restore_info.data.small_cutMWIS_weight && "ERROR: cut_vertex_reduction::apply: large_cutMWIS_weight must be larger than small_cutMWIS_weight");

    IS_status cut_status = node_status[cut_v];

    restore(br_alg);

    if (cut_status == IS_status::included)
    {
        for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_exclude)
        {
            node_status[neighbor] = IS_status::excluded;
        }
        for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_include)
        {
            node_status[neighbor] = IS_status::included;
        }
        is_weight += restore_info.data.small_cutMWIS_weight;
        is_weight += restore_info.data.cut_vertex_weight;
    }
    else
    {
        for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_exclude)
        {
            node_status[neighbor] = IS_status::excluded;
        }
        for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_include)
        {
            node_status[neighbor] = IS_status::included;
        }
        is_weight += restore_info.data.large_cutMWIS_weight;
    }
}

bool cut_vertex_reduction::get_fold_data(reduce_algorithm *br_alg, NodeID cut_v, std::vector<NodeID> &cut_v_included_i, std::vector<NodeID> &cut_v_included_e, std::vector<NodeID> &cut_v_excluded_i, std::vector<NodeID> &cut_v_excluded_e, NodeWeight &large_cutMWIS_weight, NodeWeight &small_cutMWIS_weight)
{
    auto &status = br_alg->status;
    auto &config = br_alg->config;
    auto &cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto &cut_component_set = br_alg->set_1;
    auto &cut_component = br_alg->buffers[0];
    tiny_solver *solver = br_alg->subgraph_solver;
    NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();

    if (!solve_induced_subgraph_from_set(no_limit, large_cutMWIS_weight, br_alg, cut_component, cut_component_set))
        return false;

    // save solution for later reduce exclude cut_v branch:
    cut_v_excluded_e.clear();
    cut_v_excluded_i.clear();
    cut_v_excluded_e.push_back(cut_v);
    assert(!std::all_of(cut_component.begin(), cut_component.end(), [&](NodeID v)
                        { return solver->independent_set[solver->forward_map[v]] == 1; }) &&
           "all nodes need to be unset");
    for (NodeID node : cut_component)
    {
        NodeID subgraph_node = solver->forward_map[node];
        if (solver->independent_set[subgraph_node] == 1)
            cut_v_excluded_i.push_back(node);
        else
            cut_v_excluded_e.push_back(node);
    }

    // set weight of neighborhood of cut_v to 0, ie solve G-N(cut_v)
    for (NodeID neighbor : status.graph[cut_v])
    {
        if (cut_component_set.get(neighbor))
            solver->subgraph_W[solver->forward_map[neighbor]] = 0;
    }
    tiny_solver_clear_solution(solver);
    tiny_solver_solve(solver, config.subgraph_time_limit, no_limit);
    if (solver->time_limit_exceeded || solver->node_limit_exceeded)
        return false;
    small_cutMWIS_weight = solver->independent_set_weight;

    if (status.weights[cut_v] + small_cutMWIS_weight <= large_cutMWIS_weight) // cut_vertex is excluded -> directly apply
    {
        for (NodeID node : cut_v_excluded_i)
        {
            br_alg->set(node, IS_status::included);
        }
        return false;
    }
    else
    {
        cut_v_included_e.clear();
        cut_v_included_i.clear();
        cut_v_included_i.push_back(cut_v);

        for (NodeID node : cut_component)
        {
            if (solver->independent_set[solver->forward_map[node]] == 1)
            {
                assert(!cut_v_neighbor_set.get(node));
                cut_v_included_i.push_back(node);
            }
            else
            {
                cut_v_included_e.push_back(node);
            }
        }
        return true;
    }
}
void cut_vertex_reduction::get_articulation_points(reduce_algorithm *br_alg, std::vector<NodeID> &articulation_points, std::vector<NodeID> &reverse_map, std::vector<NodeID> &map)
{
    auto &status = br_alg->status;
    auto &graph = br_alg->status.graph;
    auto &visited = br_alg->bool_buffer;
    auto &low = br_alg->buffers[0];
    auto &disc = br_alg->buffers[1];
    auto &articulation_point_set = br_alg->set_1;
    std::vector<NodeID> parent(status.n, status.n + 1);
    low.assign(status.n, status.n);
    disc.assign(status.n, status.n);
    visited.assign(status.n, false);
    articulation_point_set.clear();
    std::stack<std::pair<NodeID, NodeID>> dfsStack; // Stack to simulate recursion (pair of node and child index)

    int time = 0;         // To track discovery times
    int rootChildren = 0; // To track the number of children of the root

    // Start DFS for all unvisited nodes (in case of disconnected graph)
    for (int root = 0; root < status.remaining_nodes; ++root)
    {
        assert(status.node_status[reverse_map[root]] == IS_status::not_set);
        if (disc[root] == status.n)
        {
            // Initialize the DFS with the root
            dfsStack.push({root, 0});
            parent[root] = status.n;
            disc[root] = low[root] = time++;
            rootChildren = 0;

            while (!dfsStack.empty())
            {
                auto [u, childIndex] = dfsStack.top();
                dfsStack.pop();

                // If we haven't processed all neighbors
                if (childIndex < status.graph[reverse_map[u]].size())
                {
                    NodeID v = map[status.graph[reverse_map[u]][childIndex]];
                    assert(status.node_status[reverse_map[v]] == IS_status::not_set);
                    assert(v < status.remaining_nodes);
                    dfsStack.push({u, childIndex + 1}); // Continue to the next child in the future

                    if (disc[v] == status.n)
                    {
                        // v is an unvisited child of u
                        parent[v] = u;
                        disc[v] = low[v] = time++;
                        dfsStack.push({v, 0}); // Process v next
                        if (parent[u] == status.n)
                            rootChildren++; // If u is the root, increment rootChildren
                    }
                    else if (v != parent[u])
                    {
                        // v is a back edge
                        low[u] = std::min(low[u], disc[v]);
                    }
                }
                else
                {
                    // All children of u have been processed
                    if (parent[u] != status.n)
                    {
                        if (parent[parent[u]] == status.n)
                            continue;
                        // Check articulation point condition for non-root nodes
                        if (low[u] >= disc[parent[u]] && !articulation_point_set.get(reverse_map[parent[u]]))
                        {
                            articulation_points.push_back(reverse_map[parent[u]]);
                            articulation_point_set.add(reverse_map[parent[u]]);
                        }
                        // Update the low value of the parent
                        low[parent[u]] = std::min(low[parent[u]], low[u]);
                    }
                }
            }

            // Special case for root
            if (parent[root] == status.n && rootChildren > 1)
            {
                articulation_points.push_back(reverse_map[root]);
            }
        }
    }
    assert(articulation_points.size() == 0 || std::any_of(articulation_points.begin(), articulation_points.end(), [&](NodeID v)
                                                          { return status.node_status[v] == IS_status::not_set; }) &&
                                                  "ERROR: cut_vertex_reduction::get_articulation_points: articulation point already set");
}
void cut_vertex_reduction::generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index)
{
    auto &status = br_alg->status;
    auto &map = br_alg->buffers[2];
    auto &reverse_map = br_alg->buffers[3];
    std::vector<NodeID> articulation_points;
    get_mappings_to_remaining_graph(br_alg, map, reverse_map);
    get_articulation_points(br_alg, articulation_points, reverse_map, map);
    for (NodeID node : articulation_points)
        reduction_data[reduction_index][node] = 1;
}