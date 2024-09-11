#include "cut_vertex_reduction.h"

#include "reductions.h"
#include "general_reduction.h"
#include "branch_and_reduce_algorithm.h"

#include <stack>

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool cut_vertex_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->blowing_up) return false;
    if (br_alg->heuristically_reducing) return false;
    // if (br_alg->config.disable_cut_vertex) return false;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif
	auto& status = br_alg->status;
    auto& visited = br_alg->bool_buffer;
    auto& cut_component = br_alg->buffers[0];
    auto& cut_component_set = br_alg->set_1;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& reverse_mapping = br_alg->buffers[1];
    auto config = br_alg->config;
    size_t oldn = status.remaining_nodes;
    NodeID cut_v;
    std::vector<NodeID> unreduced_mapping(status.remaining_nodes, status.n);
    fast_set tested_cut_vertices(status.n);

    // TODO make dfs only run on reduced nodes (mappings vie buffer)
    visited.assign(status.n, true);
    NodeID local_n = 0;
    for (NodeID n = 0; n < br_alg->status.n; ++n)
    {
        if (status.node_status[n] == IS_status::not_set)
            visited[n] = false;
        if (!visited[n])
        {
            unreduced_mapping[local_n] = n;
            local_n++;
        }
        assert(local_n <= status.remaining_nodes && "ERROR: cut_vertex_reduction::reduce: local_n too small");
    }

    while (find_cut_vertex(br_alg, cut_v, cut_component, unreduced_mapping, tested_cut_vertices) && br_alg->config.reduction_time_limit > br_alg->t.elapsed())
    {
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1)
            return false; // fold1
        assert(cut_component.size() <= config.subgraph_node_limit);

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component)
        {
            cut_component_set.add(n);
        }

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut_v = false;
        for (NodeID neighbor : status.graph[cut_v])
        {
            if (!cut_component_set.get(neighbor))
            {
                real_cut_v = true;
                break;
            }
        }

        if (!real_cut_v)
        { // directly solve the component without fold
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
            graph_access component_graph;
            if (!solve_induced_subgraph_from_set(no_limit, MWIS_weight, component_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
                return false;

            for (NodeID node : cut_component)
            {
                assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
                if (component_graph.getPartitionIndex(reverse_mapping[node]) == 1)
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
        }
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    br_alg->config.disable_cut_vertex = true; // only test whole graph once
    return oldn != status.remaining_nodes;
}
bool cut_vertex_reduction::find_cut_vertex(branch_and_reduce_algorithm* br_alg, NodeID& cut_v, std::vector<NodeID>& cut_component, std::vector<NodeID>& reverse_mapping, fast_set& tested) {
    auto& visited = br_alg->bool_buffer;

    int step = 0;
    for (NodeID local_n = 0; local_n < br_alg->status.remaining_nodes; ++local_n)
    {
        if (!visited[reverse_mapping[local_n]])
        {
            if (DFS(br_alg, reverse_mapping[local_n], step, cut_v, cut_component) && !tested.get(cut_v))
            {
                cut_v = cut_v;
                tested.add(cut_v);
                return true; // New articulation point with small component found
            }
            if (br_alg->config.reduction_time_limit < br_alg->t.elapsed())
                return false;
        }
    }
    return false; // No articulation point meeting criteria found
}
bool cut_vertex_reduction::DFS(branch_and_reduce_algorithm* br_alg, NodeID u, int& step, NodeID& cut_vertex, std::vector<NodeID>& smallComponent) {
    auto& status = br_alg->status;
    auto& graph = br_alg->status.graph;
    auto& visited = br_alg->bool_buffer;
    auto& low = br_alg->buffers[0];
    auto& disc = br_alg->buffers[1];
    auto& parent = br_alg->buffers[2];
    low.assign(status.n,std::numeric_limits<NodeID>::max());
    disc.assign(status.n,0);
    parent.assign(status.n, status.n);
    std::vector<NodeID> stack;
    std::vector<std::pair<NodeID, NodeID>> edge_stack;
    stack.reserve(br_alg->status.remaining_nodes);
    edge_stack.reserve(br_alg->status.remaining_nodes);
    stack.push_back(u);

    while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed())
    {
        NodeID c = stack.back();

        disc[c] = low[c] = ++step;
        visited[c] = true;

        bool pushedChild = false;
        for (NodeID v : graph[c])
        {
            // Avoid reprocessing an edge in the undirected graph
            if (!edge_stack.empty() && edge_stack.back() == std::make_pair(v, c))
                continue;

            if (!visited[v])
            {
                stack.push_back(v);
                edge_stack.push_back(std::make_pair(c, v)); // Track the edge
                parent[v] = c;
                pushedChild = true;
                break;
            }
            else if (v != parent[c])
            {
                low[c] = std::min(low[c], disc[v]);
            }
        }
        if (!pushedChild)
        {
            stack.pop_back(); // Finish processing current node
            if (!stack.empty())
            {
                low[parent[c]] = std::min(low[parent[c]], low[c]);

                // check if current is a cut vertex
                if ((parent[c] == status.n && disc[c] != low[c]) || // root node
                    (parent[c] != status.n && low[c] >= disc[parent[c]]))
                    {
                        if (check_components(br_alg, c, cut_vertex, smallComponent))
                        {
                            return true;
                        }
                    }
            }
            if (!edge_stack.empty())
                edge_stack.pop_back(); // Pop the edge after processing
        }
    }
    return false;
}
bool cut_vertex_reduction::check_components(branch_and_reduce_algorithm* br_alg, NodeID u, NodeID& cut_vertex, std::vector<NodeID>& smallComponent) 
{
    auto& status = br_alg->status;
    auto& config = br_alg->config;

    // Mark all nodes that are not to be considered as already visited
    std::vector<bool> visited_neighbors(status.n, true);
    for (NodeID n = 0; n < status.n; ++n)
    {
        if (status.node_status[n] == IS_status::not_set)
        {
            visited_neighbors[n] = false;
        }
    }

    visited_neighbors[u] = true;
    cut_vertex = u;

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
bool cut_vertex_reduction::build_small_component(NodeID u, branch_and_reduce_algorithm* br_alg, std::vector<NodeID>& component, std::vector<bool>& visited) 
{
    auto &status = br_alg->status;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty())
    {
        NodeID c = stack.back();
        stack.pop_back();

        if (visited[c])
            continue;

        visited[c] = true;
        component.push_back(c);

        // If the component exceeds the maximum size, return false to indicate failure
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
            {
                stack.push_back(neighbor);
            }
        }
    }
    return true;
}
void cut_vertex_reduction::dfs_fill_visited(NodeID u, branch_and_reduce_algorithm *br_alg, std::vector<bool> &visited)
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
void cut_vertex_reduction::fold(branch_and_reduce_algorithm* br_alg, fold_data& data, std::vector<NodeID>& cut_v_included_i, std::vector<NodeID>& cut_v_included_e, std::vector<NodeID>& cut_v_excluded_i, std::vector<NodeID>& cut_v_excluded_e)
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
void cut_vertex_reduction::restore(branch_and_reduce_algorithm *br_alg)
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
void cut_vertex_reduction::apply(branch_and_reduce_algorithm *br_alg)
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

inline bool cut_vertex_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    if (br_alg->blowing_up)
        return false;
    // if (br_alg->config.disable_cut_vertex) return false;
    auto &status = br_alg->status;
    auto &cut_component = br_alg->buffers[0];
    auto &reverse_mapping = br_alg->buffers[1];
    auto &cut_component_set = br_alg->set_1;
    auto &cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto &config = br_alg->config;
    size_t oldn = status.remaining_nodes;
    NodeID cut_v;
    if (check_components(br_alg, v, cut_v, cut_component))
    {
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1)
            return false; // fold1
        assert(cut_component.size() <= config.subgraph_node_limit);

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component)
        {
            cut_component_set.add(n);
        }

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut_v = false;
        for (NodeID neighbor : status.graph[cut_v])
        {
            if (!cut_component_set.get(neighbor))
                real_cut_v = true;
        }

        if (!real_cut_v)
        { // directly solve the component without fold
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
            graph_access component_graph;
            if (!solve_induced_subgraph_from_set(no_limit, MWIS_weight, component_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
                return false;

            for (NodeID node : cut_component)
            {
                assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
                if (component_graph.getPartitionIndex(reverse_mapping[node]) == 1)
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
        }
    }

    return oldn != status.remaining_nodes;
}
bool cut_vertex_reduction::get_fold_data(branch_and_reduce_algorithm* br_alg, NodeID cut_v, std::vector<NodeID>& cut_v_included_i, std::vector<NodeID>& cut_v_included_e, std::vector<NodeID>& cut_v_excluded_i, std::vector<NodeID>& cut_v_excluded_e, NodeWeight& large_cutMWIS_weight, NodeWeight& small_cutMWIS_weight)
{
	auto& status = br_alg->status;
    auto& config = br_alg->config;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& cut_component_set = br_alg->set_1;
    auto& reverse_mapping = br_alg->buffers[1];
    auto& cut_component = br_alg->buffers[0];
    graph_access cut_graph;
    reverse_mapping.assign(status.n, status.n);
    NodeWeight bound = std::numeric_limits<NodeWeight>::max();

    // graph including neighborhood of cut_v
    if (!solve_induced_subgraph_from_set(bound, large_cutMWIS_weight, cut_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
        return false;

    // save solution for later reduce
    cut_v_excluded_e.clear();
    cut_v_excluded_i.clear();
    cut_v_excluded_e.push_back(cut_v);
    for (NodeID node : cut_component)
    {
        assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
        if (cut_graph.getPartitionIndex(reverse_mapping[node]) == 1)
            cut_v_excluded_i.push_back(node);
        else
            cut_v_excluded_e.push_back(node);
    }

    // set weight of neighborhood of cut_v to 0, ie solve G-N(cut_v)
    for (NodeID neighbor : status.graph[cut_v])
    {
        if (reverse_mapping[neighbor] == status.n)
            continue;
        cut_graph.setNodeWeight(reverse_mapping[neighbor], 0);
        cut_graph.setPartitionIndex(reverse_mapping[neighbor], 0);
    }
    if (!solve_graph(small_cutMWIS_weight, cut_graph, config, bound, true))
        return false;

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
            if (cut_v_neighbor_set.get(node))
            {
                cut_v_included_e.push_back(node);
                continue;
            }
            assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
            if (cut_graph.getPartitionIndex(reverse_mapping[node]) == 1)
            {
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
bool cut_vertex_reduction::generate_global_data(branch_and_reduce_algorithm* br_alg, std::vector<NodeID>& articulation_points)
{
 
    auto& status = br_alg->status;
    auto& graph = br_alg->status.graph;
    auto& visited = br_alg->bool_buffer;
    auto& low = br_alg->buffers[0];
    auto& disc = br_alg->buffers[1];
    auto& parent = br_alg->buffers[2];
    low.assign(status.n,status.n);
    disc.assign(status.n,status.n);
    parent.assign(status.n, status.n+1);
    visited.assign(status.n, false);
    std::stack<std::pair<NodeID, NodeID>> dfsStack;  // Stack to simulate recursion (pair of node and child index)
    
    int time = 0;  // To track discovery times
    int rootChildren = 0;  // To track the number of children of the root

    // Start DFS for all unvisited nodes (in case of disconnected graph)
    for (int root = 0; root < status.remaining_nodes; ++root) {
        if (disc[root] == status.n) {
            // Initialize the DFS with the root
            dfsStack.push({root, 0});
            parent[root] = status.n;
            disc[root] = low[root] = time++;
            rootChildren = 0;

            while (!dfsStack.empty()) {
                auto [u, childIndex] = dfsStack.top();
                dfsStack.pop();

                // If we haven't processed all neighbors
                if (childIndex < status.graph[u].size()) {
                    int v = status.graph[u][childIndex];
                    dfsStack.push({u, childIndex + 1});  // Continue to the next child in the future

                    if (disc[v] == status.n) {
                        // v is an unvisited child of u
                        parent[v] = u;
                        disc[v] = low[v] = time++;
                        dfsStack.push({v, 0});  // Process v next
                        if (parent[u] == status.n) rootChildren++;  // If u is the root, increment rootChildren
                    } else if (v != parent[u]) {
                        // v is a back edge
                        low[u] = std::min(low[u], disc[v]);
                    }
                } else {
                    // All children of u have been processed
                    if (parent[parent[u]] != status.n && parent[u] != status.n) {
                        // Check articulation point condition for non-root nodes
                        if (low[u] >= disc[parent[u]]) {
                            articulation_points.push_back(parent[u]);
                        }
                        // Update the low value of the parent
                        low[parent[u]] = std::min(low[parent[u]], low[u]);
                    }
                }
            }

            // Special case for root
            if (parent[root] == status.n && rootChildren > 1) {
                articulation_points.push_back(root);
            }
        }
    }
    return articulation_points.size() > 0;
}