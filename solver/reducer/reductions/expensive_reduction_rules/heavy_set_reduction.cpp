
#include "heavy_set_reduction.h"
#include "tiny_solver.h"

#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool heavy_set_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_heavy_set || br_alg->blowing_up)
        return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;

    if (marker.current_size() == 0 )
    {
        has_filtered_marker = false;
        br_alg->config.disable_heavy_set = true;
        return false;
    }
    // find common vertices from gnn suggested vertices in marker
    if (compute_common_vertices && has_filtered_marker)
    {
        if (br_alg->config.gnn_filter != ReductionConfig::GNN_Filter_Type::ALWAYS)
            compute_common_vertices = false;
        auto &vertices = br_alg->buffers[0];
        vertices.clear();
        vertices.assign(status.n, 0);
        for (NodeID v : marker.current)
            for (NodeID n : status.graph[v])
                vertices[n]++;
        marker.clear_next();
        for (NodeID v = 0 ; !is_reduced(v, br_alg) && v < vertices.size(); v++)
        {
            if (vertices[v] > 1)
                marker.add(v);
        }
        marker.get_next();
        if (marker.current_size() == 0 )
            return false;
    }
    size_t oldn = status.remaining_nodes;
    bool progress = false;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            {
        if (br_alg->t.elapsed() > br_alg->config.time_limit) return;
        // jump back to simple rules if progress and reduction not disabled
        if (reduce_vertex(br_alg, v) && !br_alg->config.disable_heavy_set) 
        {
            gnn_filter_marker(br_alg->config, v);
            progress = true;
            return;
        } });
    
    if (!progress)
        gnn_filter_marker(br_alg->config, marker.current.back());
    
#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool heavy_set_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID common_neighbor)
{
    auto &config = br_alg->config;
    if (br_alg->deg(common_neighbor) < 3)
        return false; // use other reduction

    auto &status = br_alg->status;
    assert(status.node_status[common_neighbor] == IS_status::not_set && "ERROR: heavy_set_reduction::reduce_vertex: node must be unset");
    // check if common neighbor is suitable
    if (std::all_of(status.graph[common_neighbor].begin(), status.graph[common_neighbor].end(), [&](NodeID neighbor)
                    { return br_alg->deg(neighbor) > config.subgraph_node_limit; }))
        return false;
    size_t oldn = status.remaining_nodes;
    auto &v_neighbors_set = br_alg->set_1;
    auto &candidates = status.graph[common_neighbor];
    for (size_t v_idx = 0; v_idx < candidates.size() - 1; v_idx++)
    {
        NodeID v = candidates[v_idx];
        if (is_reduced(v, br_alg))
            continue;
        if (br_alg->deg(v) > config.subgraph_node_limit)
            continue; // too many neighbors
        if (br_alg->deg(v) < 3)
            continue; // use other reduction
        NodeID deg_v = br_alg->deg(v);
        get_neighborhood_set(v, br_alg, v_neighbors_set);

        // find second heavy vertex u (not adjacent to v)
        for (size_t u_idx = v_idx + 1; u_idx < candidates.size(); u_idx++)
        {
            NodeID u = candidates[u_idx];
            assert(u != v);
            if (v_neighbors_set.get(u))
                continue; // look for non adjacent nodes
            if (is_reduced(u, br_alg))
                continue;
            if (br_alg->deg(u) + deg_v > config.subgraph_node_limit)
                continue; // subgraph too large
            if (br_alg->deg(u) < 3)
                continue; // use other reduction

            if (is_heavy_set(v, v_neighbors_set, u, br_alg))
            {
                // reduction was applied successfully
                return oldn != status.remaining_nodes;
            }
        }
    }

    return oldn != status.remaining_nodes;
}
bool heavy_set_reduction::is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, reduce_algorithm *br_alg)
{
    auto &config = br_alg->config;
    auto &status = br_alg->status;
    auto &weights = status.weights;
    auto &v_neighbors_vec = br_alg->buffers[0];
    auto &u_neighbors_vec = br_alg->buffers[1];
    auto &graph_nodes = br_alg->buffers[2];
    auto &graph_nodes_set = br_alg->double_set;
    auto solver = br_alg->subgraph_solver;
    size_t oldn = status.remaining_nodes;
    graph_nodes.clear();
    graph_nodes_set.clear();
    graph_access subgraph;

    get_neighborhood_vector(v, br_alg, v_neighbors_vec);
    get_neighborhood_vector(u, br_alg, u_neighbors_vec);

    assert(u != v && "ERROR: heavy_set_reduction::is_heavy_set: u and v must be different");
    assert(!status.graph.adjacent(v, u) && "ERROR: heavy_set_reduction::is_heavy_set: v and u must be not adjacent");
    assert(!is_reduced(v, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: v must not be reduced");
    assert(!is_reduced(u, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: u must not be reduced");

    for (auto n : v_neighbors_vec)
    {
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }
    for (auto n : u_neighbors_vec)
    {
        if (graph_nodes_set.get(n))
            continue; // only add nodes once
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }

    std::vector<NodeWeight> MWIS_weights(4, 0);
    // assert(graph_nodes.size() > 2 && "ERROR: heavy_set_reduction::is_heavy_set: graph_nodes must have at least 3 nodes");
    MWIS_weights[v_combination::uv] = weights[u] + weights[v];
    // compute MWIS in N(v) + N(u) + N(w):
    if (!solve_induced_subgraph_from_set(MWIS_weights[v_combination::uv], MWIS_weights[v_combination::oo], br_alg, graph_nodes, graph_nodes_set))
        return false;

    if (std::min(weights[v], weights[u]) >= MWIS_weights[v_combination::oo])
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        return true;
    }
    else if (MWIS_weights[v_combination::uv] < MWIS_weights[v_combination::oo])
    {
        return false;
    }
    bool original_heavy_set = false;
    NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
    // compute MWIS[uo] in N(v)\N(u): (included u)
    unset_weights(solver, u_neighbors_vec);
    tiny_solver_clear_solution(solver);
    tiny_solver_solve(solver, config.subgraph_time_limit, no_limit);
    if (solver->time_limit_exceeded || solver->node_limit_exceeded)
        return false;
    MWIS_weights[v_combination::uo] = solver->independent_set_weight;
    if (original_heavy_set && MWIS_weights[v_combination::uo] > weights[v])
        return false;
    MWIS_weights[v_combination::uo] += weights[u];

    // compute MWIS[ov] in N(u)\N(v): (included v)
    set_weights(solver, u_neighbors_vec, weights);
    unset_weights(solver, v_neighbors_vec);
    tiny_solver_clear_solution(solver);
    tiny_solver_solve(solver, config.subgraph_time_limit, no_limit);
    if (solver->time_limit_exceeded || solver->node_limit_exceeded)
        return false;
    MWIS_weights[v_combination::ov] = solver->independent_set_weight;
    if (original_heavy_set && MWIS_weights[v_combination::ov] > weights[u])
        return false;
    MWIS_weights[v_combination::ov] += weights[v];

    if (original_heavy_set)
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
    }
    else
    {
        // get max weights by weight and if equal by participating vertices in the permutation
        v_combination max_weight_combination = static_cast<v_combination>(std::max_element(MWIS_weights.begin(), MWIS_weights.end()) - MWIS_weights.begin());
        NodeWeight best_weight_excluding = 0;
        NodeWeight worst_weight_including = 0;
        switch (max_weight_combination)
        {
        case v_combination::uv:
            br_alg->set(v, IS_status::included);
            br_alg->set(u, IS_status::included);
            break;
        case v_combination::uo:
            best_weight_excluding = std::max({MWIS_weights[v_combination::ov], MWIS_weights[v_combination::oo]});
            worst_weight_including = std::min({MWIS_weights[v_combination::uv], MWIS_weights[v_combination::uo]});
            if (worst_weight_including >= best_weight_excluding)
                br_alg->set(u, IS_status::included);
            break;
        case v_combination::ov:
            best_weight_excluding = std::max({MWIS_weights[v_combination::uo], MWIS_weights[v_combination::oo]});
            worst_weight_including = std::min({MWIS_weights[v_combination::uv], MWIS_weights[v_combination::ov]});
            if (worst_weight_including >= best_weight_excluding)
                br_alg->set(v, IS_status::included);
            break;
        default:
            break;
        }
    }

    return oldn != status.remaining_nodes;
}
void heavy_set_reduction::unset_weights(tiny_solver *solver, std::vector<NodeID> &nodes)
{
    for (NodeID n : nodes)
    {
        solver->subgraph_W[solver->forward_map[n]] = 0;
    }
}
void heavy_set_reduction::set_weights(tiny_solver *solver, std::vector<NodeID> &nodes, std::vector<NodeWeight> &weights)
{
    for (NodeID n : nodes)
    {
        solver->subgraph_W[solver->forward_map[n]] = weights[n];
    }
}
inline void heavy_set_reduction::generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index)
{
    auto &status = br_alg->status;
    auto &config = br_alg->config;
    auto solver = br_alg->subgraph_solver;
    auto &v_neighbors_set = br_alg->set_1;
    for (NodeID common_neighbor = 0; common_neighbor < status.n; common_neighbor++)
    {
        if (br_alg->deg(common_neighbor) < 2)
            continue;
        // check if common neighbor is suitable
        if (std::all_of(status.graph[common_neighbor].begin(), status.graph[common_neighbor].end(), [&](NodeID neighbor)
                        { return br_alg->deg(neighbor) > config.subgraph_node_limit; }))
        {
            continue;
        }
        auto candidates = status.graph[common_neighbor];

        for (size_t v_idx = 0; v_idx < candidates.size() - 1; v_idx++)
        {
            NodeID v = candidates[v_idx];
            if (br_alg->deg(v) > config.subgraph_node_limit)
            {
                continue; // too many neighbors
            }
            get_neighborhood_set(v, br_alg, v_neighbors_set);

            // find second heavy vertex u (not adjacent to v)
            for (size_t u_idx = v_idx + 1; u_idx < candidates.size(); u_idx++)
            {
                NodeID u = candidates[u_idx];
                assert(u != v);
                if (v_neighbors_set.get(u))
                    continue; // look for non adjacent nodes
                if (br_alg->deg(u) + br_alg->deg(v) > config.subgraph_node_limit)
                {
                    continue; // subgraph too large
                }

                if (is_heavy_set(v, v_neighbors_set, u, br_alg))
                {
                    if (is_reduced(v, br_alg))
                        reduction_data[reduction_index][v] = 1;
                    if (is_reduced(u, br_alg))
                        reduction_data[reduction_index][u] = 1;
                    while (status.modified_stack.size() > 0)
                    {
                        NodeID node = status.modified_stack.back();
                        status.modified_stack.pop_back();
                        br_alg->unset(node);
                    }
                }
                else
                {
                    if (solver->time_limit_exceeded)
                    {
                        if (reduction_data[reduction_index][u] != 1)
                            reduction_data[reduction_index][u] = 2;
                        if (reduction_data[reduction_index][v] != 1)
                            reduction_data[reduction_index][v] = 2;
                    }
                }
            }
        }
    }
}
