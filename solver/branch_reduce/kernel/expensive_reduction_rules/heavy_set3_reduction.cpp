
#include "heavy_set3_reduction.h"

#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool heavy_set3_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_heavy_set) return false;
    if (br_alg->heuristically_reducing) return false;
    if (br_alg->t.elapsed() > br_alg->config.time_limit) return false;

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        if (br_alg->t.elapsed() > br_alg->config.time_limit) return;
        reduce_vertex(br_alg, v);
    });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool heavy_set3_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID common_neighbor)
{
    auto &config = br_alg->config;
    if (config.disable_heavy_set)
        return false;
    if (br_alg->blowing_up)
        return false;
    if (br_alg->deg(common_neighbor) < 3)
        return false; // no heavy_set of 3 vertives possible

    auto &status = br_alg->status;
    auto &weights = status.weights;
    assert(status.node_status[common_neighbor] == IS_status::not_set && "ERROR: heavy_set3_reduction::reduce_vertex: node must be unset");
    // check if common neighbor is suitable
    if (std::all_of(status.graph[common_neighbor].begin(), status.graph[common_neighbor].end(), [&](NodeID neighbor)
                    { return br_alg->deg(neighbor) > config.subgraph_node_limit; }))
        return false;
    size_t oldn = status.remaining_nodes;
    auto &v_neighbors_set = br_alg->set_1;
    auto &u_neighbors_set = br_alg->set_2;
    auto &candidates = status.graph[common_neighbor];
    // print all neighbors of common neighbor
    for (size_t v_idx = 0; v_idx < candidates.size() - 2; v_idx++)
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
        for (size_t u_idx = v_idx + 1; u_idx < candidates.size() - 1; u_idx++)
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
            get_neighborhood_set(u, br_alg, u_neighbors_set);

            // find third heavy vertex (not adjacent to  u and v)
            for (size_t w_idx = u_idx + 1; w_idx < candidates.size(); w_idx++)
            {
                NodeID w = candidates[w_idx];
                assert(u != w);
                assert(v != w);
                if (v_neighbors_set.get(w))
                    continue; // look for non adjacent nodes
                if (u_neighbors_set.get(w))
                    continue; // look for non adjacent nodes
                if (is_reduced(w, br_alg))
                    continue;
                if (br_alg->deg(w) + deg_v + br_alg->deg(u) > config.subgraph_node_limit)
                    continue; // subgraph too large
                if (br_alg->deg(w) < 3)
                    continue; // use other reduction
                if (weights[w] + weights[u] + weights[v] < weights[common_neighbor])
                    continue;

                if (is_heavy_set(v, v_neighbors_set, u, u_neighbors_set, w, br_alg))
                {
                    // reduction was applied
                    return oldn != status.remaining_nodes;
                }
            }
        }
    }

    return oldn != status.remaining_nodes;
} 
bool heavy_set3_reduction::is_heavy_set(NodeID v, fast_set& v_neighbors_set, NodeID u, fast_set& u_neighbors_set, NodeID w, branch_and_reduce_algorithm* br_alg) {
    auto& config = br_alg->config;
    auto& status = br_alg->status;
    auto& weights = status.weights;
	auto& v_neighbors_vec = br_alg->buffers[0];
	auto& u_neighbors_vec = br_alg->buffers[1];
    auto& graph_nodes = br_alg->buffers[2];
    auto& graph_nodes_set = br_alg->double_set;
    auto& reverse_mapping = br_alg->buffers[3];
    reverse_mapping.resize(status.n);
    graph_nodes.clear();
    graph_nodes_set.clear();
    graph_access subgraph;

    std::vector<NodeID> w_neighbors_vec(status.graph[w].size());
    get_neighborhood_vector(v, br_alg, v_neighbors_vec);
    get_neighborhood_vector(u, br_alg, u_neighbors_vec);
    get_neighborhood_vector(w, br_alg, w_neighbors_vec);

    assert(!status.graph.adjacent(v, u) && "ERROR: heavy_set_reduction::is_heavy_set: v and u must be not adjacent");
    assert(!status.graph.adjacent(v, w) && "ERROR: heavy_set_reduction::is_heavy_set: v and w must be not adjacent");
    assert(!status.graph.adjacent(u, w) && "ERROR: heavy_set_reduction::is_heavy_set: u and w must be not adjacent");
    assert(!is_reduced(v, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: v must not be reduced");
    assert(!is_reduced(u, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: u must not be reduced");
    assert(!is_reduced(w, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: w must not be reduced");

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
    for (auto n : w_neighbors_vec)
    {
        if (graph_nodes_set.get(n))
            continue; // only add nodes once
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }

    assert(graph_nodes.size() > 2 && "ERROR: heavy_set_reduction::is_heavy_set: graph_nodes must have at least 3 nodes");

    std::vector<NodeWeight> MWIS_weights(8, 0);
    NodeWeight max_weight_limit = std::numeric_limits<NodeWeight>::max();
    // compute MWIS in N(v) + N(u) + N(w):
    MWIS_weights[v_combination::uvw] = weights[u] + weights[v] + weights[w];
    if (!solve_induced_subgraph_from_set(MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::ooo], subgraph, br_alg, graph_nodes, graph_nodes_set, reverse_mapping))
        return false;

    if (std::min(std::min(weights[v], weights[u]), weights[w]) >= MWIS_weights[v_combination::ooo])
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        br_alg->set(w, IS_status::included);
        return true;
    }
    else if (MWIS_weights[v_combination::uvw] < MWIS_weights[v_combination::ooo])
    {
        return false;
    }
    bool original_heavy_set = false;
    // compute MWIS[2] in N(v)+N(w)\N(u): (included u)
    unset_weights(subgraph, u_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uoo], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::uoo] > weights[v] + weights[w])
        return false;
    MWIS_weights[v_combination::uoo] += weights[u];

    // compute MWIS[3] in N(v)\N(u)+N(w): (included u and w)
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uow], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::uow] > weights[v])
        return false;
    MWIS_weights[v_combination::uow] += weights[u] + weights[w];

    // compute MWIS[4] in N(v) + N(u) \ N(w): (included w )
    set_weights(subgraph, u_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::oow], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::oow] > weights[u] + weights[v])
        return false;
    MWIS_weights[v_combination::oow] += weights[w];

    // compute MWIS[5] in N(u)+N(w)\N(v): (included v)
    set_weights(subgraph, w_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, v_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::ovo], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::ovo] > weights[u] + weights[w])
        return false;
    MWIS_weights[v_combination::ovo] += weights[v];

    // compute MWIS in N(u)\N(v)+N(w): (included v and w)
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::ovw], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::ovw] > weights[u])
        return false;
    MWIS_weights[v_combination::ovw] += weights[v] + weights[w];

    // compute MWIS in N(w)\N(v)+N(u):
    set_weights(subgraph, w_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, v_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uvo], subgraph, config, max_weight_limit))
        return false;
    if (original_heavy_set && MWIS_weights[v_combination::uvo] > weights[w])
        return false;
    MWIS_weights[v_combination::uvo] += weights[v] + weights[u];
    size_t oldn = status.remaining_nodes;

    if (original_heavy_set)
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        br_alg->set(w, IS_status::included);
    }
    else
    {
        // get max weights by weight and if equal by participating vertices in the permutation
        v_combination max_weight_combination = static_cast<v_combination>(std::max_element(MWIS_weights.begin(), MWIS_weights.end()) - MWIS_weights.begin());
        switch (max_weight_combination)
        {
        case v_combination::uvw:
            br_alg->set(v, IS_status::included);
            br_alg->set(u, IS_status::included);
            br_alg->set(w, IS_status::included);
            break;
        case v_combination::uvo:
            // check all other combinations including vertex to be better than any excluding it
            if (check_u_combination(MWIS_weights))
                br_alg->set(u, IS_status::included);
            if (check_v_combination(MWIS_weights))
                br_alg->set(v, IS_status::included);
            break;
        case v_combination::uow:
            if (check_u_combination(MWIS_weights))
                br_alg->set(u, IS_status::included);
            if (check_w_combination(MWIS_weights))
                br_alg->set(w, IS_status::included);
            break;
        case v_combination::ovw:
            if (check_v_combination(MWIS_weights))
                br_alg->set(v, IS_status::included);
            if (check_w_combination(MWIS_weights))
                br_alg->set(w, IS_status::included);
            break;
        case v_combination::uoo:
            if (check_u_combination(MWIS_weights))
                br_alg->set(u, IS_status::included);
            break;
        case v_combination::ovo:
            if (check_v_combination(MWIS_weights))
                br_alg->set(v, IS_status::included);
            break;
        case v_combination::oow:
            if (check_w_combination(MWIS_weights))
                br_alg->set(w, IS_status::included);
            break;
        default:
            break;
        }
    }

    return oldn != status.remaining_nodes;
}
void heavy_set3_reduction::unset_weights(graph_access& graph, std::vector<NodeID>& nodes, std::vector<NodeID>& reverse_mapping) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], 0);
    }
}
void heavy_set3_reduction::set_weights(graph_access& graph, std::vector<NodeID>& nodes, std::vector<NodeID>& reverse_mapping, std::vector<NodeWeight>& weights) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], weights[n]);
    }
}
bool heavy_set3_reduction::check_u_combination(std::vector<NodeWeight> &MWIS_weights)
{
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::oow], MWIS_weights[v_combination::ovo], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::uow], MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::uoo]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight)
                       { return weight >= best_weight_excluding; });
}
bool heavy_set3_reduction::check_v_combination(std::vector<NodeWeight> &MWIS_weights)
{
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::uow], MWIS_weights[v_combination::uoo], MWIS_weights[v_combination::oow], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::ovo]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight)
                       { return weight >= best_weight_excluding; });
}
bool heavy_set3_reduction::check_w_combination(std::vector<NodeWeight> &MWIS_weights)
{
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::uoo], MWIS_weights[v_combination::ovo], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::uow], MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::oow]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight)
                       { return weight >= best_weight_excluding; });
}
