
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    size_t oldn = br_alg->status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&br_alg, &status](NodeID v)
                            {
                                NodeWeight neighbor_weights = 0;
                                for (NodeID u : status.graph[v])
                                {
                                    neighbor_weights += status.weights[u];
                                    if (status.weights[v] < neighbor_weights)
                                    {
                                        return;
                                    }
                                }
                                br_alg->set(v, IS_status::included);
                            });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != status.remaining_nodes) std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;
    return oldn != br_alg->status.remaining_nodes;
}

inline bool neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_neighborhood) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (status.weights[v] >= neighbor_weights)
    {
        br_alg->set(v, IS_status::included);
    }
    return oldn != br_alg->status.remaining_nodes;
}
bool neighborhood_reduction::is_suited(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    return (status.weights[v] >= static_cast<NodeWeight>(std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&](NodeWeight sum, NodeID neighbor)
                                                                         { 
			if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor]; 
			else return sum; })));
}

bool fold1_reduction::reduce(branch_and_reduce_algorithm *br_alg)
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
inline bool fold1_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
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
void fold1_reduction::fold(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes)
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
void fold1_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.nodes.deg1_node);

    status.weights[data.nodes.fold_node] += data.deg1_weight;
    status.reduction_offset -= data.deg1_weight;

    restore_vec.pop_back();
}
void fold1_reduction::apply(branch_and_reduce_algorithm *br_alg)
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
bool fold1_reduction::is_suited(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    return br_alg->deg(v) <= 1;
}

bool fold2_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_fold2) return false;

    auto &status = br_alg->status;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    size_t oldn = status.remaining_nodes;
    size_t oldw = status.reduction_offset;

	for_each_changed_vertex(br_alg, [this, &br_alg, &status](NodeID v) 
                            { reduce_vertex(br_alg, v); });
    
    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
	// if (oldn != status.remaining_nodes) std::cout << "fold2 redu -> " << (oldn - status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes || oldw != status.reduction_offset;
}
inline bool fold2_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_fold2) return false;
    if (br_alg->deg(v) != 2)
        return false;

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    size_t oldw = status.reduction_offset;

    // set bigger and smaller neighbor
    NodeID bigger = status.graph[v][0];
    NodeID smaller = status.graph[v][1];
    if (status.weights[v] >= status.weights[bigger] + status.weights[smaller])
    {
        br_alg->set(v, IS_status::included);
        return true;
    }
    if (status.weights[bigger] < status.weights[smaller])
    {
        smaller = status.graph[v][0];
        bigger = status.graph[v][1];
    }

    // check triangle condition
    bool triangle = false;
    for (NodeID neighbor : status.graph[smaller])
    {
        if (neighbor == bigger)
        {
            triangle = true;
            break;
        }
    }

    if (triangle && status.weights[bigger] <= status.weights[v])
    {
        br_alg->set(v, IS_status::included);
    }
    else if (triangle && status.weights[bigger] > status.weights[v] && status.weights[smaller] <= status.weights[v])
    {
        if (br_alg->config.disable_triangle_mid)
            return false;
        this->fold_triangle_mid_weight(br_alg, {v, {bigger, smaller}});
    }
    else if (triangle)
    {
        if (br_alg->config.disable_triangle_min)
            return false;
        this->fold_triangle_min_weight(br_alg, {v, {bigger, smaller}});
    }
    else if (status.weights[v] >= status.weights[bigger])
    {
        if (br_alg->config.disable_v_shape_max)
            return false;
        this->fold_v_shape_max_weight(br_alg, {v, {status.graph[v][0], status.graph[v][1]}});
    }
    else if (status.weights[v] >= status.weights[smaller])
    {
        if (br_alg->config.disable_v_shape_mid)
            return false;
        this->fold_v_shape_mid_weight(br_alg, {v, {bigger, smaller}});
    }
    else
    {
        assert(status.weights[v] < status.weights[smaller] && "v is not the smallest");
        if (br_alg->config.disable_v_shape_min)
            return false;
        // if (br_alg->blowing_up) return;
        this->fold_v_shape_min_weight(br_alg, {v, {bigger, smaller}});
    }

	return oldn != status.remaining_nodes || oldw != status.reduction_offset;
}
void fold2_reduction::fold_triangle_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes)
{
    auto &status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);
    br_alg->set(nodes.neighbors[1], IS_status::excluded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], fold2_reduction::fold2_case::triangle_mid, {}, {}});

    status.reduction_offset += status.weights[nodes.deg2_node];
    status.weights[nodes.neighbors[0]] -= status.weights[nodes.deg2_node];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.neighbors[0]);
    br_alg->add_next_level_neighborhood(nodes.neighbors[0]);
}
void fold2_reduction::fold_triangle_min_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes)
{
    auto &status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], fold2_reduction::fold2_case::triangle_min, {}, {}});

    status.reduction_offset += status.weights[nodes.deg2_node];
    status.weights[nodes.neighbors[0]] -= status.weights[nodes.deg2_node];
    status.weights[nodes.neighbors[1]] -= status.weights[nodes.deg2_node];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.neighbors[0]);
    br_alg->add_next_level_neighborhood(nodes.neighbors[0]);
}
void fold2_reduction::fold_v_shape_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes)
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->set_1;
    NodeID bigger = nodes.neighbors[0];
    NodeID smaller = nodes.neighbors[1];
    NodeWeight deg2_weight = status.weights[nodes.deg2_node];

    restore_vec.push_back({nodes, deg2_weight, fold2_reduction::fold2_case::v_shape_mid, {}, {}});
    br_alg->set(nodes.deg2_node, IS_status::folded);

    status.reduction_offset += deg2_weight;
    status.weights[bigger] -= deg2_weight;
    neighbors.clear();
    neighbors.add(nodes.deg2_node);
    neighbors.add(smaller);
    neighbors.add(bigger);

    for (auto neighbor : status.graph[smaller])
    {
        neighbors.add(neighbor);
    }

    for (auto neighbor : status.graph[bigger])
    {
        if (neighbors.add(neighbor))
        {
            status.graph.add_edge_undirected(smaller, neighbor);
            restore_vec.back().node_vecs[1].push_back(neighbor);
        }
    }

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(smaller);
    br_alg->add_next_level_node(bigger);
    br_alg->add_next_level_neighborhood(smaller);
}
void fold2_reduction::fold_v_shape_max_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes)
{
    auto &status = br_alg->status;
    auto &neighbors = br_alg->set_1;

    br_alg->set(nodes.neighbors[1], IS_status::folded, false);
    br_alg->set(nodes.neighbors[0], IS_status::folded, true);

	restore_vec.push_back({ nodes, status.weights[nodes.deg2_node], fold2_reduction::fold2_case::v_shape_max, status.graph[nodes.deg2_node], {} });

    status.reduction_offset += status.weights[nodes.deg2_node];
    status.weights[nodes.deg2_node] = status.weights[nodes.neighbors[0]] + status.weights[nodes.neighbors[1]] - status.weights[nodes.deg2_node];

    std::vector<NodeID> new_neighbors;
    neighbors.clear();
    neighbors.add(nodes.deg2_node);

    for (size_t i = 0; i < 2; i++)
    {
        for (auto neighbor : status.graph[nodes.neighbors[i]])
        {
            if (neighbors.add(neighbor))
            {
                new_neighbors.push_back(neighbor);
                status.graph.add_edge_directed(neighbor, nodes.deg2_node);
                restore_vec.back().node_vecs[i].push_back(neighbor);
            }
        }
    }

    status.graph[nodes.deg2_node] = dynamic_graph::neighbor_list(std::move(new_neighbors));
    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.deg2_node);
    br_alg->add_next_level_neighborhood(nodes.deg2_node);
}
void fold2_reduction::fold_v_shape_min_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;

    status.modified_stack.push_back(br_alg->MODIFIED_TOKEN);

    NodeWeight deg2_weight = status.weights[nodes.deg2_node];
    restore_vec.push_back({nodes, deg2_weight, fold2_reduction::fold2_case::v_shape_min, {}});
    status.reduction_offset += deg2_weight;
    status.weights[nodes.neighbors[0]] -= deg2_weight;
    status.weights[nodes.neighbors[1]] -= deg2_weight;

    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[0]);
    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[1]);

    neighbors.clear();
    neighbors.add(nodes.deg2_node);

    for (size_t i = 0; i < 2; i++)
    {
        for (auto neighbor : status.graph[nodes.neighbors[i]])
        {
            if (status.node_status[neighbor] != IS_status::not_set)
                continue;
            if (neighbors.add(neighbor))
            {
                status.graph.add_edge_undirected(neighbor, nodes.deg2_node);
                restore_vec.back().node_vecs[i].push_back(neighbor);
            }
        }
    }

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.deg2_node);
    br_alg->add_next_level_node(nodes.neighbors[0]);
    br_alg->add_next_level_node(nodes.neighbors[1]);
    br_alg->add_next_level_neighborhood(nodes.neighbors[0]);
    br_alg->add_next_level_neighborhood(nodes.neighbors[1]);
}
void fold2_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();
    auto &deg2_node = data.nodes.deg2_node;
    auto &nb0 = data.nodes.neighbors[0];
    auto &nb1 = data.nodes.neighbors[1];

    status.reduction_offset -= data.deg2_weight;
    switch (data.fold_case)
    {
        case fold2_reduction::fold2_case::triangle_mid:
            br_alg->unset(deg2_node);
	        status.weights[nb0] += data.deg2_weight;
            break;

        case fold2_reduction::fold2_case::triangle_min:
            br_alg->unset(deg2_node);
            status.weights[nb0] += data.deg2_weight;
            status.weights[nb1] += data.deg2_weight;
            break;

        case fold2_reduction::fold2_case::v_shape_max:
            status.graph.hide_node(data.nodes.deg2_node);
            status.graph[data.nodes.deg2_node] = std::move(data.neighbor_list);
            for (int i = 0; i < 2; i++) {
                br_alg->unset(data.nodes.neighbors[i]);

                for (NodeID neighbor : data.node_vecs[i]) {
                    status.graph.replace_last_restored_edge(neighbor, data.nodes.neighbors[i]);
                }
            }
            status.weights[deg2_node] = data.deg2_weight;
            break;

        case fold2_reduction::fold2_case::v_shape_mid:
            br_alg->unset(deg2_node);
            br_alg->unset(nb0, false);
            br_alg->unset(nb1, false);

            // restore neighbors of smaller node
            for (auto neighbor : data.node_vecs[1]) {
                status.graph.hide_edge_undirected(neighbor, nb1);
            }

            status.weights[nb0] += data.deg2_weight;
            break;

        case fold2_reduction::fold2_case::v_shape_min:
            br_alg->unset(deg2_node, false);
            br_alg->unset(nb0, false);
            br_alg->unset(nb1, false);

	        for (size_t i = 0; i < 2; i++) {
                for (auto second_neighbor : data.node_vecs[i]) {
                    status.graph.hide_edge_undirected(second_neighbor, deg2_node);
                }
                status.graph.add_edge_undirected(data.nodes.neighbors[i], deg2_node);
                status.weights[data.nodes.neighbors[i]] += data.deg2_weight;
	        }  
            status.weights[deg2_node] = data.deg2_weight;
            break;
    }
    restore_vec.pop_back();
}
void fold2_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto nodes = restore_vec.back().nodes;
    int fold_case = restore_vec.back().fold_case;

    auto deg_2_status = status.node_status[nodes.deg2_node];
    auto bigger_status = status.node_status[nodes.neighbors[0]];
    auto smaller_status = status.node_status[nodes.neighbors[1]];

	restore(br_alg);

    status.node_status[nodes.deg2_node]    = IS_status::excluded;
    status.node_status[nodes.neighbors[0]] = IS_status::excluded;
    status.node_status[nodes.neighbors[1]] = IS_status::excluded;

    switch (fold_case)
    {
        case fold2_reduction::fold2_case::v_shape_max:
            if (deg_2_status == IS_status::included) {
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.node_status[nodes.neighbors[1]] = IS_status::included;
            } else {
                status.node_status[nodes.deg2_node]    = IS_status::included;
            }
            status.is_weight += status.weights[nodes.deg2_node];
            break;
        case fold2_reduction::fold2_case::v_shape_mid:

            if (smaller_status == IS_status::included) { 
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.node_status[nodes.neighbors[1]] = IS_status::included;
                status.is_weight += status.weights[nodes.neighbors[0]];
                status.is_weight += status.weights[nodes.neighbors[1]];
            } else if (bigger_status == IS_status::included) {
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.is_weight += status.weights[nodes.neighbors[0]];
            } else {
                status.node_status[nodes.deg2_node] = IS_status::included;
                status.is_weight += status.weights[nodes.deg2_node];
            }
            break;
        case fold2_reduction::fold2_case::v_shape_min:

            if (smaller_status == IS_status::included) {
           	    status.node_status[nodes.neighbors[1]] = IS_status::included;
          		status.is_weight += status.weights[nodes.neighbors[1]];
            }
            if (bigger_status == IS_status::included) {
           	    status.node_status[nodes.neighbors[0]] = IS_status::included;
          		status.is_weight += status.weights[nodes.neighbors[0]];
            }
            if (smaller_status == IS_status::excluded && bigger_status == IS_status::excluded) {
	            status.node_status[nodes.deg2_node]    = IS_status::included;
	            status.is_weight += status.weights[nodes.deg2_node];
            }
            break;

        default: // triangle_mid and triangle_min
            if (bigger_status == IS_status::included) {
            	status.node_status[nodes.neighbors[0]] = IS_status::included;
            	status.is_weight += status.weights[nodes.deg2_node];

            } else if (smaller_status == IS_status::included){
            	status.node_status[nodes.neighbors[1]] = IS_status::included;
            	status.is_weight += status.weights[nodes.deg2_node];

            } else {
            	status.node_status[nodes.deg2_node] = IS_status::included;
            	status.is_weight += status.weights[nodes.deg2_node];
            }
            break;
    }
}
bool fold2_reduction::is_suited(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    return br_alg->deg(v) == 2;
}

bool single_edge_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_basic_se) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_basic_se) return false;
    if (br_alg->deg(v) == 0)
        return false;

    auto &status = br_alg->status;
    auto &graph = br_alg->status.graph;
    auto &weights = br_alg->status.weights;
    auto &neighbors = br_alg->set_1;
    size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, neighbors);

    for (NodeID neighbor : graph[v])
    {
        if (is_reduced(neighbor, br_alg))
            continue;
        if (weights[v] <= weights[neighbor])
        { // otherwise not applicable to this edge
            NodeWeight partial_neighbor_sum = 0;
            for (NodeID second_neighbor : status.graph[neighbor])
            {
                if (!is_reduced(second_neighbor, br_alg) && !neighbors.get(second_neighbor))
                {
                    partial_neighbor_sum += status.weights[second_neighbor];
                    if (partial_neighbor_sum > weights[neighbor])
                        break;
                }
            }

// note: weight of v is in partial_neighbor_sum included
// if N(neighbor) \subset N(v) partial_neighbor_sum = weights[v]
#ifdef gen_training_data
            if (partial_neighbor_sum == weights[v])
                break; // should be domination
#endif
            if (partial_neighbor_sum <= weights[neighbor])
            {
                br_alg->set(v, IS_status::excluded);
                break;
            }
        }
    }

    return oldn != status.remaining_nodes;
}

bool extended_single_edge_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_extended_se) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}
inline bool extended_single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_extended_se) return false;
    if (br_alg->deg(v) == 0)
        return false;

	auto& status = br_alg->status;
    auto& weights = status.weights;
	size_t oldn = status.remaining_nodes;
	auto& neighbors = br_alg->set_1;
	auto& neighbors_vec = br_alg->buffers[0];


    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    NodeID max_weight_neighbor = get_max_weight_neighbor(v, br_alg);
    assert(max_weight_neighbor != v && "ERROR: max_weight_neighbor == v");
    if (weights[v] < neighbors_weight - weights[max_weight_neighbor])
    {
        return false;
    }
#ifndef gen_training_data
    if (try_neighborhood_reduction(v, br_alg, neighbors_weight))
    {
        return oldn != br_alg->status.remaining_nodes;
    }
#endif
    get_neighborhood_set(v, br_alg, neighbors);
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    std::sort(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b)
              { return weights[a] > weights[b]; });
    for (NodeID max_neighbor : neighbors_vec)
    {
#ifndef gen_training_data
        if (v > max_neighbor)
            continue;
#endif
        if (weights[v] < neighbors_weight - weights[max_neighbor])
            break;

        bool progress = false;
        for (NodeID neighbor : status.graph[max_neighbor])
        {
            if (neighbor == v)
                continue;
            if (is_reduced(neighbor, br_alg))
                continue;
            // exclude neighborhood intersection and update neighborhood
            if (neighbors.get(neighbor))
            {
                br_alg->set(neighbor, IS_status::excluded);
                progress = true;
            }
        }

        if (progress)
            break;
    }

    return oldn != status.remaining_nodes;
}

bool domination_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    if (br_alg->config.disable_domination)
        return false;
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
	// if (oldn != br_alg->status.remaining_nodes) std::cout << "domination redu -> " << (oldn - br_alg->status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes;
}
inline bool domination_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    if (br_alg->config.disable_domination)
        return false;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->set_1;
    size_t oldn = status.remaining_nodes;
#ifndef gen_training_data
    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (try_neighborhood_reduction(v, br_alg, neighbor_weights))
    {
        return oldn != status.remaining_nodes;
    }
#endif

    size_t neighbors_count = 0;
    bool is_subset;

    neighbors.clear();
    neighbors.add(v);
    for (NodeID neighbor : status.graph[v])
    {
        neighbors.add(neighbor);
        neighbors_count++;
    }

    for (NodeID neighbor : status.graph[v])
    {
        if (br_alg->deg(neighbor) > neighbors_count)
            continue;
        if (status.weights[neighbor] < status.weights[v])
            continue;

        is_subset = true;

        for (NodeID neighbor2 : status.graph[neighbor])
        {
            if (!neighbors.get(neighbor2))
            {
                is_subset = false;
                break;
            }
        }

        // if (is_subset && status.weights[neighbor] >= status.weights[v]) {
        if (is_subset)
        {
            br_alg->set(v, IS_status::excluded);
            break;
            }
        }

	return oldn != status.remaining_nodes;
}

bool extended_domination_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if(br_alg->config.disable_extended_domination) return false;
    if(br_alg->blowing_up) return false;
	auto& status = br_alg->status;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif
	size_t oldn = status.remaining_nodes;
    bool progress = false;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        progress = reduce_vertex(br_alg, v);
    });

    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
	return progress;
}
inline bool extended_domination_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if(br_alg->config.disable_domination) return false;
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
	#ifndef generate_training_data
        NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
        if (try_neighborhood_reduction(v, br_alg, neighbor_weights)) {
            return oldn != status.remaining_nodes;
        }
    #endif

    size_t neighbors_count = status.graph[v].size();
    bool is_subset;
    bool progress = false;

    neighbors.clear();
    neighbors.add(v);
    for (NodeID neighbor : status.graph[v]) {
        neighbors.add(neighbor);
    }

    for (NodeID neighbor : status.graph[v]) {
        if (br_alg->deg(neighbor) > neighbors_count)
            continue;

        is_subset = true;

        for (NodeID neighbor2 : status.graph[neighbor]) {
            if (!neighbors.get(neighbor2)) {
                is_subset = false;
                break;
            }
        }

        if (is_subset) {
            progress = true;
            if (status.weights[neighbor] < status.weights[v])
            {
                // // remove edge
                fold(br_alg, v, neighbor);
                neighbors.remove(neighbor);
                neighbors_count--;
            }
            else 
            {
                br_alg->set(v, IS_status::excluded);
                break;
            }
            assert(neighbors.get(neighbor) == false && "ERROR: neighbor not removed");
        }
    }
    
	return progress;
}
void extended_domination_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID v, NodeID neighbor) {

    auto& status = br_alg->status;

    status.modified_stack.push_back(br_alg->MODIFIED_TOKEN);
	restore_vec.push_back({v, neighbor, status.weights[neighbor]});

    status.weights[v] -= status.weights[neighbor];
    status.graph.hide_edge_undirected(v, neighbor);

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(v);
    br_alg->add_next_level_node(neighbor);
    br_alg->add_next_level_neighborhood(v);
}
void extended_domination_reduction::restore(branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    auto& data = restore_vec.back();

    status.weights[data.v] += data.neighborWeight;
    status.graph.add_edge_undirected(data.v, data.neighbor);
    restore_vec.pop_back();
}
void extended_domination_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
    assert(restore_vec.size() > 0 && "restore vec empty");
	auto v = restore_vec.back().v;
	auto neighbor = restore_vec.back().neighbor;
	auto neighbor_status = status.node_status[neighbor];
	auto v_status = status.node_status[v];
    if (v_status == IS_status::included ) 
        assert(neighbor_status == IS_status::included && "ERROR: domination_reduction::apply: v included neighbor not included");
	restore(br_alg);

    status.node_status[neighbor] = neighbor_status;
    status.node_status[v] = v_status;
	if (v_status == IS_status::included) {
	    status.is_weight += status.weights[neighbor];
        if (neighbor_status == IS_status::included) {
            status.is_weight -= status.weights[neighbor];
		    status.node_status[neighbor] = IS_status::excluded;
        }
    }
}

bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_clique_neighborhood_fast) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

	// if (oldn != status.remaining_nodes) std::cout << "clique_neighborhood fast redu -> " << (oldn - status.remaining_nodes) << std::endl;
	return oldn != status.remaining_nodes;
}
inline bool clique_neighborhood_reduction_fast::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
// if (br_alg->config.disable_clique_neighborhood_fast) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    auto &neighbors = br_alg->buffers[0];
    auto &neighborhood = br_alg->set_1;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (this->try_neighborhood_reduction(v, br_alg, neighbor_weights))
    {
        return oldn != br_alg->status.remaining_nodes;
    }
    this->get_neighborhood_vector(v, br_alg, neighbors);
    this->get_neighborhood_set(v, br_alg, neighborhood);

    std::sort(neighbors.begin(), neighbors.end(), [&status](const NodeID lhs, const NodeID rhs)
              { return status.weights[lhs] > status.weights[rhs]; });

    bool is_reducible = false;

    for (size_t i = 0; i < neighbors.size() && !is_reducible; i++)
    {
        NodeID neighbor1 = neighbors[i];

        if (!neighborhood.get(neighbor1))
            continue;

        for (NodeID neighbor2 : status.graph[neighbor1])
        {
            if (neighbor2 != neighbor1 && neighborhood.get(neighbor2))
            {
                // triangle [v, neighbor1, neighbor2] found
                neighbor_weights -= std::min(status.weights[neighbor1], status.weights[neighbor2]);
                neighborhood.remove(neighbor1);
                neighborhood.remove(neighbor2);

                if (status.weights[v] >= neighbor_weights)
                {
                    is_reducible = true;
                }

                break;
            }
        }
    }

    if (is_reducible)
    {
        br_alg->set(v, IS_status::included);
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_clique_neighborhood) return false;
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
inline bool clique_neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_clique_neighborhood) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    if (partition_into_cliques(v, br_alg))
    {
        br_alg->set(v, IS_status::included);
    }

    return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction::partition_into_cliques(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &weights = br_alg->status.weights;
    auto &neighbors_vec = br_alg->buffers[0];
    auto &clique_neighbors_set = br_alg->set_1;
    target_weight = weights[v];

    neighbor_weights = get_neighborhood_weight(v, br_alg);
#ifdef gen_training_data
    if (neighbor_weights <= target_weight)
        return false;
#endif
#ifndef gen_training_data
    if (neighbor_weights <= target_weight)
        return true;
#endif
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    // partition neigbors of v into cliques
    while (neighbors_vec.size() >= 2 && br_alg->config.reduction_time_limit > br_alg->t.elapsed())
    {
        NodeID max_neighbor = *std::max_element(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b)
                                                { return weights[a] < weights[b]; });
        clique_neighbors_set.clear();
        for (auto neighbor : neighbors_vec)
        {
            clique_neighbors_set.add(neighbor);
        }
        clique_neighbors_set.remove(max_neighbor);

        if (expand_clique(max_neighbor, neighbors_vec, clique_neighbors_set, br_alg))
            return true;
    }

    return false;
}
bool clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, std::vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set, branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& temp_set = br_alg->set_2;
	auto& clique_set = br_alg->double_set;
	clique_set.clear();

    size_t local_max = status.n;
    NodeWeight local_max_weight = 0;
    bool intersection_empty = true;

    while (true)
    {
        // intersect neighbors of clique with neighbors of max_neighbor
        intersection_empty = true;
        local_max_weight = 0;
        clique_set.add(max_neighbor);
        temp_set.clear();

        for (auto neighbor : status.graph[max_neighbor])
        {
            if (clique_neighbors_set.get(neighbor))
            {
                temp_set.add(neighbor);
                intersection_empty = false;

                if (status.weights[neighbor] > local_max_weight)
                {
                    local_max = neighbor;
                    local_max_weight = status.weights[neighbor];
                }
            }
        }
        if (intersection_empty)
            break;
        // if (local_max == status.n) break;

        // add local_max to current clique
        neighbor_weights -= local_max_weight;
        if (neighbor_weights <= target_weight)
            return true;

        std::swap(clique_neighbors_set, temp_set);
        clique_neighbors_set.remove(local_max);
        max_neighbor = local_max;
    }

    auto &reamining_neighbors = br_alg->buffers[1];
    reamining_neighbors.clear();

    // adjust neigbors_vec
    for (auto neighbor : neighbors_vec)
    {
        if (!clique_set.get(neighbor))
            reamining_neighbors.push_back(neighbor);
    }

    std::swap(reamining_neighbors, neighbors_vec);
    return false;
}

bool clique_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_clique) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    auto &set_1 = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    auto &isolated = br_alg->buffers[1];
    std::vector<NodeID> non_isolated;

    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            {
                                neighbors.clear();
                                set_1.clear();
                                set_1.add(v);

                                // find potential clique
                                for (NodeID neighbor : status.graph[v])
                                {
                                    neighbors.push_back(neighbor);
                                    set_1.add(neighbor);
                                }

                                // check if clique
                                isolated.clear();
                                isolated.push_back(v);
                                non_isolated.clear();

                                size_t max_isolated_idx = 0;
                                weighted_node max_isolated{v, status.weights[v]};
                                weighted_node max_non_isolated{0, 0};

                                for (auto neighbor : neighbors)
                                {
                                    size_t count = 0;
                                    bool is_isolated = true;

                                    for (NodeID neighbor_2nd : status.graph[neighbor])
                                    {
                                        if (set_1.get(neighbor_2nd))
                                            count++;
                                        else
                                            is_isolated = false;
                                    }
                                    if (count != neighbors.size())
                                        return;

                                    if (is_isolated)
                                    {
                                        isolated.push_back(neighbor);
                                        if (status.weights[neighbor] > max_isolated.weight)
                                        {
                                            max_isolated = {neighbor, status.weights[neighbor]};
                                            max_isolated_idx = isolated.size() - 1;
                                        }
                                    }
                                    else
                                    {
                                        non_isolated.push_back(neighbor);
                                        if (status.weights[neighbor] > max_non_isolated.weight)
                                        {
                                            max_non_isolated = {neighbor, status.weights[neighbor]};
                                        }
                                    }
                                }

                                // one of "isolated" members has highest weight of clique: Add to IS
                                // also handles completely isolated cliques
                                if (max_isolated.weight >= max_non_isolated.weight)
                                {
                                    br_alg->set(max_isolated.node, IS_status::included);
                                    return;
                                }

                                // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
                                isolated[max_isolated_idx] = isolated.back();
                                isolated.pop_back();

                                for (auto neighbor : isolated)
                                {
                                    br_alg->set(neighbor, IS_status::excluded);
                                }

                                for (size_t i = 0; i < non_isolated.size(); i++)
                                {
                                    NodeID neighbor = non_isolated[i];
                                    if (status.weights[neighbor] <= max_isolated.weight)
                                    {
                                        br_alg->set(neighbor, IS_status::excluded);
                                        non_isolated[i] = non_isolated.back();
                                        non_isolated.pop_back();
                                        i--;
                                    }
                                }

                                fold(br_alg, std::move(max_isolated), std::move(non_isolated)); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != status.remaining_nodes) std::cout << "clique redu -> " << (oldn - status.remaining_nodes) << std::endl;
    return oldn != status.remaining_nodes;
}

inline bool clique_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    auto &status = br_alg->status;
    auto &set_1 = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    auto &isolated = br_alg->buffers[1];
    std::vector<NodeID> non_isolated;
    size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, set_1);
    get_neighborhood_vector(v, br_alg, neighbors);
    set_1.add(v);

    // check if clique
    // non_isolated.clear();
    isolated.clear();
    isolated.push_back(v);

    size_t max_isolated_idx = 0;
    weighted_node max_isolated{v, status.weights[v]};
    weighted_node max_non_isolated{0, 0};

    for (auto neighbor : neighbors)
    {
        size_t count = 0;
        bool is_isolated = true;

        for (NodeID neighbor_2nd : status.graph[neighbor])
        {
            if (set_1.get(neighbor_2nd))
                count++;
            else
                is_isolated = false;
        }

        if (is_isolated)
        {
            isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_isolated.weight)
            {
                max_isolated = {neighbor, status.weights[neighbor]};
                max_isolated_idx = isolated.size() - 1;
            }
        }
        else
        {
            non_isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_non_isolated.weight)
            {
                max_non_isolated = {neighbor, status.weights[neighbor]};
            }
        }

        if (count != neighbors.size())
            return false;
    }

    // one of "isolated" members has highest weight of clique: Add to IS
    // also handles completely isolated cliques
    if (max_isolated.weight >= max_non_isolated.weight)
    {
        br_alg->set(max_isolated.node, IS_status::included);
        return oldn != br_alg->status.remaining_nodes;
    }

    // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
    isolated[max_isolated_idx] = isolated.back();
    isolated.pop_back();

    for (auto neighbor : isolated)
    {
        br_alg->set(neighbor, IS_status::excluded);
    }

    for (size_t i = 0; i < non_isolated.size(); i++)
    {
        NodeID neighbor = non_isolated[i];
        if (status.weights[neighbor] <= max_isolated.weight)
        {
            br_alg->set(neighbor, IS_status::excluded);
            non_isolated[i] = non_isolated.back();
            non_isolated.pop_back();
            i--;
        }
    }

    fold(br_alg, std::move(max_isolated), std::move(non_isolated));

    return oldn != status.remaining_nodes;
}
void clique_reduction::fold(branch_and_reduce_algorithm *br_alg, const weighted_node &isolated, std::vector<NodeID> &&non_isolated)
{
    auto &status = br_alg->status;

    br_alg->set(isolated.node, IS_status::folded);
    status.reduction_offset += isolated.weight;

    for (auto node : non_isolated)
    {
        status.weights[node] -= isolated.weight;
        br_alg->add_next_level_neighborhood(node);
    }

    status.folded_stack.push_back(get_reduction_type());
    br_alg->add_next_level_neighborhood(non_isolated);

    restore_vec.emplace_back(isolated, std::move(non_isolated));
}
void clique_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.isolated.node);
    status.reduction_offset -= data.isolated.weight;

    for (auto node : data.non_isolated)
    {
        status.weights[node] += data.isolated.weight;
    }

    restore_vec.pop_back();
}
void clique_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto isolated = restore_vec.back().isolated.node;

    bool set_isolated = true;

    for (auto node : restore_vec.back().non_isolated)
    {
        if (status.node_status[node] == IS_status::included)
        {
            set_isolated = false;
            break;
        }
    }

    status.is_weight += restore_vec.back().isolated.weight;

    restore(br_alg);

    if (set_isolated)
    {
        status.node_status[isolated] = IS_status::included;
    }
    else
    {
        status.node_status[isolated] = IS_status::excluded;
    }
}

bool funnel_reduction::reduce(branch_and_reduce_algorithm *br_alg)
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
inline bool funnel_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
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
bool funnel_reduction::is_funnel(NodeID node, NodeID& funnel_neighbor, branch_and_reduce_algorithm* br_alg, fast_set& funnel_set, std::vector<NodeID>& funnel_nodes) 
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
bool funnel_reduction::is_clique(branch_and_reduce_algorithm* br_alg, fast_set& clique_set, std::vector<NodeID>& clique_nodes) {
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
void funnel_reduction::fold(const fold_data &data, fast_set &funnel_set, branch_and_reduce_algorithm *br_alg)
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
void funnel_reduction::restore(branch_and_reduce_algorithm *br_alg)
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
void funnel_reduction::apply(branch_and_reduce_algorithm *br_alg)
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

bool funnel_fold_reduction::reduce(branch_and_reduce_algorithm *br_alg)
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
inline bool funnel_fold_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_funnel_fold) return false;
    if (br_alg->deg(v) <= 2)
        return false; // fold1 or 2
    auto &weights = br_alg->status.weights;
    auto &remaining_n = br_alg->status.remaining_nodes;
    auto &funnel_set = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    size_t oldn = remaining_n;
    NodeID funnel_neighbor = br_alg->status.n;
    assert(!is_reduced(v, br_alg) && "ERROR: funnel_reduction::reduce_vertex: node must be unset");

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    /* bool skip = false; */

    if (std::none_of(neighbors.begin(), neighbors.end(), [&](NodeID neighbor)
                     { return weights[neighbor] > weights[v]; }))
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
bool funnel_fold_reduction::is_funnel(NodeID v, NodeID& funnel_neighbor, branch_and_reduce_algorithm* br_alg, fast_set& funnel_set, std::vector<NodeID>& funnel_nodes) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    funnel_neighbor = status.n;
    if (is_clique(br_alg, funnel_set, funnel_nodes))
        return false;

    for (NodeID u : funnel_nodes)
    {
        assert(u != v && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(u, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");
        assert(std::find(funnel_nodes.begin(), funnel_nodes.end(), u)!=funnel_nodes.end() && "ERROR: funnel_reduction::is_funnel: node must be in funnel_nodes");
        auto iter = std::find(funnel_nodes.begin(), funnel_nodes.end(), u);
        std::swap(*iter, funnel_nodes.back());
        funnel_nodes.pop_back();
        if (std::any_of(funnel_nodes.begin(), funnel_nodes.end(), [&](NodeID neighbor) { return weights[neighbor] + weights[u] < weights[v]; })) {
            funnel_nodes.push_back(u);
            continue;
        }

        funnel_set.remove(u);
        if (is_clique(br_alg, funnel_set, funnel_nodes))
        {
            funnel_neighbor = u;
            funnel_nodes.push_back(funnel_neighbor);
            funnel_set.add(funnel_neighbor);
            return true;
        }
        else
        {
            funnel_nodes.push_back(u);
            funnel_set.add(u);
        }
    }
    return false;
}
bool funnel_fold_reduction::is_clique(branch_and_reduce_algorithm* br_alg, fast_set& clique_set, std::vector<NodeID>& clique_nodes) {
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
void funnel_fold_reduction::fold(const fold_data &data, fast_set &funnel_set, branch_and_reduce_algorithm *br_alg)
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
    remaining_neighbors.resize(graph[node].size());

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
        br_alg->set(node, IS_status::included);
        return;
    }
    return;
    status.folded_stack.push_back(get_reduction_type());
    status.reduction_offset += weights[node];

    // add additional edges connecting the remaining neighbors with funnel neighors outside the funnel set
    restore_vec.push_back({data.node, data.funnel_neighbor, remaining_neighbors, {}});
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++)
    {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];
        for (auto neighbor_2nd : outside_funnel_neighbors)
        {
            if (!graph.adjacent(neighbor, neighbor_2nd))
            {
                graph.add_edge_undirected(neighbor, neighbor_2nd);
                restore_vec.back().node_vecs[idx].push_back(neighbor_2nd);
            }
        }
        weights[neighbor] += weights[funnel_neighbor];
        weights[neighbor] -= weights[node];
        br_alg->add_next_level_node(neighbor);
    }

    br_alg->set(funnel_neighbor, IS_status::folded, false);
    br_alg->set(node, IS_status::folded, true);
    br_alg->add_next_level_neighborhood(funnel_neighbor);
}
void funnel_fold_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    // restore nodes
    br_alg->unset(data.node);
    br_alg->unset(data.funnel_neighbor);

    // remove added edges:
    for (size_t idx = 0; idx < data.remaining_neighbors.size(); idx++)
    {
        NodeID neighbor = data.remaining_neighbors[idx];
        for (auto neighbor_2nd : restore_vec.back().node_vecs[idx])
        {
            status.graph.hide_edge_undirected(neighbor, neighbor_2nd);
        }
        status.weights[neighbor] += status.weights[data.node];
        status.weights[neighbor] -= status.weights[data.funnel_neighbor];
    }
    status.reduction_offset -= status.weights[data.node];

    restore_vec.pop_back();
}
void funnel_fold_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &node_status = br_alg->status.node_status;
    auto &is_weight = br_alg->status.is_weight;
    auto &weights = br_alg->status.weights;
    auto data = restore_vec.back();
    auto &remaining = data.remaining_neighbors;

    bool include_node = std::none_of(remaining.begin(), remaining.end(), [&](NodeID neighbor)
                                     { return node_status[neighbor] == IS_status::included; });

    restore(br_alg);

    if (include_node)
    {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
        is_weight += weights[data.node];
    }
    else
    {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
        for (auto neighbor : remaining)
        {
            if (node_status[neighbor] == IS_status::included)
            { // add original weight for neighbors to IS
                is_weight += weights[data.node];
                is_weight -= weights[data.funnel_neighbor];
            }
        }
        is_weight += weights[data.funnel_neighbor];
    }
}

bool twin_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_twin) return false;
    size_t oldn = br_alg->status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            { reduce_vertex(br_alg, v); });

#ifdef REDUCTION_INFO
        reduced_nodes += (oldn - br_alg->status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != br_alg->status.remaining_nodes) std::cout << "twin redu -> " << (oldn - br_alg->status.remaining_nodes) << std::endl;
	return oldn != br_alg->status.remaining_nodes;
}
inline bool twin_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_twin) return false;
    auto &status = br_alg->status;
    auto &neighbors = br_alg->buffers[0];
    auto &twin_candidates_set = br_alg->set_1;
    auto &tmp_set = br_alg->set_2;

    size_t oldn = status.remaining_nodes;
    NodeID twin;
    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    if (try_neighborhood_reduction(v, br_alg, neighbors_weight))
    {
        return oldn != status.remaining_nodes;
    }
    get_neighborhood_vector(v, br_alg, neighbors);

    twin_candidates_set.clear();
    bool candidates_empty = true;

    for (NodeID neighbor : status.graph[neighbors[0]])
    {
        if (is_reduced(neighbor, br_alg))
            continue;
        if (neighbor != v && br_alg->deg(neighbor) == neighbors.size())
        {
            twin_candidates_set.add(neighbor);
            candidates_empty = false;
            twin = neighbor;
        }
    }

    for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++)
    {
        NodeID neighbor = neighbors[i];
        tmp_set.clear();
        candidates_empty = true;

        for (NodeID candidate : status.graph[neighbor])
        {
            if (is_reduced(candidate, br_alg))
                continue;
            if (twin_candidates_set.get(candidate))
            {
                tmp_set.add(candidate);
                candidates_empty = false;
                twin = candidate;
            }
        }

        std::swap(twin_candidates_set, tmp_set);
    }

    if (candidates_empty)
        return false;

    if (status.weights[v] + status.weights[twin] >= neighbors_weight)
    {
        br_alg->set(v, IS_status::included);
        br_alg->set(twin, IS_status::included);
    }
    else
    {
        fold(br_alg, v, twin);
    }

    return oldn != status.remaining_nodes;
}
void twin_reduction::fold(branch_and_reduce_algorithm *br_alg, NodeID main, NodeID twin)
{
    auto &status = br_alg->status;

    restore_vec.push_back({main, twin});

    br_alg->set(twin, IS_status::folded, true);
    status.weights[main] += status.weights[twin];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(main);
    br_alg->add_next_level_neighborhood(main);
}
void twin_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.twin);
    status.weights[data.main] -= status.weights[data.twin];

    restore_vec.pop_back();
}
void twin_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto main = restore_vec.back().main;
    auto twin = restore_vec.back().twin;

    restore(br_alg);

    if (status.node_status[main] == IS_status::included)
    {
        status.node_status[twin] = IS_status::included;
    }
    else
    {
        status.node_status[twin] = IS_status::excluded;
    }
}

