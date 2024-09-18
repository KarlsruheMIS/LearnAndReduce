#include "fold2_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;
bool fold2_reduction::reduce(reduce_algorithm *br_alg)
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
inline bool fold2_reduction::reduce_vertex(reduce_algorithm *br_alg, NodeID v)
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
void fold2_reduction::fold_triangle_mid_weight(reduce_algorithm *br_alg, const fold_nodes &nodes)
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
void fold2_reduction::fold_triangle_min_weight(reduce_algorithm *br_alg, const fold_nodes &nodes)
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
void fold2_reduction::fold_v_shape_mid_weight(reduce_algorithm *br_alg, const fold_nodes &nodes)
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
void fold2_reduction::fold_v_shape_max_weight(reduce_algorithm *br_alg, const fold_nodes &nodes)
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
void fold2_reduction::fold_v_shape_min_weight(reduce_algorithm* br_alg, const fold_nodes& nodes) {
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
void fold2_reduction::restore(reduce_algorithm *br_alg)
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
void fold2_reduction::apply(reduce_algorithm *br_alg)
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
bool fold2_reduction::is_suited(NodeID v, reduce_algorithm *br_alg)
{
    return br_alg->deg(v) == 2;
}

inline int fold2_reduction::generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
     if (br_alg->deg(v) == 2)
     {
        label.push_back(v);
        return true;
     }
     return false;
}