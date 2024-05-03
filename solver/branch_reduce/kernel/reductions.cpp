/******************************************************************************
* reductions.cpp
*
* Copyright (C) 2015-2018 Robert Williger
*
* This program is free software: you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the Free
* Software Foundation, either version 2 of the License, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along with
* this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/


#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "flow_graph.h"
#include "push_relabel.h"
#include "key_functions.h"
#include "solution_check.h"

#include <utility>


typedef branch_and_reduce_algorithm::IS_status IS_status;

// general reduction
template<typename F>
void general_reduction::for_each_changed_vertex(branch_and_reduce_algorithm* br_alg, F f) {
    auto& status = br_alg->status;
    for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
        if (br_alg->config.reduction_time_limit < br_alg->t.elapsed()) return;
        NodeID v = marker.current_vertex(v_idx);

        if (v < status.n && !is_reduced(v, br_alg)) {
                f(v);
        }
    }
}
NodeID general_reduction::get_max_weight_neighbor(NodeID v, branch_and_reduce_algorithm* br_alg) { 
    auto& status = br_alg->status;
    NodeID max_neighbor = status.graph[v][0];
    NodeWeight max_weight = status.weights[max_neighbor];

    for (NodeID neighbor : status.graph[v]) {
        if (status.weights[neighbor] > max_weight && status.node_status[neighbor] == IS_status::not_set) {
            max_weight = status.weights[neighbor];
            max_neighbor = neighbor;
        }
    }

    return max_neighbor;
}
NodeWeight general_reduction::get_neighborhood_weight(NodeID v, branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    return std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&status](NodeWeight sum, NodeID neighbor) { 
        if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor];
        else return sum;
    });
}
void general_reduction::get_neighborhood_set(NodeID v, branch_and_reduce_algorithm* br_alg, fast_set& set) {
    auto& status = br_alg->status;
    set.clear();
    for (NodeID neighbor : status.graph[v] ) {
        if (status.node_status[neighbor] == IS_status::not_set)
            set.add(neighbor);
    }
}
void general_reduction::get_neighborhood_vector(NodeID v, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& vec) {
    auto& status = br_alg->status;
    vec.clear();
    for (NodeID neighbor : status.graph[v]) {
        if (status.node_status[neighbor] == IS_status::not_set)
            vec.push_back(neighbor);
    }
}
bool general_reduction::try_neighborhood_reduction(NodeID v, branch_and_reduce_algorithm* br_alg, NodeWeight neighbors_weight) {
    if (br_alg->config.generate_training_data) return false;
    auto& status = br_alg->status;
    if (status.weights[v] >= neighbors_weight) {
        br_alg->set(v, IS_status::included);
        return true;
    }
    return false;
}
bool general_reduction::solve_induced_subgraph_from_set(NodeWeight& solution, graph_access& graph, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& nodes_vec, const fast_set& nodes_set, sized_vector<NodeID>& reverse_mapping, bool apply_solution) {
    if (nodes_vec.size() == 0) {
        solution = 0;
        return true;
    }
    if (nodes_vec.size() == 1) {
        solution = br_alg->status.weights[nodes_vec[0]];
        return true;
    }
    reverse_mapping.set_size(br_alg->status.n);
    br_alg->build_induced_subgraph(graph, nodes_vec, nodes_set, reverse_mapping);
    return solve_graph(solution, graph, br_alg->config, apply_solution);
}
bool general_reduction::solve_induced_neighborhood_subgraph(NodeWeight& solution, graph_access& neighborhood_graph, branch_and_reduce_algorithm* br_alg, NodeID v, bool apply_solution) {
    br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
    return solve_graph(solution, neighborhood_graph, br_alg->config, apply_solution);
}
bool general_reduction::solve_graph(NodeWeight& solution, graph_access& graph, ReductionConfig &config, bool apply_solution) {
    if (graph.number_of_nodes() == 0 ){
        solution = 0;
        return true;
    } 
    if (graph.number_of_edges() == 0 ) {
        solution = 0;
        forall_nodes(graph, node) {
            if (graph.getNodeWeight(node)>0)
            {
                graph.setPartitionIndex(node, 1);
                solution += graph.getNodeWeight(node);
            }
        } endfor
        return true;
    }
    auto c = config;
    c.time_limit = graph.number_of_nodes() / 10.0;
    // c.time_limit = config.reduction_time_limit*0.1;

    branch_and_reduce_algorithm solver(graph, config, true);
	solver.ch.disable_cout();
    if (!solver.run_branch_reduce()) {
        std::cerr << "%br_call time out" << std::endl;
        return false;
    }
    if (apply_solution) {
        solver.apply_branch_reduce_solution(graph);
    }
	solver.ch.enable_cout();
    solution = solver.get_current_is_weight();
    return true;
}
bool general_reduction::is_reduced(NodeID v, branch_and_reduce_algorithm* br_alg) {
    return br_alg->status.node_status[v] != IS_status::not_set;
}

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_neighborhood) return false;
    size_t oldn = br_alg->status.remaining_nodes;
    br_alg->reduction_timer.restart();

	for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });
    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != br_alg->status.remaining_nodes;
}
bool neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_neighborhood) return false;
    auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    NodeWeight neighbor_weights = this->get_neighborhood_weight(v, br_alg);
    if (status.weights[v] >= neighbor_weights) {
        br_alg->set(v, IS_status::included);
    }
	return oldn != br_alg->status.remaining_nodes;
}
bool neighborhood_reduction::is_suited(NodeID v, branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
	return (status.weights[v] >= std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&](NodeWeight sum, NodeID neighbor) { 
			if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor]; 
			else return sum;}));
}

bool fold1_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_fold1) return false;
	auto& status = br_alg->status;
    br_alg->reduction_timer.restart();
	size_t oldn = status.remaining_nodes;

	for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool fold1_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_fold1) return false;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
    assert(status.node_status[v] == IS_status::not_set);

    if (status.weights[v] == 0) 
    {
        br_alg->set(v, IS_status::excluded);
    } else if (br_alg->deg(v) == 0) 
    {
        br_alg->set(v, IS_status::included);
    } else if (br_alg->deg(v) == 1) 
    {
        NodeID neighbor = status.graph[v][0];
        if (status.weights[neighbor] > status.weights[v]) 
            fold(br_alg, {v, neighbor});
        else 
            br_alg->set(v, IS_status::included, true);
    }

	return oldn != status.remaining_nodes;
}
void fold1_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {

    auto& status = br_alg->status;

	restore_vec.push_back({nodes, status.weights[nodes.deg1_node]});
    br_alg->set(nodes.deg1_node, IS_status::folded);

    status.reduction_offset += status.weights[nodes.deg1_node];
    status.weights[nodes.fold_node] -= status.weights[nodes.deg1_node];

    status.folded_stack.push_back(get_reduction_type());

    br_alg->add_next_level_node(nodes.fold_node);
    br_alg->add_next_level_neighborhood(nodes.fold_node);
}
void fold1_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

    br_alg->unset(data.nodes.deg1_node);

	status.weights[data.nodes.fold_node] += data.deg1_weight;
	status.reduction_offset -= data.deg1_weight; 

	restore_vec.pop_back();
}
void fold1_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto main_status = status.node_status[nodes.fold_node];
	restore(br_alg);

	if (main_status == IS_status::included) {
		status.node_status[nodes.fold_node] = IS_status::included;
		status.node_status[nodes.deg1_node] = IS_status::excluded;

		status.is_weight += status.weights[nodes.fold_node];

	} else if (main_status == IS_status::excluded) {
		status.node_status[nodes.fold_node] = IS_status::excluded;
		status.node_status[nodes.deg1_node] = IS_status::included;

		status.is_weight += status.weights[nodes.deg1_node];
	}
}
bool fold1_reduction::is_suited(NodeID v, branch_and_reduce_algorithm* br_alg) {
    return br_alg->deg(v) <= 1;
}

bool fold2_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_fold2) return false;

	auto& status = br_alg->status;
    br_alg->reduction_timer.restart();
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for_each_changed_vertex(br_alg, [this, &br_alg, &status](NodeID v) {
        reduce_vertex(br_alg, v);
	});
    
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes || oldw != status.reduction_offset;
}
bool fold2_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_fold2) return false;
	if (br_alg->deg(v) != 2) return false;

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

    // set bigger and smaller neighbor
    NodeID bigger  = status.graph[v][0];
    NodeID smaller = status.graph[v][1];
    if (try_neighborhood_reduction(v, br_alg, status.weights[bigger] + status.weights[smaller])) {
        return oldn != status.remaining_nodes;
    }
    if (status.weights[bigger] < status.weights[smaller]) {
        smaller = status.graph[v][0];
        bigger  = status.graph[v][1];
    }

    // check triangle condition
	bool triangle = false;
    for (NodeID neighbor : status.graph[smaller]) {
        if (neighbor == bigger)
        { 
            triangle = true;
            break;
        }
    }
    
    if (triangle && status.weights[bigger] <= status.weights[v]) {
          br_alg->set(v, IS_status::included);
    } else if (triangle && status.weights[bigger] > status.weights[v] && status.weights[smaller] <= status.weights[v]) {
        if (br_alg->config.disable_triangle_mid) return false;
        this->fold_triangle_mid_weight(br_alg, {v, {bigger, smaller}});
    } else if (triangle) {
        if (br_alg->config.disable_triangle_min) return false;
        this->fold_triangle_min_weight(br_alg, {v, {bigger, smaller}});
    } else if (status.weights[v] >= status.weights[bigger]) {
        if (br_alg->config.disable_v_shape_max) return false;
        this->fold_v_shape_max_weight(br_alg, {v, {bigger, smaller}});
    } else if (status.weights[v] >= status.weights[smaller]) 
    {
        if (br_alg->config.disable_v_shape_mid) return false;
        this->fold_v_shape_mid_weight(br_alg, {v, {bigger, smaller}});
    } else {
        assert(status.weights[v] < status.weights[smaller] && "v is not the smallest");
        if (br_alg->config.disable_v_shape_min) return false;
        // if (br_alg->blowing_up) return;
        this->fold_v_shape_min_weight(br_alg, {v, {bigger, smaller}});
    }

    if (v_shape_min_count >= status.modified_stack.capacity() - status.n) { // initial cpapcity is 2n (one n for v_shape min)
        assert(status.modified_stack.capacity() >= status.modified_stack.size() && "stack size error");
        assert(status.folded_stack.capacity() >= status.folded_stack.size() && "stack size error");
        status.modified_stack.resize(status.modified_stack.capacity() + v_shape_min_count);
        status.folded_stack.resize(status.folded_stack.capacity() + v_shape_min_count);
    }
	return oldn != status.remaining_nodes || oldw != status.reduction_offset;
}
void fold2_reduction::fold_triangle_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);
	br_alg->set(nodes.neighbors[1], IS_status::excluded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], fold2_reduction::fold_case::triangle_mid, {} });
    
	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.neighbors[0]] -= status.weights[nodes.deg2_node];

	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.neighbors[0]);
	br_alg->add_next_level_neighborhood(nodes.neighbors[0]);
}
void fold2_reduction::fold_triangle_min_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;

    br_alg->set(nodes.deg2_node, IS_status::folded);

    restore_vec.push_back({nodes, status.weights[nodes.deg2_node], fold2_reduction::fold_case::triangle_min, {} });
    
	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.neighbors[0]] -= status.weights[nodes.deg2_node];
	status.weights[nodes.neighbors[1]] -= status.weights[nodes.deg2_node];

	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(nodes.neighbors[0]);
	br_alg->add_next_level_neighborhood(nodes.neighbors[0]);

}
void fold2_reduction::fold_v_shape_mid_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
    auto oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;
    NodeID bigger  = nodes.neighbors[0];
    NodeID smaller = nodes.neighbors[1];
    NodeWeight deg2_weight = status.weights[nodes.deg2_node];

	restore_vec.push_back({nodes, deg2_weight, fold2_reduction::fold_case::v_shape_mid,{}});
	br_alg->set(nodes.deg2_node, IS_status::folded);

	status.reduction_offset += deg2_weight;
	status.weights[bigger]  -= deg2_weight;
	neighbors.clear();
	neighbors.add(nodes.deg2_node);
	neighbors.add(smaller);
	neighbors.add(bigger);

    std::vector<NodeID> new_neighbors;

    for (auto neighbor : status.graph[smaller]) {
        neighbors.add(neighbor);
    }

    for (auto neighbor : status.graph[bigger]) {
        if (neighbors.add(neighbor)) {
            status.graph.add_edge_undirected(smaller, neighbor);
            restore_vec.back().node_vecs[1].push_back(neighbor);
        }
    }

	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(smaller);
	br_alg->add_next_level_node(bigger);
	br_alg->add_next_level_neighborhood(smaller);

}
void fold2_reduction::fold_v_shape_max_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
    auto oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	restore_vec.push_back({ nodes, status.weights[nodes.deg2_node], fold2_reduction::fold_case::v_shape_max, {}});
	status.reduction_offset += status.weights[nodes.deg2_node];
	status.weights[nodes.deg2_node] = status.weights[nodes.neighbors[0]] + status.weights[nodes.neighbors[1]] - status.weights[nodes.deg2_node];

	neighbors.clear();
	neighbors.add(nodes.deg2_node);

	for (size_t i = 0; i < 2; i++) {
		for (auto neighbor : status.graph[nodes.neighbors[i]]) {
			if (neighbors.add(neighbor)) {
                // TODO: check relink
                // status.graph.relink_directed(neighbor, nodes.neighbors[i], nodes.deg2_node);
                // status.graph.add_edge_directed(nodes.deg2_node, neighbor);
                status.graph.add_edge_undirected(nodes.deg2_node, neighbor);
			    restore_vec.back().node_vecs[i].push_back(neighbor);
			}
		}
	}
	br_alg->set(nodes.neighbors[1], IS_status::folded,true);
	br_alg->set(nodes.neighbors[0], IS_status::folded,false);

	status.folded_stack.push_back(get_reduction_type());
	br_alg->add_next_level_node(nodes.neighbors[0]);
}
void fold2_reduction::fold_v_shape_min_weight(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
    v_shape_min_count++;

    if (v_shape_min_count >= status.modified_stack.capacity() - status.n) { // initial cpapcity is 2n (one n for v_shape min)
        assert(status.modified_stack.capacity() >= status.modified_stack.size() && "stack size error");
        assert(status.folded_stack.capacity() >= status.folded_stack.size() && "stack size error");
        // std::cout << "resize modified stack " << status.modified_stack.capacity(); 
        status.modified_stack.resize(status.modified_stack.capacity() + v_shape_min_count);
        status.folded_stack.resize(status.folded_stack.capacity() + v_shape_min_count);
        // std::cout << " to " << status.modified_stack.capacity();
        // std::cout << " v_shape_min_count " << v_shape_min_count << std::endl;
    }

    status.modified_stack.push_back(br_alg->MODIFIED_TOKEN);

    NodeWeight deg2_weight = status.weights[nodes.deg2_node];
	restore_vec.push_back({ nodes, deg2_weight, fold2_reduction::fold_case::v_shape_min, {}});
	status.reduction_offset += deg2_weight; 
	status.weights[nodes.neighbors[0]] -= deg2_weight;
	status.weights[nodes.neighbors[1]] -= deg2_weight; 

    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[0]);
    status.graph.hide_edge_undirected(nodes.deg2_node, nodes.neighbors[1]);

	neighbors.clear();
    neighbors.add(nodes.deg2_node);

	for (size_t i = 0; i < 2; i++) {
		for (auto neighbor : status.graph[nodes.neighbors[i]]) {
            if (status.node_status[neighbor] != IS_status::not_set) continue;
			if (neighbors.add(neighbor)) {
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
void fold2_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();
    auto& deg2_node = data.nodes.deg2_node;
    auto& nb0 = data.nodes.neighbors[0];
    auto& nb1 = data.nodes.neighbors[1];

	status.reduction_offset -= data.deg2_weight;
    switch (data.fold_case)
    {
        case fold2_reduction::fold_case::triangle_mid:
            br_alg->unset(deg2_node);
	        status.weights[nb0] += data.deg2_weight;
            break;
        case fold2_reduction::fold_case::triangle_min:
            br_alg->unset(deg2_node);
            status.weights[nb0] += data.deg2_weight;
            status.weights[nb1] += data.deg2_weight;
            break;
        case fold2_reduction::fold_case::v_shape_max:
            br_alg->unset(nb1);
            br_alg->unset(nb0);
            for (int i = 0; i < 2; i++) {
                for (NodeID neighbor : data.node_vecs[i]) {
                    // status.graph.relink_directed(neighbor, deg2_node, data.nodes.neighbors[i]);
                    // status.graph.hide_edge(deg2_node, neighbor);
                    status.graph.hide_edge_undirected(deg2_node, neighbor);
                }
            }
            status.weights[deg2_node] = data.deg2_weight;
            break;
        case fold2_reduction::fold_case::v_shape_mid:
            br_alg->unset(deg2_node);
            if (status.node_status[nb0] != IS_status::not_set) br_alg->unset(nb0, false);
            if (status.node_status[nb1] != IS_status::not_set) br_alg->unset(nb1, false);

            // restore neighbors of smaller node
            for (auto neighbor : data.node_vecs[1]) {
                status.graph.hide_edge_undirected(neighbor, nb1);
            }
            status.weights[nb0] += data.deg2_weight;
            break;
        case fold2_reduction::fold_case::v_shape_min:
            if (status.node_status[deg2_node] != IS_status::not_set) br_alg->unset(deg2_node, false);
            if (status.node_status[nb0] != IS_status::not_set) br_alg->unset(nb0, false);
            if (status.node_status[nb1] != IS_status::not_set) br_alg->unset(nb1, false);

	        for (size_t i = 0; i < 2; i++) {
                for (auto second_neighbor : data.node_vecs[i]) {
                    status.graph.hide_edge_undirected(second_neighbor, deg2_node);
                }
                status.graph.add_edge_undirected(data.nodes.neighbors[i], deg2_node);
                status.weights[data.nodes.neighbors[i]] += data.deg2_weight;
	        }  
            v_shape_min_count--;
        }
	restore_vec.pop_back();
}
void fold2_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
    int fold_case = restore_vec.back().fold_case;
    
	auto deg_2_status = status.node_status[nodes.deg2_node];
    auto bigger_status = status.node_status[nodes.neighbors[0]];
    auto smaller_status = status.node_status[nodes.neighbors[1]];
	restore(br_alg);

    switch (fold_case)
    {
        case fold2_reduction::fold_case::v_shape_max:
            if (deg_2_status == IS_status::included) {
                status.node_status[nodes.deg2_node]    = IS_status::excluded;
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.node_status[nodes.neighbors[1]] = IS_status::included;

                status.is_weight += status.weights[nodes.neighbors[0]]; 
                status.is_weight += status.weights[nodes.neighbors[1]];

            } else { // (deg_2_status == IS_status::excluded)
                status.node_status[nodes.deg2_node]    = IS_status::included;
                status.node_status[nodes.neighbors[0]] = IS_status::excluded;
                status.node_status[nodes.neighbors[1]] = IS_status::excluded;

                status.is_weight += status.weights[nodes.deg2_node];
            }
            break;
        case fold2_reduction::fold_case::v_shape_mid:
            status.node_status[nodes.deg2_node]    = IS_status::excluded;
            status.node_status[nodes.neighbors[0]] = IS_status::excluded;
            status.node_status[nodes.neighbors[1]] = IS_status::excluded;

            if (smaller_status == IS_status::included) { // smaller included
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.node_status[nodes.neighbors[1]] = IS_status::included;
                status.is_weight += status.weights[nodes.neighbors[0]];
                status.is_weight += status.weights[nodes.neighbors[1]];
            } else if (bigger_status == IS_status::included) {
                status.node_status[nodes.neighbors[0]] = IS_status::included;
                status.is_weight += status.weights[nodes.neighbors[0]];
            } else if (bigger_status == IS_status::excluded && smaller_status == IS_status::excluded) {
                status.node_status[nodes.deg2_node] = IS_status::included;
                status.is_weight += status.weights[nodes.deg2_node];
            } else if (bigger_status == IS_status::not_set || smaller_status == IS_status::not_set) {
                std::cerr << "ERROR restore mid nodes not set" << std::endl;
            }
            break;
        case fold2_reduction::fold_case::v_shape_min:
	        if (deg_2_status == IS_status::included) {
               status.node_status[nodes.deg2_node]    = IS_status::excluded;
               status.node_status[nodes.neighbors[0]] = IS_status::included;
               status.node_status[nodes.neighbors[1]] = IS_status::included;
               status.is_weight += status.weights[nodes.neighbors[0]];
               status.is_weight += status.weights[nodes.neighbors[1]];
	        } else if (smaller_status == IS_status::included) {
                status.node_status[nodes.deg2_node]    = IS_status::excluded;
           	    status.node_status[nodes.neighbors[1]] = IS_status::included;
          		status.is_weight += status.weights[nodes.neighbors[1]];
                if (bigger_status == IS_status::included) {
           	        status.node_status[nodes.neighbors[0]] = IS_status::included;
               	    status.is_weight += status.weights[nodes.neighbors[0]];
           	    } else if (bigger_status == IS_status::excluded){
           	        status.node_status[nodes.neighbors[0]] = IS_status::excluded;
           	    }
            } else if (smaller_status == IS_status::excluded) {
               status.node_status[nodes.neighbors[1]] = IS_status::excluded;
               if (bigger_status == IS_status::included) {
                  	status.node_status[nodes.deg2_node]    = IS_status::excluded;
           	        status.node_status[nodes.neighbors[0]] = IS_status::included;
                    status.is_weight += status.weights[nodes.neighbors[0]];
                } else if (bigger_status == IS_status::excluded){
               	    status.node_status[nodes.neighbors[0]] = IS_status::excluded;
	                status.node_status[nodes.deg2_node]    = IS_status::included;
	                status.is_weight += status.weights[nodes.deg2_node];
                }
            } else { std::cerr << "Error included v_shape_min_reduction::apply" <<  std::endl; }
            break;

        default: // triangle_mid and triangle_min
            if (bigger_status == IS_status::included) {
            	status.node_status[nodes.neighbors[0]] = IS_status::included;
            	status.node_status[nodes.neighbors[1]] = IS_status::excluded;
            	status.node_status[nodes.deg2_node] = IS_status::excluded;
            	status.is_weight += status.weights[nodes.neighbors[0]];

            } else if (smaller_status == IS_status::included){
            	status.node_status[nodes.neighbors[0]] = IS_status::excluded;
            	status.node_status[nodes.neighbors[1]] = IS_status::included;
            	status.node_status[nodes.deg2_node] = IS_status::excluded;
            	status.is_weight += status.weights[nodes.neighbors[1]];

            } else {
            	status.node_status[nodes.neighbors[0]] = IS_status::excluded;
            	status.node_status[nodes.neighbors[1]] = IS_status::excluded;
            	status.node_status[nodes.deg2_node] = IS_status::included;
            	status.is_weight += status.weights[nodes.deg2_node];
            }
            break;
    }
}
bool fold2_reduction::is_suited(NodeID v, branch_and_reduce_algorithm* br_alg) {
    return br_alg->deg(v) == 2 ;
}

bool single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_basic_se) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_basic_se) return false;
    if (br_alg->deg(v) == 0 ) return false;

	auto& status = br_alg->status;
    auto& graph = br_alg->status.graph;    
    auto& weights = br_alg->status.weights;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, neighbors);

    for (NodeID neighbor : graph[v]) {
        if (is_reduced(neighbor, br_alg)) continue;
        if (weights[v] <= weights[neighbor]) { // otherwise not applicable to this edge
            NodeWeight partial_neighbor_sum = 0;
            for (NodeID second_neighbor : status.graph[neighbor]) {
                if (!is_reduced(second_neighbor, br_alg) && !neighbors.get(second_neighbor)) {
                    partial_neighbor_sum += status.weights[second_neighbor];
                    if (partial_neighbor_sum > weights[neighbor]) 
                        break;
                }
            }
         
            // note: weight of v is in partial_neighbor_sum included
            // if N(neighbor) \subset N(v) partial_neighbor_sum = weights[v]
            if (br_alg->config.generate_training_data && partial_neighbor_sum == weights[v]) break; // should be domination
            if (partial_neighbor_sum <= weights[neighbor]) { 
                br_alg->set(v, IS_status::excluded);
                break;
            }
        }
    }

	return oldn != status.remaining_nodes;
}

bool extended_single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_extended_se) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool extended_single_edge_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_extended_se) return false;
    if (br_alg->deg(v) == 0 ) return false;

    auto& config = br_alg->config;
	auto& status = br_alg->status;
    auto& weights = status.weights;
	size_t oldn = status.remaining_nodes;
	auto& neighbors = br_alg->set_1;
	auto& neighbors_vec = br_alg->buffers[0];


    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    NodeID max_weight_neighbor = get_max_weight_neighbor(v, br_alg);
    assert(max_weight_neighbor != v && "ERROR: max_weight_neighbor == v");
    if  (weights[v] < neighbors_weight - weights[max_weight_neighbor]) {
        return false;
    }
    if (!config.generate_training_data) {
        if (try_neighborhood_reduction(v, br_alg, neighbors_weight)) {
            return oldn != br_alg->status.remaining_nodes;
        }
    }
    get_neighborhood_set(v, br_alg, neighbors);
    get_neighborhood_vector(v, br_alg, neighbors_vec);

    std::sort(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b) { return weights[a] > weights[b]; });
    for (NodeID max_neighbor : neighbors_vec) 
    {
        if (!config.generate_training_data) 
        {
            if (v > max_neighbor) continue;
        }
        if (weights[v] < neighbors_weight - weights[max_neighbor]) 
            break;
        
        bool progress = false;
        for (NodeID neighbor : status.graph[max_neighbor]) 
        {
            if (neighbor == v)  continue; 
            if (is_reduced(neighbor, br_alg)) continue;
            // exclude neighborhood intersection and update neighborhood
            if (neighbors.get(neighbor)) { 
                br_alg->set(neighbor, IS_status::excluded);
                progress = true;
            }
        }

        if (progress) break;
    }

	return oldn != status.remaining_nodes;
}

bool domination_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if(br_alg->config.disable_domination) return false;
	auto& status = br_alg->status;
    br_alg->reduction_timer.restart();
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool domination_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if(br_alg->config.disable_domination) return false;
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
    if (!br_alg->config.generate_training_data) {
        NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
        if (try_neighborhood_reduction(v, br_alg, neighbor_weights)) {
            return oldn != status.remaining_nodes;
        }
    }

    size_t neighbors_count = 0;
    bool is_subset;

    neighbors.clear();
    neighbors.add(v);
    for (NodeID neighbor : status.graph[v]) {
        neighbors.add(neighbor);
        neighbors_count++;
    }

    for (NodeID neighbor : status.graph[v]) {
        if (br_alg->deg(neighbor) > neighbors_count)
            continue;
        if (status.weights[neighbor] < status.weights[v])
            continue;

        is_subset = true;

        for (NodeID neighbor2 : status.graph[neighbor]) {
            if (!neighbors.get(neighbor2)) {
                is_subset = false;
                break;
            }
        }

        // if (is_subset && status.weights[neighbor] >= status.weights[v]) {
        if (is_subset) {
            br_alg->set(v, IS_status::excluded);
            break;
        }
    }
    
	return oldn != status.remaining_nodes;
}
bool cut_vertex_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->blowing_up) return false;
    if (br_alg->config.disable_cut_vertex) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
    auto& visited = br_alg->bool_buffer;
    auto& cut_component = br_alg->buffers[0];
    auto& cut_component_set = br_alg->set_1;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& reverse_mapping = br_alg->buffers[1];
    auto config = br_alg->config;
    graph_access cut_graph;
	size_t oldn = status.remaining_nodes;
    NodeID cut_v;
    std::vector<NodeID> unreduced_mapping(status.remaining_nodes, status.n);
    fast_set tested_cut_vertices(status.n);

    //TODO make dfs only run on reduced nodes (mappings vie buffer)
    visited.assign(status.n, true);
    NodeID local_n = 0;
    for (NodeID n = 0; n < br_alg->status.n; ++n) {
        if (status.node_status[n] == IS_status::not_set) visited[n] = false;
        if (!visited[n]) {
            unreduced_mapping[local_n] = n;
            local_n++;
        }
        assert(local_n <= status.remaining_nodes && "ERROR: cut_vertex_reduction::reduce: local_n too small");
    }

    while (find_cut_vertex(br_alg, cut_v, cut_component, unreduced_mapping, tested_cut_vertices) && br_alg->config.reduction_time_limit > br_alg->t.elapsed())    {
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1) return false; //fold1
        assert(cut_component.size() <= config.subgraph_node_limit);

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component) {
            cut_component_set.add(n);
        }

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut_v = false;
        for (NodeID neighbor : status.graph[cut_v]) {
            if (!cut_component_set.get(neighbor))
                real_cut_v = true;
        }

        if (!real_cut_v) { //directly solve the component without fold
            if (config.generate_training_data)
                return false;
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            graph_access component_graph;
            if (!solve_induced_subgraph_from_set(MWIS_weight, component_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
                return false;

            for (NodeID node : cut_component)
            {
                assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
                if (component_graph.getPartitionIndex(reverse_mapping[node]) == 1) 
                    br_alg->set(node, IS_status::included);
            }
        } else {
            NodeWeight small_cutMWIS_weight = 0;
            NodeWeight large_cutMWIS_weight = 0;
            sized_vector<NodeID> cut_v_included_i(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_included_e(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_excluded_i(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_excluded_e(config.subgraph_node_limit);
            if (get_fold_data(br_alg, cut_v, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e, large_cutMWIS_weight, small_cutMWIS_weight)) 
            {
                fold_data data = {cut_v, status.weights[cut_v], large_cutMWIS_weight, small_cutMWIS_weight, cut_component};
                fold(br_alg, data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e);
            }
        }
    }

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
    br_alg->config.disable_cut_vertex = true; // only test whole graph once
	return oldn != status.remaining_nodes;
}
bool cut_vertex_reduction::find_cut_vertex(branch_and_reduce_algorithm* br_alg, NodeID& cut_v, sized_vector<NodeID>& cut_component, std::vector<NodeID>& reverse_mapping, fast_set& tested) {
    auto& status = br_alg->status;
    auto& visited = br_alg->bool_buffer;

    int step = 0;
    for (NodeID local_n = 0; local_n < br_alg->status.remaining_nodes; ++local_n) {
            if (!visited[reverse_mapping[local_n]]) {
                if (DFS(br_alg, reverse_mapping[local_n], step, cut_v, cut_component) && !tested.get(cut_v)) {
                        cut_v = cut_v;
                        tested.add(cut_v);
                        return true; // New articulation point with small component found
                }
                if (br_alg->config.reduction_time_limit < br_alg->t.elapsed()) return false;
            }
        }
    return false; // No articulation point meeting criteria found
}
bool cut_vertex_reduction::DFS(branch_and_reduce_algorithm* br_alg, NodeID u, int& step, NodeID& cut_vertex, sized_vector<NodeID>& smallComponent) {
    auto& status = br_alg->status;
    auto& graph = br_alg->status.graph;
    auto& visited = br_alg->bool_buffer;
    auto& low = br_alg->buffers[0];
    auto& disc = br_alg->buffers[1];
    auto& parent = br_alg->buffers[2];
    low.assign(status.n,INT_MAX);
    disc.assign(status.n,0);
    parent.assign(status.n, status.n);
    std::vector<NodeID> stack;
    std::vector<std::pair<NodeID,NodeID>> edge_stack;
    stack.reserve(br_alg->status.remaining_nodes);
    edge_stack.reserve(br_alg->status.remaining_nodes);
    stack.push_back(u);

    while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
        NodeID current = stack.back();

        int children = 0;
        disc[current] = low[current] = ++step;
        visited[current] = true;

        bool pushedChild = false;
        for (NodeID v : graph[current]) {
            // Avoid reprocessing an edge in the undirected graph
            if (!edge_stack.empty() && edge_stack.back() == std::make_pair(v, current)) continue;

            if (!visited[v]) {
                stack.push_back(v);
                edge_stack.push_back(std::make_pair(current, v)); // Track the edge
                children++;
                parent[v] = current;
                pushedChild = true;
                break;
            } else if (v != parent[current]) {
                low[current] = std::min(low[current], disc[v]);
            }
        }
        if (!pushedChild) 
        {
           stack.pop_back(); // Finish processing current node 
           if (!stack.empty()) {
                low[parent[current]] = std::min(low[parent[current]], low[current]);
            
                // check if current is a cut vertex
                if ((parent[current] == status.n && disc[current] != low[current]) ||  // root node
                    (parent[current] != status.n && low[current] >= disc[parent[current]])) 
                    if (check_components(br_alg, current, cut_vertex, smallComponent)) {
                        return true;
                    }
           }
            if (!edge_stack.empty()) edge_stack.pop_back(); // Pop the edge after processing
        }
    }
    return false;
}
bool cut_vertex_reduction::check_components(branch_and_reduce_algorithm* br_alg, NodeID u, NodeID& cut_vertex, sized_vector<NodeID>& smallComponent) {
    auto& status = br_alg->status;
    auto& config = br_alg->config;

    // Mark all nodes that are not to be considered as already visited
    std::vector<bool> visited_neighbors(status.n, true);
    for (NodeID n = 0; n < status.n; ++n) {
        if (status.node_status[n] == IS_status::not_set) {
            visited_neighbors[n] = false;
        }
    }

    visited_neighbors[u] = true;
    cut_vertex = u;

    for (NodeID n : status.graph[u]){
        if (!visited_neighbors[n]) {
            if (br_alg->config.reduction_time_limit < br_alg->t.elapsed()) return false;
            smallComponent.clear();

            if (!build_small_component(n, br_alg, smallComponent, visited_neighbors)) {
                // If component is too large, early terminate and set the visited status for this path
                dfs_fill_visited(n, br_alg, visited_neighbors);
            } else if (smallComponent.size() <= config.subgraph_node_limit) {
                // A small enough component was found, no need to check further
                return true;
            }
        }
    }
    // After exploring all components and no valid small component was found
    return false;
}
bool cut_vertex_reduction::build_small_component(NodeID u, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& component, std::vector<bool>& visited) 
{
    auto& status = br_alg->status;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty()) {
        NodeID current = stack.back();
        stack.pop_back();

        if (visited[current]) continue;

        visited[current] = true;
        component.push_back(current);

        // If the component exceeds the maximum size, return false to indicate failure
        if (component.size() > br_alg->config.subgraph_node_limit) {
            for (NodeID n : component) visited[n] = false;
            component.clear();
            return false;
        }

        for (NodeID neighbor : status.graph[current]) {
            if (!visited[neighbor]) {
                stack.push_back(neighbor);
            }
        }
    }
    return true;
}
void cut_vertex_reduction::dfs_fill_visited(NodeID u, branch_and_reduce_algorithm* br_alg, std::vector<bool>& visited) {
    auto& status = br_alg->status;
    auto& visited_dfs = br_alg->bool_buffer;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
        NodeID current = stack.back();
        stack.pop_back();

        if (visited[current]) {
            continue;
        }

        // Mark the current node as visited
        visited[current] = true;
        visited_dfs[current] = true;

        for (NodeID neighbor : br_alg->status.graph[current]) {
            if (!visited[neighbor]) {
                stack.push_back(neighbor);
            }
        }
    }
}
void cut_vertex_reduction::fold(branch_and_reduce_algorithm* br_alg, fold_data& data, sized_vector<NodeID>& cut_v_included_i, sized_vector<NodeID>& cut_v_included_e, sized_vector<NodeID>& cut_v_excluded_i, sized_vector<NodeID>& cut_v_excluded_e)
{
    auto& status = br_alg->status;
	restore_vec.push_back({data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e});
	for (int i = data.cut_component.size() - 1; i >= 1; i--) {
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
void cut_vertex_reduction::restore(branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    auto& data = restore_vec.back().data;

    status.reduction_offset -= data.large_cutMWIS_weight;
    status.weights[data.cut_vertex] = data.cut_vertex_weight;

    for (int i = 0; i < data.cut_component.size(); i++) {
        br_alg->unset(data.cut_component[i]);
    }

    restore_vec.pop_back();
}
void cut_vertex_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& node_status = br_alg->status.node_status;
    NodeWeight& is_weight = br_alg->status.is_weight;
    auto restore_info = restore_vec.back();
    NodeID cut_v = restore_info.data.cut_vertex;
    assert(restore_info.data.large_cutMWIS_weight < restore_info.data.cut_vertex_weight + restore_info.data.small_cutMWIS_weight && "ERROR: cut_vertex_reduction::apply: large_cutMWIS_weight must be larger than small_cutMWIS_weight");

    IS_status cut_status = node_status[cut_v];

	restore(br_alg);

    if (cut_status == IS_status::included) {
        for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_exclude) {
            node_status[neighbor] = IS_status::excluded;
        }
        for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_include) {
            node_status[neighbor] = IS_status::included;
        }
        is_weight += restore_info.data.small_cutMWIS_weight;
        is_weight += restore_info.data.cut_vertex_weight;
    } else {
        for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_exclude) {
            node_status[neighbor] = IS_status::excluded;
        }
        for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_include) {
            node_status[neighbor] = IS_status::included;
        }
        is_weight += restore_info.data.large_cutMWIS_weight;
    }
}

bool cut_vertex_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->blowing_up) return false;
    if (br_alg->config.disable_cut_vertex) return false;
	auto& status = br_alg->status;
    auto& cut_component = br_alg->buffers[0];
    auto& reverse_mapping = br_alg->buffers[1];
    auto& cut_component_set = br_alg->set_1;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& config = br_alg->config;
	size_t oldn = status.remaining_nodes;
    NodeID cut_v;
    if (check_components(br_alg, v, cut_v, cut_component)) {
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1) return false; //fold1
        assert(cut_component.size() <= config.subgraph_node_limit);

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component) {
            cut_component_set.add(n);
        }

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut_v = false;
        for (NodeID neighbor : status.graph[cut_v]) {
            if (!cut_component_set.get(neighbor))
                real_cut_v = true;
        }

        if (!real_cut_v) { //directly solve the component without fold
            if (config.generate_training_data)
                return false;
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            graph_access component_graph;
            if (!solve_induced_subgraph_from_set(MWIS_weight, component_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
                return false;

            for (NodeID node : cut_component)
            {
                assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
                if (component_graph.getPartitionIndex(reverse_mapping[node]) == 1) 
                    br_alg->set(node, IS_status::included);
            }
        } else {
            NodeWeight small_cutMWIS_weight = 0;
            NodeWeight large_cutMWIS_weight = 0;
            sized_vector<NodeID> cut_v_included_i(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_included_e(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_excluded_i(config.subgraph_node_limit);
            sized_vector<NodeID> cut_v_excluded_e(config.subgraph_node_limit);
            if (get_fold_data(br_alg, cut_v, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e, large_cutMWIS_weight, small_cutMWIS_weight)) 
            {
                fold_data data = {cut_v, status.weights[cut_v], large_cutMWIS_weight, small_cutMWIS_weight, cut_component};
                fold(br_alg, data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e);
            }
        }
    }

	return oldn != status.remaining_nodes;
}
bool cut_vertex_reduction::get_fold_data(branch_and_reduce_algorithm* br_alg, NodeID cut_v, sized_vector<NodeID>& cut_v_included_i, sized_vector<NodeID>& cut_v_included_e, sized_vector<NodeID>& cut_v_excluded_i, sized_vector<NodeID>& cut_v_excluded_e, NodeWeight& large_cutMWIS_weight, NodeWeight& small_cutMWIS_weight)
{
	auto& status = br_alg->status;
    auto& config = br_alg->config;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& cut_component_set = br_alg->set_1;
    auto& visited = br_alg->bool_buffer;
    auto& reverse_mapping = br_alg->buffers[1];
    auto& cut_component = br_alg->buffers[0];
    graph_access cut_graph;
    reverse_mapping.assign(status.n, status.n);

    // graph including neighborhood of cut_v
    if (!solve_induced_subgraph_from_set(large_cutMWIS_weight, cut_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
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
        if (reverse_mapping[neighbor] == status.n) continue;
        cut_graph.setNodeWeight(reverse_mapping[neighbor], 0);
        cut_graph.setPartitionIndex(reverse_mapping[neighbor], 0);
    } 
    if (!solve_graph(small_cutMWIS_weight, cut_graph, config, true))
        return false;

    if (status.weights[cut_v] + small_cutMWIS_weight <= large_cutMWIS_weight) // cut_vertex is excluded -> directly apply
    {
        for (NodeID node : cut_v_excluded_i)
        {
            br_alg->set(node, IS_status::included);
        }
        return false;

    } else {
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
            if (cut_graph.getPartitionIndex(reverse_mapping[node]) == 1) {
                cut_v_included_i.push_back(node);
            } else {
                cut_v_included_e.push_back(node);
            }
        }
        return true;
    }
}
// bool cut_vertex_reduction::check_components(branch_and_reduce_algorithm* br_alg, NodeID u, NodeID& cut_vertex, sized_vector<NodeID>& smallComponent) {
//     auto& status = br_alg->status;
//     auto& config = br_alg->config;
//     auto& component_visited = br_alg->bool_buffer;

//     // Mark all nodes that are not to be considered as already visited
//     for (NodeID n = 0; n < status.n; ++n) {
//         if (status.node_status[n] == IS_status::not_set) {
//             component_visited[n] = false;
//         } else {
//             component_visited[n] = true;
//         }
//     }

//     component_visited[u] = true;
//     cut_vertex = u;

//     for (NodeID n : status.graph[u]){
//         if (!component_visited[n]) {
//             if (br_alg->config.reduction_time_limit < br_alg->t.elapsed()) return false;
//             smallComponent.clear();

//             if (!build_small_component(n, br_alg, smallComponent, component_visited)) {
//                 // If component is too large, early terminate and set the visited status for this part
//                 dfs_fill_visited(n, br_alg, component_visited);
//             } else if (smallComponent.size() <= config.subgraph_node_limit) {
//                 // A small enough component was found, no need to check further
//                 return true;
//             }
//         }
//     }
//     // After exploring all components and no valid small component was found
//     return false;
// }
// bool cut_vertex_reduction::build_small_component(NodeID u, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& component, std::vector<bool>& component_visited) 
// {
//     auto& status = br_alg->status;
//     std::vector<NodeID> stack;
//     stack.reserve(status.n);
//     stack.push_back(u);

//     while (!stack.empty()) {
//         NodeID current = stack.back();
//         stack.pop_back();

//         if (component_visited[current]) continue;

//         component_visited[current] = true;
//         component.push_back(current);

//         // If the component exceeds the maximum size, return false to indicate failure
//         if (component.size() > br_alg->config.subgraph_node_limit) {
//             for (NodeID n : component) component_visited[n] = false;
//             component.clear();
//             return false;
//         }

//         for (NodeID neighbor : status.graph[current]) {
//             if (!component_visited[neighbor]) {
//                 stack.push_back(neighbor);
//             }
//         }
//     }
//     return true;
// }
// void cut_vertex_reduction::dfs_fill_visited(NodeID u, branch_and_reduce_algorithm* br_alg, std::vector<bool>& component_visited) {
//     auto& status = br_alg->status;
//     std::vector<NodeID> stack;
//     stack.reserve(status.n);
//     stack.push_back(u);

//     while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
//         NodeID current = stack.back();
//         stack.pop_back();

//         if (component_visited[current]) {
//             continue;
//         }

//         component_visited[current] = true;
//         // if (global_visited.size() > 0) global_visited[current] = true; // not used in reduce by vertex

//         for (NodeID neighbor : br_alg->status.graph[current]) {
//             if (!component_visited[neighbor]) {
//                 stack.push_back(neighbor);
//             }
//         }
//     }
// }
// void cut_vertex_reduction::fold(branch_and_reduce_algorithm* br_alg, fold_data& data, sized_vector<NodeID>& cut_v_included_i, sized_vector<NodeID>& cut_v_included_e, sized_vector<NodeID>& cut_v_excluded_i, sized_vector<NodeID>& cut_v_excluded_e)
// {
//     auto& status = br_alg->status;
// 	restore_vec.push_back({data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e});
// 	for (int i = data.cut_component.size() - 1; i >= 1; i--) {
//         br_alg->set(data.cut_component[i], IS_status::folded, false);
//     }
//     br_alg->set(data.cut_component[0], IS_status::folded);

//     status.reduction_offset += data.large_cutMWIS_weight;
//     status.weights[data.cut_vertex] += data.small_cutMWIS_weight;
//     status.weights[data.cut_vertex] -= data.large_cutMWIS_weight;

//     status.folded_stack.push_back(get_reduction_type());

//     br_alg->add_next_level_node(data.cut_vertex);
//     br_alg->add_next_level_neighborhood(data.cut_vertex);
// }
// void cut_vertex_reduction::restore(branch_and_reduce_algorithm* br_alg) {
//     auto& status = br_alg->status;
//     auto& data = restore_vec.back().data;

//     status.reduction_offset -= data.large_cutMWIS_weight;
//     status.weights[data.cut_vertex] = data.cut_vertex_weight;

//     for (int i = 0; i < data.cut_component.size(); i++) {
//         br_alg->unset(data.cut_component[i]);
//     }

//     restore_vec.pop_back();
// }
// void cut_vertex_reduction::apply(branch_and_reduce_algorithm* br_alg) {
// 	auto& node_status = br_alg->status.node_status;
//     NodeWeight& is_weight = br_alg->status.is_weight;
//     auto restore_info = restore_vec.back();
//     NodeID cut_v = restore_info.data.cut_vertex;
//     assert(restore_info.data.large_cutMWIS_weight < restore_info.data.cut_vertex_weight + restore_info.data.small_cutMWIS_weight && "ERROR: cut_vertex_reduction::apply: large_cutMWIS_weight must be larger than small_cutMWIS_weight");

//     IS_status cut_status = node_status[cut_v];

// 	restore(br_alg);

//     if (cut_status == IS_status::included) {
//         for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_exclude) {
//             node_status[neighbor] = IS_status::excluded;
//         }
//         for (NodeID neighbor : restore_info.case_cut_v_included_nodes_to_include) {
//             node_status[neighbor] = IS_status::included;
//         }
//         is_weight += restore_info.data.small_cutMWIS_weight;
//         is_weight += restore_info.data.cut_vertex_weight;
//     } else {
//         for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_exclude) {
//             node_status[neighbor] = IS_status::excluded;
//         }
//         for (NodeID neighbor : restore_info.case_cut_v_excluded_nodes_to_include) {
//             node_status[neighbor] = IS_status::included;
//         }
//         is_weight += restore_info.data.large_cutMWIS_weight;
//     }
// }


bool component_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->blowing_up) return false;
    if (br_alg->config.disable_cut_vertex) return false;
    br_alg->reduction_timer.restart();
    size_t oldn = br_alg->status.remaining_nodes;
    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != br_alg->status.remaining_nodes;
}
bool component_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->blowing_up) return false;
    if (br_alg->config.disable_cut_vertex) return false;
    if (!br_alg->config.generate_training_data) return false; // only used for training data (otherwise included in cut_vertex)
	auto& status = br_alg->status;
    auto& cut_component = br_alg->buffers[0];
    auto& reverse_mapping = br_alg->buffers[1];
    auto& cut_component_set = br_alg->set_1;
    auto& cut_v_neighbor_set = br_alg->double_set; // since set_2 used in later iterative function call
    auto& config = br_alg->config;
	size_t oldn = status.remaining_nodes;
    NodeID cut_v;
    if (check_components(br_alg, v, cut_v, cut_component)) {
        assert(cut_component.size() <= config.subgraph_node_limit && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1) return false; //fold1
        assert(cut_component.size() <= config.subgraph_node_limit);

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto n : cut_component) {
            cut_component_set.add(n);
        }

        // check if  v actual cut vertex (i.e. has neighbors outside the component)
        bool real_cut_v = false;
        for (NodeID neighbor : status.graph[cut_v]) {
            if (!cut_component_set.get(neighbor))
                real_cut_v = true;
        }

        if (!real_cut_v) { //directly solve the component without fold
            if (config.generate_training_data)
                return false;
            cut_component.push_back(cut_v);
            cut_component_set.add(cut_v);
            NodeWeight MWIS_weight = 0;
            graph_access component_graph;
            if (!solve_induced_subgraph_from_set(MWIS_weight, component_graph, br_alg, cut_component, cut_component_set, reverse_mapping, true))
                return false;

            for (NodeID node : cut_component)
            {
                assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
                if (component_graph.getPartitionIndex(reverse_mapping[node]) == 1) 
                    br_alg->set(node, IS_status::included);
            }
        }
    }

	return oldn != status.remaining_nodes;
}
bool component_reduction::check_components(branch_and_reduce_algorithm* br_alg, NodeID u, NodeID& cut_vertex, sized_vector<NodeID>& smallComponent) {
    auto& status = br_alg->status;
    auto& config = br_alg->config;
    auto& component_visited = br_alg->bool_buffer;

    // Mark all nodes that are not to be considered as already visited
    for (NodeID n = 0; n < status.n; ++n) {
        if (status.node_status[n] == IS_status::not_set) {
            component_visited[n] = false;
        } else {
            component_visited[n] = true;
        }
    }

    component_visited[u] = true;
    cut_vertex = u;

    for (NodeID n : status.graph[u]){
        if (!component_visited[n]) {
            if (br_alg->config.reduction_time_limit < br_alg->t.elapsed()) return false;
            smallComponent.clear();

            if (!build_small_component(n, br_alg, smallComponent, component_visited)) {
                // If component is too large, early terminate and set the visited status for this part
                dfs_fill_visited(n, br_alg, component_visited);
            } else if (smallComponent.size() <= config.subgraph_node_limit) {
                // A small enough component was found, no need to check further
                return true;
            }
        }
    }
    // After exploring all components and no valid small component was found
    return false;
}
bool component_reduction::build_small_component(NodeID u, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& component, std::vector<bool>& component_visited) 
{
    auto& status = br_alg->status;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty()) {
        NodeID current = stack.back();
        stack.pop_back();

        if (component_visited[current]) continue;

        component_visited[current] = true;
        component.push_back(current);

        // If the component exceeds the maximum size, return false to indicate failure
        if (component.size() > br_alg->config.subgraph_node_limit) {
            for (NodeID n : component) component_visited[n] = false;
            component.clear();
            return false;
        }

        for (NodeID neighbor : status.graph[current]) {
            if (!component_visited[neighbor]) {
                stack.push_back(neighbor);
            }
        }
    }
    return true;
}
void component_reduction::dfs_fill_visited(NodeID u, branch_and_reduce_algorithm* br_alg, std::vector<bool>& component_visited) {
    auto& status = br_alg->status;
    std::vector<NodeID> stack;
    stack.reserve(status.n);
    stack.push_back(u);

    while (!stack.empty() && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
        NodeID current = stack.back();
        stack.pop_back();

        if (component_visited[current]) {
            continue;
        }

        component_visited[current] = true;

        for (NodeID neighbor : br_alg->status.graph[current]) {
            if (!component_visited[neighbor]) {
                stack.push_back(neighbor);
            }
        }
    }
}


bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_clique_neighborhood_fast) return false;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
	});

	return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction_fast::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_clique_neighborhood_fast) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
    auto& neighbors = br_alg->buffers[0];
    auto& neighborhood = br_alg->set_1;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (this->try_neighborhood_reduction(v, br_alg, neighbor_weights)) {
        return oldn != br_alg->status.remaining_nodes;
    } 
    this->get_neighborhood_vector(v, br_alg, neighbors);
    this->get_neighborhood_set(v, br_alg, neighborhood);

    std::sort(neighbors.begin(), neighbors.end(), [&status](const NodeID lhs, const NodeID rhs) { return status.weights[lhs] > status.weights[rhs]; });

    bool is_reducible = false;

    for (size_t i = 0; i < neighbors.size() && !is_reducible; i++) {
        NodeID neighbor1 = neighbors[i];

        if (!neighborhood.get(neighbor1)) continue;

        for (NodeID neighbor2 : status.graph[neighbor1]) {
            if (neighbor2 != neighbor1 && neighborhood.get(neighbor2)) {
                // triangle [v, neighbor1, neighbor2] found
                neighbor_weights -= std::min(status.weights[neighbor1], status.weights[neighbor2]);
                neighborhood.remove(neighbor1);
                neighborhood.remove(neighbor2);

                if (status.weights[v] >= neighbor_weights) {
                    is_reducible = true;
                }

                break;
            }
        }
    }

    if (is_reducible) {
        br_alg->set(v, IS_status::included);
    }

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_clique_neighborhood) return false;
	auto& status = br_alg->status;
    br_alg->reduction_timer.restart();
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_clique_neighborhood) return false;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

    if (partition_into_cliques(v, br_alg)) 
    {
		br_alg->set(v, IS_status::included);
    }

	return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction::partition_into_cliques(NodeID v, branch_and_reduce_algorithm* br_alg) {
	auto& weights= br_alg->status.weights;
	auto& neighbors_vec = br_alg->buffers[0];
	auto& clique_neighbors_set = br_alg->set_1;
    target_weight = weights[v];

    neighbor_weights = get_neighborhood_weight(v, br_alg);
	if (neighbor_weights <= target_weight && br_alg->config.generate_training_data)  return false; 
	if (neighbor_weights <= target_weight && !br_alg->config.generate_training_data)  return true; 
    get_neighborhood_vector(v, br_alg, neighbors_vec);

	// partition neigbors of v into cliques
	while (neighbors_vec.size() >= 2 && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
        NodeID max_neighbor = *std::max_element(neighbors_vec.begin(), neighbors_vec.end(), [&](NodeID a, NodeID b) { return weights[a] < weights[b]; });
        clique_neighbors_set.clear();
        for (auto neighbor : neighbors_vec) {
            clique_neighbors_set.add(neighbor);
        }
		clique_neighbors_set.remove(max_neighbor);

		if (expand_clique(max_neighbor, neighbors_vec, clique_neighbors_set, br_alg))
			return true;
	}

	return false;
}
bool clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set, branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& temp_set = br_alg->set_2;
	auto& clique_set = br_alg->double_set;
	clique_set.clear();

	size_t local_max = status.n;
	NodeWeight local_max_weight = 0;
	bool intersection_empty = true;

	while (true) {
		// intersect neighbors of clique with neighbors of max_neighbor
		intersection_empty = true;
		local_max_weight = 0;
		clique_set.add(max_neighbor);
		temp_set.clear();

		for (auto neighbor : status.graph[max_neighbor]) {
			if (clique_neighbors_set.get(neighbor)) {
				temp_set.add(neighbor);
				intersection_empty = false;

				if (status.weights[neighbor] > local_max_weight) {
					local_max = neighbor;
					local_max_weight = status.weights[neighbor];
				}
			}
		}
		if (intersection_empty)	break;
        // if (local_max == status.n) break;

		// add local_max to current clique
		neighbor_weights -= local_max_weight;
		if (neighbor_weights <= target_weight)	return true;

		std::swap(clique_neighbors_set, temp_set);
		clique_neighbors_set.remove(local_max);
		max_neighbor = local_max;
	}

	auto& reamining_neighbors = br_alg->buffers[1];
	reamining_neighbors.clear();

	// adjust neigbors_vec
	for (auto neighbor : neighbors_vec) {
		if (!clique_set.get(neighbor))
			reamining_neighbors.push_back(neighbor);
	}

	std::swap(reamining_neighbors, neighbors_vec);
	return false;
}

bool critical_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->blowing_up) return false;
    if (br_alg->config.disable_critical_set) return false;
    br_alg->reduction_timer.restart();

    auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    auto &mapping = br_alg->buffers[0];
    auto &inverse_mapping = br_alg->buffers[1];

    NodeID n = 0;
    for (NodeID node = 0; node < br_alg->status.n; node++) {
        if (status.node_status[node] == IS_status::not_set) {
            mapping[node] = n;
            inverse_mapping[n] = node;
            ++n;
        }
    }

    // build bipartite flow graph
    // node '+ n' shows that we refer to the node in the rest[1] partition
    flow_graph fg;
    fg.start_construction(2 * n + 2);

    const NodeID source = 2 * n;
    const NodeID sink = source + 1;

    for (NodeID id = 0; id < n; id++) {
        NodeID node = inverse_mapping[id];
        // add source and target edges
        fg.new_edge(source, id, status.weights[node]);
        fg.new_edge(id + n, sink, status.weights[node]);

        // add edges between node and its neighbors in the other partition
        for (NodeID neighbor : status.graph[node]) {
            // each outgoing edge has enough capacity to support full flow of the single incoming edge into 'node'
            fg.new_edge(id, mapping[neighbor] + n, status.weights[node]);
        }
    }
    fg.finish_construction();


    // solve max-flow problem
    push_relabel flow_solver;
    std::vector<NodeID> dummy_vec;
    flow_solver.solve_max_flow_min_cut(fg, source, sink, false, dummy_vec);

    auto& max_cs_set = br_alg->double_set;

    max_cs_set.clear();
    // (source, node) edges where flow < capacity indicate that node is in the maximum critical set
    forall_out_edges(fg, edge, source) {
        NodeID id = fg.getEdgeTarget(source, edge);
        if (fg.getEdgeFlow(source, edge) < fg.getEdgeCapacity(source, edge)) {
            max_cs_set.add(inverse_mapping[id]);
        }
    } endfor

    // isolated nodes in the maximum critical set form the maximum independent critical set
    for (NodeID id = 0; id < n; id++) {
        NodeID node = inverse_mapping[id];
        if (max_cs_set.get(node)) {
            bool isolated = true;

            for (NodeID neighbor : status.graph[node]) {
                if (max_cs_set.get(neighbor)) {
                    isolated = false;
                    break;
                }
            }

            if (isolated) {
                // found isolated node
                br_alg->set(node, IS_status::included);
            }
        }
    }

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
    return oldn != status.remaining_nodes;
}

bool clique_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_clique) return false;
	auto& status = br_alg->status;
    br_alg->reduction_timer.restart();
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID node) {
        reduce_vertex(br_alg, node);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool clique_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_clique) return false;
	auto& status = br_alg->status;
	auto& set_1 = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	auto& isolated = br_alg->buffers[1];
	std::vector<NodeID> non_isolated;
	size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, set_1);
    get_neighborhood_vector(v, br_alg, neighbors);
    set_1.add(v);

    // check if clique
    non_isolated.clear();
    isolated.clear();
    isolated.push_back(v);

    size_t max_isolated_idx = 0;
    weighted_node max_isolated{ v, status.weights[v] };
    weighted_node max_non_isolated{ 0, 0 };

    for (auto neighbor : neighbors) {
        size_t count = 0;
        bool is_isolated = true;

        for (NodeID neighbor_2nd : status.graph[neighbor]) {
            if (set_1.get(neighbor_2nd)) count++;
            else is_isolated = false;
        }

        if (is_isolated) {
            isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_isolated.weight) {
                max_isolated = { neighbor, status.weights[neighbor] };
                max_isolated_idx = isolated.size() - 1;
            }
        }
        else {
            non_isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_non_isolated.weight) {
                max_non_isolated = { neighbor, status.weights[neighbor] };
            }
        }

        if (count != neighbors.size()) return false;
    }

    // one of "isolated" members has highest weight of clique: Add to IS
    // also handles completely isolated cliques
    if (max_isolated.weight >= max_non_isolated.weight) {
        br_alg->set(max_isolated.node, IS_status::included);
        return oldn != br_alg->status.remaining_nodes;
    } 


    // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
    isolated[max_isolated_idx] = isolated.back();
    isolated.pop_back();

    for (auto neighbor : isolated) {
        br_alg->set(neighbor, IS_status::excluded);
    }

    for (size_t i = 0; i < non_isolated.size(); i++) {
        NodeID neighbor = non_isolated[i];
        if (status.weights[neighbor] <= max_isolated.weight) {
            br_alg->set(neighbor, IS_status::excluded);
            non_isolated[i] = non_isolated.back();
            non_isolated.pop_back();
            i--;
        }
    }

    fold(br_alg, std::move(max_isolated), std::move(non_isolated));

	return oldn != status.remaining_nodes;
}
void clique_reduction::fold(branch_and_reduce_algorithm* br_alg, const weighted_node& isolated, std::vector<NodeID>&& non_isolated) {
	auto& status = br_alg->status;

	br_alg->set(isolated.node, IS_status::folded);
	status.reduction_offset += isolated.weight;

	for (auto node : non_isolated) {
		status.weights[node] -= isolated.weight;
		br_alg->add_next_level_neighborhood(node);
	}

	status.folded_stack.push_back(get_reduction_type());
	br_alg->add_next_level_neighborhood(non_isolated);

	restore_vec.emplace_back(isolated, std::move(non_isolated));
}
void clique_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.isolated.node);
	status.reduction_offset -= data.isolated.weight;

	for (auto node : data.non_isolated) {
		status.weights[node] += data.isolated.weight;
	}

	restore_vec.pop_back();
}
void clique_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto isolated = restore_vec.back().isolated.node;

	bool set_isolated = true;

	for (auto node : restore_vec.back().non_isolated) {
		if (status.node_status[node] == IS_status::included) {
			set_isolated = false;
			break;
		}
	}

	status.is_weight += restore_vec.back().isolated.weight;

	restore(br_alg);

	if (set_isolated) {
		status.node_status[isolated] = IS_status::included;
	}
	else {
		status.node_status[isolated] = IS_status::excluded;
	}
}

bool funnel_reduction::reduce(branch_and_reduce_algorithm* br_alg)
{
    if (br_alg->config.disable_funnel) return false;
	size_t oldn = br_alg->status.remaining_nodes;
    br_alg->reduction_timer.restart();

    for_each_changed_vertex(br_alg, [&](NodeID node) {
        reduce_vertex(br_alg, node);
    });

    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != br_alg->status.remaining_nodes;
}
bool funnel_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v)
{
    if (br_alg->config.disable_funnel) return false;
    if (br_alg->deg(v) <= 2) return false; // fold1 or 2
	auto& weights = br_alg->status.weights;
	auto& graph = br_alg->status.graph;
	auto& remaining_n = br_alg->status.remaining_nodes;
	auto& funnel_set = br_alg->set_1;
    auto& f_neighbors = br_alg->set_2; 
	auto& neighbors = br_alg->buffers[0];
	size_t oldn = remaining_n;
    NodeID funnel_neighbor = br_alg->status.n;

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    // need one vertex of weight >= v
    if (std::any_of(neighbors.begin(), neighbors.end(), [&](NodeID neighbor) { return weights[neighbor] >= weights[v]; })) {
        get_neighborhood_set(v, br_alg, funnel_set);
        funnel_set.add(v);

        if (is_funnel(v, funnel_neighbor, br_alg, funnel_set, neighbors)) {
            fold({v, funnel_neighbor}, funnel_set, br_alg);
        }
    }

	return oldn != remaining_n;
}
bool funnel_reduction::is_funnel(NodeID node, NodeID& funnel_neighbor, branch_and_reduce_algorithm* br_alg, fast_set& funnel_set, sized_vector<NodeID>& funnel_nodes) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    auto& neighbors = br_alg->buffers[1];
    if (is_clique(br_alg, funnel_set, funnel_nodes)) 
        return false;
    funnel_neighbor = status.n;
    for (NodeID v : funnel_nodes) {
        assert(v != node && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(v, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");
        if (weights[v] >= weights[node]) {

            bool skip = false;
            for (auto neighbor :graph[node])
            {
                if (is_reduced(neighbor, br_alg)) continue;
                if (neighbor == v) continue;
                if (weights[neighbor] > weights[node]) { 
                    skip = true;
                    break;
                }
            }
            if(skip) continue;

            funnel_nodes.remove(v);
            funnel_set.remove(v);
            if (is_clique(br_alg, funnel_set, funnel_nodes)) {
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
bool funnel_reduction::is_clique(branch_and_reduce_algorithm* br_alg, fast_set& clique_set, sized_vector<NodeID>& clique_nodes) {
    auto& status = br_alg->status;
    auto& graph = status.graph;

    for (auto neighbor : clique_nodes) {
        size_t count = 0;

        for (NodeID neighbor_2nd : graph[neighbor]) {
            if (status.node_status[neighbor_2nd] != IS_status::not_set) continue;
            if (clique_set.get(neighbor_2nd)) count++;
        }

        if (count != clique_nodes.size()) return false; // not all neighbors are in the clique
    }
    return true;
}
void funnel_reduction::fold(const fold_data& data, fast_set& funnel_set, branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    auto& f_neighbors = br_alg->set_2; 
    NodeID node = data.node;
    NodeID funnel_neighbor = data.funnel_neighbor;

    f_neighbors.clear();
    sized_vector<NodeID> outside_funnel_neighbors(status.graph[funnel_neighbor].size());
    sized_vector<NodeID> remaining_neighbors(graph[node].size());

    //collect neighbors of funnel neighbor and those outside the funnel set
    for (auto neighbor : graph[funnel_neighbor]) {
        if (!is_reduced(neighbor, br_alg)) 
        {
            f_neighbors.add(neighbor);
            if (!funnel_set.get(neighbor)) 
                outside_funnel_neighbors.push_back(neighbor);
        }
    }

    // common neighbors are not part of the IS in either case 
    for (auto neighbor : graph[node]) {
        if (is_reduced(neighbor, br_alg)) continue;
        if (neighbor == funnel_neighbor) continue;
        if (f_neighbors.get(neighbor)) 
        {
            br_alg->set(neighbor, IS_status::excluded);
        }
        else
        {
            remaining_neighbors.push_back(neighbor); 
        }
    }

    if (remaining_neighbors.size() == 0) {
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
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++) {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];
        assert(!is_reduced(neighbor, br_alg) && "ERROR: funnel_reduction::fold: remaining neighbors must be unset");
        assert( node != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        assert( funnel_neighbor != neighbor && "ERROR: funnel_reduction::restore: node must not be in remaining_neighbors");
        for (auto neighbor_2nd : outside_funnel_neighbors) {
            if (!graph.adjacent(neighbor, neighbor_2nd)) {
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
void funnel_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
    auto& data = restore_vec.back();

    br_alg->unset(data.node);
    status.weights[data.funnel_neighbor] += status.weights[data.node];

    // remove added edges:
    for (size_t idx = 0; idx < data.remaining_neighbors.size(); idx++) {
        NodeID neighbor = data.remaining_neighbors[idx];
        for (auto neighbor_2nd : restore_vec.back().node_vecs[idx]) {
            status.graph.hide_edge_undirected(neighbor, neighbor_2nd);
        }
    }
    status.reduction_offset -= status.weights[data.node];

	restore_vec.pop_back();
}
void funnel_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& node_status = br_alg->status.node_status;
	auto& is_weight = br_alg->status.is_weight;
	auto& weights = br_alg->status.weights;
	auto data = restore_vec.back();
	auto& remaining = data.remaining_neighbors;

    bool include_node = node_status[data.funnel_neighbor] == IS_status::excluded;
    if (include_node) 
        include_node = std::none_of(remaining.begin(), remaining.end(), [&](NodeID neighbor) { return node_status[neighbor] == IS_status::included; });

    restore(br_alg);

    if (include_node) {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
        is_weight += weights[data.node];
    } else {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
        is_weight += weights[data.node];
    }
    // std::cout << "node: " << data.node << " (" << node_status[data.node] << ") funnel_neighbor: " << data.funnel_neighbor << " (" << node_status[data.funnel_neighbor] << ")\n";
}

bool funnel_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg)
{
    if (br_alg->config.disable_funnel) return false;
	size_t oldn = br_alg->status.remaining_nodes;
    br_alg->reduction_timer.restart();

    for_each_changed_vertex(br_alg, [&](NodeID node) {
        reduce_vertex(br_alg, node);
    });

    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != br_alg->status.remaining_nodes;
}
bool funnel_fold_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v)
{
    if (br_alg->config.disable_funnel_fold) return false;
    if (br_alg->deg(v) <= 2) return false; // fold1 or 2
	auto& weights = br_alg->status.weights;
	auto& remaining_n = br_alg->status.remaining_nodes;
	auto& funnel_set = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	size_t oldn = remaining_n;
    NodeID funnel_neighbor = br_alg->status.n;
    assert(!is_reduced(v, br_alg) && "ERROR: funnel_reduction::reduce_vertex: node must be unset");

    get_neighborhood_vector(v, br_alg, neighbors);
    funnel_set.clear();
    bool skip = false;

    if (std::none_of(neighbors.begin(), neighbors.end(), [&](NodeID neighbor) { return weights[neighbor] > weights[v]; })) {

        get_neighborhood_set(v, br_alg, funnel_set);
        funnel_set.add(v);


        if (is_funnel(v, funnel_neighbor, br_alg, funnel_set, neighbors)) {
            fold({v, funnel_neighbor}, funnel_set, br_alg);
        }
    }

	return oldn != remaining_n;
}
bool funnel_fold_reduction::is_funnel(NodeID v, NodeID& funnel_neighbor, branch_and_reduce_algorithm* br_alg, fast_set& funnel_set, sized_vector<NodeID>& funnel_nodes) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    auto& neighbors = br_alg->buffers[1]; // neighbors of v (set before)
    funnel_neighbor = status.n;
    if (is_clique(br_alg, funnel_set, funnel_nodes)) 
        return false;
    
    for (NodeID u : funnel_nodes) {
        assert(u != v && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(u, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");
        bool skip = false;
        funnel_nodes.remove(std::find(funnel_nodes.begin(), funnel_nodes.end(), u));
        if (std::any_of(funnel_nodes.begin(), funnel_nodes.end(), [&](NodeID neighbor) { return weights[neighbor] + weights[u] < weights[v]; })) {
            funnel_nodes.push_back(u);
            continue;
        }
        #ifdef DEBUG

        for (auto neighbor : funnel_nodes) {
            assert(weights[neighbor] + weights[u] >= weights[v] && "ERROR: funnel_reduction::fold: weight must be larger than node weight");
        }
        #endif

        funnel_set.remove(u);
        if (is_clique(br_alg, funnel_set, funnel_nodes)) {
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
bool funnel_fold_reduction::is_clique(branch_and_reduce_algorithm* br_alg, fast_set& clique_set, sized_vector<NodeID>& clique_nodes) {
    auto& status = br_alg->status;
    auto& graph = status.graph;

    for (auto neighbor : clique_nodes) {
        size_t count = 0;

        for (NodeID neighbor_2nd : graph[neighbor]) {
            if (status.node_status[neighbor_2nd] != IS_status::not_set) continue;
            if (clique_set.get(neighbor_2nd)) count++;
        }

        if (count != clique_nodes.size()) return false; // not all neighbors are in the clique
    }
    return true;
}
void funnel_fold_reduction::fold(const fold_data& data, fast_set& funnel_set, branch_and_reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    auto& f_neighbors = br_alg->set_2; 
    NodeID node = data.node;
    NodeID funnel_neighbor = data.funnel_neighbor;

    f_neighbors.clear();
    sized_vector<NodeID> outside_funnel_neighbors(status.graph[funnel_neighbor].size());
    sized_vector<NodeID> remaining_neighbors(graph[node].size());

    //collect neighbors of funnel neighbor and those outside the funnel set
    for (auto neighbor : graph[funnel_neighbor]) {
        if (!is_reduced(neighbor, br_alg)) 
        {
            f_neighbors.add(neighbor);
            if (!funnel_set.get(neighbor)) 
                outside_funnel_neighbors.push_back(neighbor);
        }
    }

    // common neighbors are not part of the IS in either case 
    for (auto neighbor : graph[node]) {
        if (is_reduced(neighbor, br_alg)) continue;
        if (neighbor == funnel_neighbor) continue;
        if (f_neighbors.get(neighbor)) 
        {
            br_alg->set(neighbor, IS_status::excluded);
        }
        else
        {
            remaining_neighbors.push_back(neighbor); 
        }
    }

    if (remaining_neighbors.size() == 0) {
        assert(br_alg->deg(node) == 1 && "ERROR: funnel_reduction::fold: remaining neighbors must be non-empty");
        br_alg->set(node, IS_status::included);
        return;
    }
    return;
    status.folded_stack.push_back(get_reduction_type());
    status.reduction_offset += weights[node];

    // add additional edges connecting the remaining neighbors with funnel neighors outside the funnel set  
    restore_vec.push_back({data.node, data.funnel_neighbor, remaining_neighbors, {}});
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++) {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];
        for (auto neighbor_2nd : outside_funnel_neighbors) {
            if (!graph.adjacent(neighbor, neighbor_2nd)) {
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
void funnel_fold_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
    auto& data = restore_vec.back();

    // restore nodes 
    br_alg->unset(data.node);
    br_alg->unset(data.funnel_neighbor);

    // remove added edges:
    for (size_t idx = 0; idx < data.remaining_neighbors.size(); idx++) {
        NodeID neighbor = data.remaining_neighbors[idx];
        for (auto neighbor_2nd : restore_vec.back().node_vecs[idx]) {
            status.graph.hide_edge_undirected(neighbor, neighbor_2nd);
        }
        status.weights[neighbor] += status.weights[data.node];
        status.weights[neighbor] -= status.weights[data.funnel_neighbor];
    }
    status.reduction_offset -= status.weights[data.node];

	restore_vec.pop_back();
}
void funnel_fold_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& node_status = br_alg->status.node_status;
	auto& is_weight = br_alg->status.is_weight;
	auto& weights = br_alg->status.weights;
	auto data = restore_vec.back();
	auto& remaining = data.remaining_neighbors;

    bool include_node = std::none_of(remaining.begin(), remaining.end(), [&](NodeID neighbor) { return node_status[neighbor] == IS_status::included; });

    restore(br_alg);

    if (include_node) {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
        is_weight += weights[data.node];
    } else {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
        for (auto neighbor : remaining) {
            if (node_status[neighbor] == IS_status::included)
                is_weight += weights[neighbor];
        }
        is_weight += weights[data.funnel_neighbor];
    }
    // std::cout << "node: " << data.node << " (" << node_status[data.node] << ") funnel_neighbor: " << data.funnel_neighbor << " (" << node_status[data.funnel_neighbor] << ")\n";
}


bool twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_twin) return false;
	size_t oldn = br_alg->status.remaining_nodes;
    br_alg->reduction_timer.restart();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != br_alg->status.remaining_nodes;
}
bool twin_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_twin) return false;
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;

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

    for (NodeID neighbor : status.graph[neighbors[0]]) {
        if (is_reduced(neighbor, br_alg)) continue;
        if (neighbor != v && br_alg->deg(neighbor) == neighbors.size()) {
            twin_candidates_set.add(neighbor);
            candidates_empty = false;
            twin = neighbor;
        }
    }

    for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++) {
        NodeID neighbor = neighbors[i];
        tmp_set.clear();
        candidates_empty = true;

        for (NodeID candidate : status.graph[neighbor]) {
            if (is_reduced(candidate, br_alg)) continue;
            if (twin_candidates_set.get(candidate)) {
                tmp_set.add(candidate);
                candidates_empty = false;
                twin = candidate;
            }
        }

        std::swap(twin_candidates_set, tmp_set);
    }

    if (candidates_empty) return false;

    if (status.weights[v] + status.weights[twin] >= neighbors_weight) {
        br_alg->set(v, IS_status::included);
        br_alg->set(twin, IS_status::included);
    } else {
        fold(br_alg, v, twin);
    }

	return oldn != status.remaining_nodes;
}
void twin_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin) {
	auto& status = br_alg->status;

	restore_vec.push_back({ main, twin });

	br_alg->set(twin, IS_status::folded, true);
	status.weights[main] += status.weights[twin];

	status.folded_stack.push_back(get_reduction_type());

	br_alg->add_next_level_node(main);
	br_alg->add_next_level_neighborhood(main);
}
void twin_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.twin);
	status.weights[data.main] -= status.weights[data.twin];

	restore_vec.pop_back();
}
void twin_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto main = restore_vec.back().main;
	auto twin = restore_vec.back().twin;

	restore(br_alg);

	if (status.node_status[main] == IS_status::included) {
		status.node_status[twin] = IS_status::included;
	} else {
		status.node_status[twin] = IS_status::excluded;
	}
}

bool heavy_vertex_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_heavy_vertex) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool heavy_vertex_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    auto& config = br_alg->config;
    if (config.disable_heavy_vertex) return false;
    if (br_alg->deg(v) < 3) return false; // reduce with fold1 or 2

	auto& status = br_alg->status;
	auto& reverse_mapping = br_alg->buffers[1];
	auto& build_graph_neighbors = br_alg->buffers[0];
	auto& build_graph_neighbors_set = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

    if (status.weights[v] < status.weights[get_max_weight_neighbor(v, br_alg)]) 
        return false;

	graph_access subgraph;

    NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
    if (this->try_neighborhood_reduction(v, br_alg, neighbors_weight)) {
        return true; 
    }

    if (br_alg->deg(v) > config.subgraph_node_limit) return false; //subgraph too large

    get_neighborhood_set(v, br_alg, build_graph_neighbors_set);
    get_neighborhood_vector(v, br_alg, build_graph_neighbors);

    // compute MWIS in N(v) 
    ::NodeWeight MWIS_weight = 0;
    bool solved_exact = solve_induced_subgraph_from_set(MWIS_weight, subgraph, br_alg, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);
    if (!solved_exact) return false; 

    if (status.weights[v] >= MWIS_weight) {
        br_alg->set(v, IS_status::included);
    }
    return oldn != br_alg->status.remaining_nodes;
}

bool heavy_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_heavy_set) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool heavy_set_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID common_neighbor) {
   auto& config = br_alg->config;
    if (config.disable_heavy_set) return false;
    if (br_alg->blowing_up) return false;
    if (br_alg->deg(common_neighbor) < 3) return false; // use other reduction

	auto& status = br_alg->status;
    auto& weights = status.weights;
    assert(status.node_status[common_neighbor] == IS_status::not_set && "ERROR: heavy_set_reduction::reduce_vertex: node must be unset");
    // check if common neighbor is suitable
    if (std::all_of(status.graph[common_neighbor].begin(), status.graph[common_neighbor].end(), [&](NodeID neighbor) { return br_alg->deg(neighbor) > config.subgraph_node_limit; })) 
        return false;
	size_t oldn = status.remaining_nodes;
	auto& v_neighbors_set = br_alg->set_1;
    auto& candidates = status.graph[common_neighbor];
    for (size_t v_idx = 0; v_idx < candidates.size()-1; v_idx++)
    {
        NodeID v = candidates[v_idx];
        if (is_reduced(v, br_alg)) continue;
        if (br_alg->deg(v) > config.subgraph_node_limit) continue; // too many neighbors
        if (br_alg->deg(v) < 3) continue; // use other reduction 
        NodeID deg_v = br_alg->deg(v);
        get_neighborhood_set(v, br_alg, v_neighbors_set);

	    //find second heavy vertex u (not adjacent to v)
	    for (size_t u_idx = v_idx + 1; u_idx < candidates.size(); u_idx++)
        {
            NodeID u = candidates[u_idx];
            assert(u != v);
            if (v_neighbors_set.get(u)) continue; // look for non adjacent nodes
	        if (is_reduced(u, br_alg)) continue;
            if (br_alg->deg(u) + deg_v > config.subgraph_node_limit) continue; //subgraph too large
            if (br_alg->deg(u) < 3) continue; // use other reduction 

            if (is_heavy_set(v, v_neighbors_set, u, br_alg)) {
                    // reduction was applied
                    return oldn != status.remaining_nodes;
            }
        }
    }
    
    return oldn != status.remaining_nodes; 
} 
bool heavy_set_reduction::is_heavy_set(NodeID v, fast_set& v_neighbors_set, NodeID u, branch_and_reduce_algorithm* br_alg) {
    auto& config = br_alg->config;
    auto& status = br_alg->status;
    auto& weights = status.weights;
	auto& v_neighbors_vec = br_alg->buffers[0];
	auto& u_neighbors_vec = br_alg->buffers[1];
    auto& graph_nodes = br_alg->buffers[2];
    auto& graph_nodes_set = br_alg->double_set;
    auto& reverse_mapping = br_alg->buffers[3];
    size_t oldn = status.remaining_nodes;
    reverse_mapping.set_size(status.n);
    graph_nodes.clear();
    graph_nodes_set.clear();
    graph_access subgraph;

    get_neighborhood_vector(v, br_alg,v_neighbors_vec);
    get_neighborhood_vector(u, br_alg,u_neighbors_vec);


    assert(!status.graph.adjacent(v,u) && "ERROR: heavy_set_reduction::is_heavy_set: v and u must be not adjacent");
    assert(!is_reduced(v, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: v must not be reduced");
    assert(!is_reduced(u, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: u must not be reduced");

    for (auto n : v_neighbors_vec) {
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }
    for (auto n : u_neighbors_vec) {
        if (graph_nodes_set.get(n)) continue; // only add nodes once
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }

    std::vector<NodeWeight> MWIS_weights(4,0);
    assert(graph_nodes.size() > 2 && "ERROR: heavy_set_reduction::is_heavy_set: graph_nodes must have at least 3 nodes");
    // compute MWIS in N(v) + N(u) + N(w):
    if (!solve_induced_subgraph_from_set(MWIS_weights[v_combination::oo], subgraph, br_alg, graph_nodes, graph_nodes_set, reverse_mapping)) return false;
    MWIS_weights[v_combination::uv] = weights[u] + weights[v] ;

    if (std::min(weights[v], weights[u]) >= MWIS_weights[v_combination::oo]) {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        return true;
    } else if (MWIS_weights[v_combination::uv] < MWIS_weights[v_combination::oo]) {
        return false;
    }
    bool original_heavy_set = false;
    // compute MWIS[uo] in N(v)\N(u): (included u)
    unset_weights(subgraph, u_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uo], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::uo] > weights[v] ) return false;
    MWIS_weights[v_combination::uo] += weights[u];

    // compute MWIS[5] in N(u)\N(v): (included v)
    set_weights(subgraph, u_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, v_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::ov], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::ov] > weights[u]) return false;
    MWIS_weights[v_combination::ov] += weights[v];

    if (original_heavy_set) {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
    }
    else {
        // get max weights by weight and if equal by participating vertices in the permutation
        v_combination max_weight_combination = static_cast<v_combination>(std::max_element(MWIS_weights.begin(), MWIS_weights.end()) - MWIS_weights.begin());
        NodeWeight best_weight_excluding = 0;
        std::vector<NodeWeight> weights_including(2,0);
        switch (max_weight_combination) 
        {
            case v_combination::uv:
                br_alg->set(v, IS_status::included);
                br_alg->set(u, IS_status::included);
                break;
            case v_combination::uo:
                best_weight_excluding = std::max({MWIS_weights[v_combination::ov], MWIS_weights[v_combination::oo]});
                weights_including = {MWIS_weights[v_combination::uv], MWIS_weights[v_combination::uo]};
                if (std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight) { return weight >= best_weight_excluding; }))
                    br_alg->set(u, IS_status::included);
                break;
            case v_combination::ov:
                best_weight_excluding = std::max({MWIS_weights[v_combination::uo], MWIS_weights[v_combination::oo]});
                weights_including = {MWIS_weights[v_combination::uv], MWIS_weights[v_combination::ov]};
                if (std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight) { return weight >= best_weight_excluding; }))
                    br_alg->set(v, IS_status::included);
                break;
            default:
                break;
        }
    }

    return oldn != status.remaining_nodes;
}
void heavy_set_reduction::unset_weights(graph_access& graph, sized_vector<NodeID>& nodes, sized_vector<NodeID>& reverse_mapping) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], 0);
    }
}
void heavy_set_reduction::set_weights(graph_access& graph, sized_vector<NodeID>& nodes, sized_vector<NodeID>& reverse_mapping, std::vector<NodeWeight>& weights) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], weights[n]);
    }
}

bool heavy_set3_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_heavy_set) return false;

	auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;
    br_alg->reduction_timer.restart();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool heavy_set3_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID common_neighbor) {
    auto& config = br_alg->config;
    if (config.disable_heavy_set) return false;
    if (br_alg->blowing_up) return false;
    if (br_alg->deg(common_neighbor) < 3) return false; // no heavy_set of 3 vertives possible 

	auto& status = br_alg->status;
    auto& weights = status.weights;
    assert(status.node_status[common_neighbor] == IS_status::not_set && "ERROR: heavy_set3_reduction::reduce_vertex: node must be unset");
    // check if common neighbor is suitable
    if (std::all_of(status.graph[common_neighbor].begin(), status.graph[common_neighbor].end(), [&](NodeID neighbor) { return br_alg->deg(neighbor) > config.subgraph_node_limit; })) 
        return false;
	size_t oldn = status.remaining_nodes;
	auto& v_neighbors_set = br_alg->set_1;
	auto& u_neighbors_set = br_alg->set_2;
    auto& candidates = status.graph[common_neighbor];
    // print all neighbors of common neighbor
    for (size_t v_idx = 0; v_idx < candidates.size()-2; v_idx++)
    {
        NodeID v = candidates[v_idx];
        if (is_reduced(v, br_alg)) continue;
        if (br_alg->deg(v) > config.subgraph_node_limit) continue; // too many neighbors
        if (br_alg->deg(v) < 3) continue; // use other reduction 
        NodeID deg_v = br_alg->deg(v);
        get_neighborhood_set(v, br_alg, v_neighbors_set);

	    //find second heavy vertex u (not adjacent to v)
	    for (size_t u_idx = v_idx + 1; u_idx < candidates.size()-1; u_idx++)
        {
            NodeID u = candidates[u_idx];
            assert(u != v);
            if (v_neighbors_set.get(u)) continue; // look for non adjacent nodes
	        if (is_reduced(u, br_alg)) continue;
            if (br_alg->deg(u) + deg_v > config.subgraph_node_limit) continue; //subgraph too large
            if (br_alg->deg(u) < 3) continue; // use other reduction 
            get_neighborhood_set(u, br_alg, u_neighbors_set);

	        //find third heavy vertex (not adjacent to  u and v)
	        for (size_t w_idx = u_idx + 1; w_idx < candidates.size(); w_idx++)
            {
                NodeID w = candidates[w_idx];
                assert(u != w);
                assert(v != w);
                if (v_neighbors_set.get(w)) continue; // look for non adjacent nodes
                if (u_neighbors_set.get(w)) continue; // look for non adjacent nodes
                if (is_reduced(w, br_alg)) continue;
                if (br_alg->deg(w) + deg_v + br_alg->deg(u) > config.subgraph_node_limit) continue; //subgraph too large
                if (br_alg->deg(w) < 3) continue; // use other reduction 
                if (weights[w] + weights[u] + weights[v] < weights[common_neighbor]) continue;

                if (is_heavy_set(v, v_neighbors_set, u, u_neighbors_set, w, br_alg)) {
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
    reverse_mapping.set_size(status.n);
    graph_nodes.clear();
    graph_nodes_set.clear();
    graph_access subgraph;

    sized_vector<NodeID> w_neighbors_vec(status.graph[w].size());
    get_neighborhood_vector(v, br_alg,v_neighbors_vec);
    get_neighborhood_vector(u, br_alg,u_neighbors_vec);
    get_neighborhood_vector(w, br_alg,w_neighbors_vec);


    assert(!status.graph.adjacent(v,u) && "ERROR: heavy_set_reduction::is_heavy_set: v and u must be not adjacent");
    assert(!status.graph.adjacent(v,w) && "ERROR: heavy_set_reduction::is_heavy_set: v and w must be not adjacent");
    assert(!status.graph.adjacent(u,w) && "ERROR: heavy_set_reduction::is_heavy_set: u and w must be not adjacent");
    assert(!is_reduced(v, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: v must not be reduced");
    assert(!is_reduced(u, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: u must not be reduced");
    assert(!is_reduced(w, br_alg) && "ERROR: heavy_set_reduction::is_heavy_set: w must not be reduced");

    for (auto n : v_neighbors_vec) {
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }
    for (auto n : u_neighbors_vec) {
        if (graph_nodes_set.get(n)) continue; // only add nodes once
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }
    for (auto n : w_neighbors_vec) {
        if (graph_nodes_set.get(n)) continue; // only add nodes once
        graph_nodes.push_back(n);
        graph_nodes_set.add(n);
    }

    assert(graph_nodes.size() > 2 && "ERROR: heavy_set_reduction::is_heavy_set: graph_nodes must have at least 3 nodes");

    std::vector<NodeWeight> MWIS_weights(8,0);

    // compute MWIS in N(v) + N(u) + N(w):
    if (!solve_induced_subgraph_from_set(MWIS_weights[v_combination::ooo], subgraph, br_alg, graph_nodes, graph_nodes_set, reverse_mapping)) return false;
    MWIS_weights[v_combination::uvw] = weights[u] + weights[v] + weights[w];

    if (std::min(std::min(weights[v], weights[u]), weights[w]) >= MWIS_weights[v_combination::ooo]) {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        br_alg->set(w, IS_status::included);
        return true;
    } else if (MWIS_weights[v_combination::uvw] < MWIS_weights[v_combination::ooo]) {
        return false;
    }
    bool original_heavy_set = false;
    // compute MWIS[2] in N(v)+N(w)\N(u): (included u)
    unset_weights(subgraph, u_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uoo], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::uoo] > weights[v] + weights[w]) return false;
    MWIS_weights[v_combination::uoo] += weights[u];

    // compute MWIS[3] in N(v)\N(u)+N(w): (included u and w)
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uow], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::uow] > weights[v]) return false;
    MWIS_weights[v_combination::uow] += weights[u] + weights[w];

    // compute MWIS[4] in N(v) + N(u) \ N(w): (included w )
    set_weights(subgraph, u_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::oow], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::oow] > weights[u] + weights[v]) return false;
    MWIS_weights[v_combination::oow] += weights[w];

    // compute MWIS[5] in N(u)+N(w)\N(v): (included v)
    set_weights(subgraph, w_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, v_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::ovo], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::ovo] > weights[u] + weights[w]) return false;
    MWIS_weights[v_combination::ovo] += weights[v];

    // compute MWIS in N(u)\N(v)+N(w): (included v and w)
    unset_weights(subgraph, w_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::ovw], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::ovw] > weights[u]) return false;
    MWIS_weights[v_combination::ovw] += weights[v] + weights[w];

    // compute MWIS in N(w)\N(v)+N(u): 
    set_weights(subgraph, w_neighbors_vec, reverse_mapping, weights);
    unset_weights(subgraph, v_neighbors_vec, reverse_mapping);
    if (!solve_graph(MWIS_weights[v_combination::uvo], subgraph, config)) return false;
    if (original_heavy_set && MWIS_weights[v_combination::uvo] > weights[w]) return false;
    MWIS_weights[v_combination::uvo] += weights[v] + weights[u];
    size_t oldn = status.remaining_nodes;

    if (original_heavy_set) {
        br_alg->set(v, IS_status::included);
        br_alg->set(u, IS_status::included);
        br_alg->set(w, IS_status::included);
    }
    else {
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
                if (check_u_combination(MWIS_weights)) br_alg->set(u, IS_status::included);
                if (check_v_combination(MWIS_weights)) br_alg->set(v, IS_status::included);
                break;
            case v_combination::uow:
                if (check_u_combination(MWIS_weights)) br_alg->set(u, IS_status::included);
                if (check_w_combination(MWIS_weights)) br_alg->set(w, IS_status::included);
                break;
            case v_combination::ovw:
                if (check_v_combination(MWIS_weights)) br_alg->set(v, IS_status::included);
                if (check_w_combination(MWIS_weights)) br_alg->set(w, IS_status::included);
                break;
            case v_combination::uoo:
                if (check_u_combination(MWIS_weights)) br_alg->set(u, IS_status::included);
                break;
            case v_combination::ovo:
                if (check_v_combination(MWIS_weights)) br_alg->set(v, IS_status::included);
                break;
            case v_combination::oow:
                if (check_w_combination(MWIS_weights)) br_alg->set(w, IS_status::included);
                break;
            default:
                break;
        }
    }

    return oldn != status.remaining_nodes;

}
void heavy_set3_reduction::unset_weights(graph_access& graph, sized_vector<NodeID>& nodes, sized_vector<NodeID>& reverse_mapping) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], 0);
    }
}
void heavy_set3_reduction::set_weights(graph_access& graph, sized_vector<NodeID>& nodes, sized_vector<NodeID>& reverse_mapping, std::vector<NodeWeight>& weights) {
    for (NodeID n : nodes) {
        graph.setNodeWeight(reverse_mapping[n], weights[n]);
    }
}
bool heavy_set3_reduction::check_u_combination(std::vector<NodeWeight>& MWIS_weights) {
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::oow], MWIS_weights[v_combination::ovo], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::uow], MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::uoo]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight) { return weight >= best_weight_excluding; });
}
bool heavy_set3_reduction::check_v_combination(std::vector<NodeWeight>& MWIS_weights) {
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::uow], MWIS_weights[v_combination::uoo], MWIS_weights[v_combination::oow], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::ovo]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight) { return weight >= best_weight_excluding; });
}
bool heavy_set3_reduction::check_w_combination(std::vector<NodeWeight>& MWIS_weights) {
    NodeWeight best_weight_excluding = std::max({MWIS_weights[v_combination::uvo], MWIS_weights[v_combination::uoo], MWIS_weights[v_combination::ovo], MWIS_weights[v_combination::ooo]});
    std::vector<NodeWeight> weights_including = {MWIS_weights[v_combination::uvw], MWIS_weights[v_combination::uow], MWIS_weights[v_combination::ovw], MWIS_weights[v_combination::oow]};
    return std::all_of(weights_including.begin(), weights_including.end(), [&](NodeWeight weight) { return weight >= best_weight_excluding; });
}

bool generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_generalized_fold) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
        if(try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;
        reduce_vertex(br_alg, v);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool generalized_fold_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    if (br_alg->config.disable_generalized_fold) return false;
    if (br_alg->deg(v) <= 1) return false;
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
	auto& MWIS_set = br_alg->set_2;
	auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;

    get_neighborhood_set(v, br_alg, neighbors_set);
    get_neighborhood_vector(v, br_alg, neighbors);

    NodeID max_neighbor = get_max_weight_neighbor(v, br_alg);
    NodeWeight max_neighbor_weight = status.weights[max_neighbor];

    if (status.graph[v].size() > br_alg->config.subgraph_node_limit) {
        return false;
    }

    // compute MWIS in N(v)

    NodeWeight MWIS_weight = 0;
    bool solved_exact = solve_induced_subgraph_from_set(MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping, true); 
    if (!solved_exact)  {
        reduction_time += br_alg->reduction_timer.elapsed();
        return false;
    }

    if (status.weights[v] >= MWIS_weight) {
        if (br_alg->config.generate_training_data) {
            reduction_time += br_alg->reduction_timer.elapsed();
            return false; // is heavy vertex reduction
        }
        // same as in generalized_neighborhood_reduction
        br_alg->set(v, IS_status::included);
        return oldn != status.remaining_nodes;
    }
    NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();

    MWIS_set.clear();

    forall_nodes(neighborhood_graph, node) {
        if (neighborhood_graph.getPartitionIndex(node) == 1) {
            const NodeID neighbor = neighbors[node];
            MWIS_set.add(neighbor);

            if (status.weights[neighbor] < min_MWIS_neighbor_weight)
                min_MWIS_neighbor_weight = status.weights[neighbor];
        }
    } endfor

    if (status.weights[v] < MWIS_weight - min_MWIS_neighbor_weight) {
        // multiple IS exist that have bigger weight than v
        return false;
    }

    bool check_failed = false;

    // check that no other IS in N(v) exists with weight greater than v
    for (const NodeID neighbor : status.graph[v]) {
        if (!MWIS_set.get(neighbor))
            continue;

        neighbors.remove(std::find(neighbors.begin(), neighbors.end(), neighbor));
        neighbors_set.remove(neighbor);

        NodeWeight MWIS_weight = 0;
        bool solved_exact = solve_induced_subgraph_from_set(MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping);
        if (!solved_exact) {
            check_failed = true;
        } else if (MWIS_weight >= status.weights[v]) {
            check_failed = true;
        }
        // neighborhood_br_alg.ch.enable_cout();

        neighbors.push_back(neighbor);
        neighbors_set.add(neighbor);

        if (check_failed)
            break;
    }

    if (!check_failed) {
        fold(br_alg, v, MWIS_set, MWIS_weight);
        return oldn != status.remaining_nodes;
    }

    auto& neighborhood_intersection_set = MWIS_set;
    bool remove_node;

    // we can't fold but we can possibly remove some neighbors of v
    do {
        for (const NodeID node : status.graph[v]) {
            neighborhood_intersection_set.clear();

            for (const NodeID neighbor : status.graph[node]) {
                if (neighbors_set.get(neighbor)) {
                    neighborhood_intersection_set.add(neighbor);
                }
            }

            // "force" node into an IS (= remove node and its neighbors from N(v) and compute MWIS in remaining N(v))
            neighbors.remove(std::find(neighbors.begin(), neighbors.end(), node));
            neighbors_set.remove(node);

            for (const NodeID neighbor : status.graph[node]) {
                if (neighborhood_intersection_set.get(neighbor)) {
                    neighbors.remove(std::find(neighbors.begin(), neighbors.end(), neighbor));
                    neighbors_set.remove(neighbor);
                }
            }

            NodeWeight MWIS_weight = 0;
            bool solved_exact = solve_induced_subgraph_from_set(MWIS_weight, neighborhood_graph, br_alg, neighbors, neighbors_set, reverse_mapping);
            if (!solved_exact) {
                remove_node = false;
            } else {
                // if the weight of every MWIS in N(v) which contains "node" is smaller than w(v) then we can remove "node"
                remove_node = MWIS_weight + status.weights[node] <= status.weights[v];
            }

            for (const NodeID neighbor : status.graph[node]) {
                if (neighborhood_intersection_set.get(neighbor)) {
                    neighbors.push_back(neighbor);
                    neighbors_set.add(neighbor);
                }
            }

            if (remove_node) {
                br_alg->set(node, IS_status::excluded);
                break; // break and restart loop because set(..) modifies the range which we currently iterate
            }

            neighbors.push_back(node);
            neighbors_set.add(node);
        }
    } while (remove_node);


	return oldn != status.remaining_nodes;
}
void generalized_fold_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main_node, fast_set& MWIS_set, NodeWeight MWIS_weight) {
	auto& status = br_alg->status;

	restore_vec.emplace_back();
	restore_data& data = restore_vec.back();
	data.main_weight = status.weights[main_node];
	data.MWIS_weight = MWIS_weight;

	auto& nodes = data.nodes;
	nodes.main = main_node;

	// temporary copy for iteration
	data.main_neighbor_list = status.graph[main_node];

	for (auto neighbor : data.main_neighbor_list) {
		if (MWIS_set.get(neighbor))
			nodes.MWIS.push_back(neighbor);
		else
			br_alg->set(neighbor, IS_status::excluded);
	}

	// reverse order because of later "restore_edge_and_replace"
	for (int i = nodes.MWIS.size() - 1; i >= 1; i--) {
		br_alg->set(nodes.MWIS[i], IS_status::folded, false);
	}

	br_alg->set(nodes.MWIS[0], IS_status::folded, true);

	data.main_neighbor_list = status.graph[main_node];

	// "move" weight into redu offset
	status.reduction_offset += data.main_weight;

	status.weights[nodes.main] = MWIS_weight - data.main_weight;

	std::vector<NodeID> new_neighbors;
	auto& neighbors = MWIS_set;
	neighbors.clear();
	neighbors.add(main_node);

	for (NodeID MWIS_node : nodes.MWIS) {
		std::vector<NodeID> node_vec;

		for (auto neighbor : status.graph[MWIS_node]) {
			if (neighbors.add(neighbor)) {
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
void generalized_fold_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	// is "restored" in following loops
	status.graph.hide_node(data.nodes.main);
	status.graph[data.nodes.main] = std::move(data.main_neighbor_list);

	for (size_t i = 0; i < data.nodes.MWIS.size(); i++) {
		br_alg->unset(data.nodes.MWIS[i]);

		for (auto neighbor : data.MWIS_node_vecs[i]) {
			status.graph.replace_last_restored_edge(neighbor, data.nodes.MWIS[i]);
		}
	}

	status.weights[data.nodes.main] = data.main_weight;
	status.reduction_offset -= data.main_weight;

	restore_vec.pop_back();
}
void generalized_fold_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto nodes = restore_vec.back().nodes;
	auto MWIS_weight = restore_vec.back().MWIS_weight;
	auto main_status = status.node_status[nodes.main];
	restore(br_alg);

	if (main_status == IS_status::included) {
		status.node_status[nodes.main] = IS_status::excluded;

		for (auto node : nodes.MWIS) {
			status.node_status[node] = IS_status::included;
		}

		status.is_weight += MWIS_weight;
	} else {
		status.node_status[nodes.main] = IS_status::included;

		for (auto node : nodes.MWIS) {
			status.node_status[node] = IS_status::excluded;
		}

		status.is_weight += status.weights[nodes.main];
	}
}

template<typename struction_type, reduction_type type, int vertex_increase>
bool iterative_struction<struction_type, type, vertex_increase>::reduce(branch_and_reduce_algorithm* br_alg) {
    auto &status = br_alg->status;
    br_alg->reduction_timer.restart();
    NodeID oldn = br_alg->status.remaining_nodes;
    bool applied = false;
    for_each_changed_vertex(br_alg, [&](NodeID v) {
        if (reduce_vertex(br_alg, v))
            applied = true;
    });

    reduced_nodes += oldn-br_alg->status.remaining_nodes;
    reduction_time += br_alg->reduction_timer.elapsed();
    return applied;
}
template<typename struction_type, reduction_type type, int vertex_increase>
bool iterative_struction<struction_type, type, vertex_increase>::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    bool applied = false;
    if (br_alg->deg(v) > br_alg->config.struction_degree || !s.reduce(br_alg, v, s.removed_vertices(br_alg, v) + vertex_increase)) 
        return false;

    br_alg->status.folded_stack.push_back(get_reduction_type());
    applied = true;

    return applied;
}

template class iterative_struction<extended_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<true>, reduction_type::struction_plateau, 0>;

template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::reduce(branch_and_reduce_algorithm *br_alg) {
    this->br_alg = br_alg;
    auto &status = br_alg->status;
    br_alg->reduction_timer.restart();
    init_blow_up_phase();

    while (!is_done() && clean_up_queue()) {
        NodeID n = queue.maxElement();
        Gain key = denoise(queue.maxValue());
        size_t set_limit_by_key = f.set_limit(br_alg, n, key);
        size_t plain_set_limit = br_alg->config.set_limit;
        if (f.set_estimate(br_alg, n, key) > plain_set_limit)
            break;

        queue.deleteMax();
        if (s.reduce(br_alg, n, std::min(set_limit_by_key, plain_set_limit))) {
            //struction was successfully executed
            if (br_alg->config.backtrack_style != ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE && update_set.add(n))
                update_list.emplace_back(n, key);
            status.folded_stack.push_back(get_reduction_type());
            ++blow_ups;
            //"blow up" actually can reduce remaining nodes if multiple blow ups per phase are enabled.
            br_alg->min_kernel = std::min(status.remaining_nodes, br_alg->min_kernel);
        } else if (plain_set_limit <= set_limit_by_key) {
            break;
        } else {
            //reinsert with tighter bound
            update_queue_by_key(n, f.key_by_set_estimate(br_alg, n, set_limit_by_key + 1));
        }
    }
    reduction_time += br_alg->reduction_timer.elapsed();
    return blow_ups != 0;
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::update_queue(NodeID n) {
    if (br_alg->deg(n) > br_alg->config.struction_degree) return;
    update_queue_by_key(n, f.key(br_alg, n));
}

template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::clean_up_queue() {
    auto &status = br_alg->status;
    update_set.resize(status.n);

    //Update candidate set
    for (NodeID n : marker.next)
        update_queue(n);
    marker.clear_next();

    while (queue.size()) {
        NodeID n = queue.maxElement();
        if (status.node_status[n] == IS_status::not_set && br_alg->deg(n) <= br_alg->config.struction_degree) return true;
        if (update_set.add(n))
            update_list.emplace_back(n, denoise(queue.maxValue()));
        queue.deleteMax();
    }
    return false;
}

template<typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::is_done() {
    auto &config = br_alg->config;
    auto &status = br_alg->status;
    size_t cur_kernel = status.remaining_nodes;
    size_t max_kernel = config.phase_blow_up_factor * phase_start_kernel;
    return max_kernel < cur_kernel || blow_ups == config.phase_blow_ups;
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::restore(branch_and_reduce_algorithm *br_alg) {
    s.restore(br_alg);
    restored = true;
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::reset(branch_and_reduce_algorithm* br_alg, size_t comp_size) {
    restored = false;
    queue.clear();
}

template<typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::init_blow_up_phase() {
    auto &status = br_alg->status;

    update_set.resize(status.n);
    blow_ups = 0;
    if (restored) {
        for (NodeID n : added_list)
            if (queue.contains(n))
                queue.deleteNode(n);
        for (auto &e : update_list)
            update_queue_by_key(e.first, e.second);
    } else {
        for_each_changed_vertex(br_alg, [&](NodeID n) {
            update_queue(n);
        });
    }
    phase_start_kernel = status.remaining_nodes;
    update_list.clear();
    added_list.clear();
    update_set.clear();

    restored = false;
}

template class blow_up_struction<DegreeKey>;
template class blow_up_struction<IncreaseKey>;
template class blow_up_struction<ApproximateIncreaseKey>;
template class blow_up_struction<RandomKey>;


bool path_reduction::reduce(branch_and_reduce_algorithm *br_alg) {
    
    this->br_alg = br_alg;
    auto &status = br_alg->status;
    size_t old_n = status.remaining_nodes;

    br_alg->buffers[0].clear();
    br_alg->buffers[1].clear();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        reduce_vertex(br_alg, v);
    });

    while (reduce_degree_one_node() || reduce_path());
    reduced_nodes += (old_n - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
    return old_n != status.remaining_nodes;
}
bool path_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) {

    auto& status = br_alg->status;
    size_t old_n = status.remaining_nodes;
    enqueue_node<false>(v);
    reduced_nodes += (old_n - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
    return old_n != status.remaining_nodes;
}
void path_reduction::restore(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_path) return;
    auto &status = br_alg->status;
    auto &restore_data = restore_vec.back();
    auto &path = restore_data.path;
    NodeID v = path[restore_data.start], w = path[restore_data.end - 1];
    for (auto &n_w : restore_data.node_weights)
        status.weights[n_w.first] = n_w.second;
    for (size_t i = restore_data.start; i < restore_data.end; ++i) {
        status.node_status[path[i]] = IS_status::not_set;
        status.remaining_nodes++;
    }
    if (restore_data.relink) {
        status.graph.relink_directed(v, w, path[restore_data.start + 1]);
        status.graph.relink_directed(w, v, path[restore_data.end - 2]);
    } else {
        status.graph.restore_node(w);
        status.graph.restore_node(v);
    }
    status.reduction_offset -= restore_data.offset;
    restore_vec.pop_back();
}
void path_reduction::apply(branch_and_reduce_algorithm* br_alg) {
    auto &restore_data = restore_vec.back();
    auto path = restore_data.path;
    restore(br_alg);
    apply(path);
}
bool path_reduction::dequeue_node(sized_vector<NodeID> &queue, NodeID &n, size_t degree) {
    for (;;) {
        if (queue.empty()) return false;
        n = queue.back();
        queue.pop_back();
        if (br_alg->status.node_status[n] == IS_status::not_set && br_alg->deg(n) == degree) return true;
    }
}
bool path_reduction::reduce_degree_one_node() {
    NodeID n;
    if (!dequeue_node(br_alg->buffers[0], n, 1)) return false;

    auto &status = br_alg->status;

    NodeID u = status.graph[n][0];
    if (status.weights[n] >= status.weights[u]) {
        status.node_status[n] = IS_status::included;
        status.node_status[u] = IS_status::excluded;
        status.graph.hide_node(n);
        status.graph.hide_node(u);
        status.remaining_nodes -= 2;
        status.is_weight += status.weights[n];
        status.modified_stack.push_back(n);
        status.modified_stack.push_back(u);
        for (NodeID v : status.graph[u])
            enqueue_node(v);
    } else {
        restore_vec.emplace_back();
        auto &restore_data = restore_vec.back();

        br_alg->set(n, IS_status::folded);
        add_reduction_offset(status.weights[n], restore_data);
        reassign_weight(u, status.weights[u] - status.weights[n], restore_data);
        enqueue_node(u);
    }
    return true;
}
bool path_reduction::reduce_path() {
    NodeID n;
    if (!dequeue_node(br_alg->buffers[1], n, 2)) return false;

    auto &status = br_alg->status;
    auto &path = br_alg->buffers[2];
    find_max_path(n, path);
    find_MIS_on_path(path);

    NodeID v = path.front(), w = path.back();
    restore_vec.emplace_back();
    auto &restore_data = restore_vec.back();
    if (v == w) {
        NodeWeight w_e = w_e_e;
        NodeWeight w_i = w_i_i - status.weights[v];
        if (br_alg->deg(v) != 2 && w_i > w_e) {
            //fold
            add_reduction_offset(w_e, restore_data);
            fold_path(1, path.size() - 1, restore_data);
            reassign_weight(v, w_i - w_e, restore_data);
            enqueue_node(v);
        } else {
            //include/exclude v
            fold_path(0, path.size() - 1, restore_data);
            status.node_status[v] = w_i > w_e ? IS_status::folded : IS_status::excluded; //set node state to folded instead of included to be consistent with mis_weight
            add_reduction_offset(std::max(w_i, w_e), restore_data);
        }
    } else {
        bool connected = are_connected(v, w);
        if (connected || w_i_i <= w_i_e || w_i_i <= w_e_i) {
            bool keep_v = w_i_e > w_e_e;
            bool keep_w = w_e_i > w_e_e;

            add_reduction_offset(w_e_e, restore_data);
            fold_path(keep_v, path.size() - keep_w,restore_data);
            reassign_weight(v, w_i_e - w_e_e, restore_data);
            reassign_weight(w,  w_e_i - w_e_e, restore_data);
            if (!connected && keep_v && keep_w)
                reconnect(v, w, restore_data);
        } else {
            int gap = static_cast<int>(w_e_e - w_e_i - w_i_e + w_i_i);
            NodeID c = path[1];
            if (gap >= 0) {
                if (path.size() <= 3) return true;
                add_reduction_offset(w_e_e - gap,restore_data);
                fold_path(2, path.size() - 1,restore_data);
                reassign_weight(v, w_i_i - w_e_i, restore_data);
                reassign_weight(c, gap, restore_data);
                reassign_weight(w, w_i_i - w_i_e, restore_data);
                reconnect(c, w, restore_data);
            } else {
                if (path.size() <= 4) return true;
                NodeID d = path[path.size() - 2];
                add_reduction_offset(w_e_e + gap, restore_data);
                fold_path(2, path.size() - 2,restore_data);
                reassign_weight(v, w_i_e - w_e_e, restore_data);
                reassign_weight(c, -gap, restore_data);
                reassign_weight(d, -gap, restore_data);
                reassign_weight(w, w_e_i - w_e_e, restore_data);
                reconnect(c, d, restore_data);
            }
        }
        enqueue_node(v);
        enqueue_node(w);
    }

    restore_vec.back().path.resize(path.size());
    for (NodeID n : path)
        restore_vec.back().path.push_back(n);
    return true;
}
void path_reduction::reassign_weight(NodeID n, NodeWeight w, restore_data &data) {
    br_alg->status.weights[n] = w;
    data.node_weights.emplace_back(n, w);
}
void path_reduction::reconnect(NodeID v, NodeID w, restore_data &data) {
    auto &status = br_alg->status;
    status.graph.add_edge_undirected(v, w);
    data.relink = true;
}
void path_reduction::fold_path(size_t start, size_t end, restore_data &data) {
    auto &path = br_alg->buffers[2];
    auto &status = br_alg->status;
    for (size_t i = start; i < end; ++i) {
        status.node_status[path[i]] = IS_status::folded;
        status.remaining_nodes--;
    }
    status.graph.hide_node(path[start]);
    status.graph.hide_node(path[end - 1]);
    br_alg->add_next_level_neighborhood(path[start]);
    br_alg->add_next_level_neighborhood(path[end - 1]);
    status.modified_stack.push_back(path[end - 1]);
    data.start = start;
    data.end = end;
}
void path_reduction::add_reduction_offset(size_t offset, restore_data &data) {
    br_alg->status.reduction_offset += offset;
    data.offset = offset;
}
bool path_reduction::are_connected(NodeID v, NodeID w) const {
    auto &neighbors = br_alg->status.graph[v];
    return std::find(neighbors.begin(), neighbors.end(), w) != neighbors.end();
}
void path_reduction::find_MIS_on_path(sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    NodeWeight w_i = status.weights[path[1]], w_e = 0;
    find_MIS_on_path(w_i, w_e, path);
    w_e_i = w_i;
    w_e_e = w_e;
    w_e = status.weights[path[0]];
    w_i = 0;
    find_MIS_on_path(w_i, w_e, path);
    w_i_i = w_i;
    w_i_e = w_e;
}
template<bool track_choices>
void path_reduction::find_MIS_on_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    for (size_t i = 2; i < path.size(); ++i) {
        if (track_choices)
            br_alg->bool_buffer[i - 1] = w_i > w_e;
        NodeWeight next_e = std::max(w_e, w_i);
        w_i = w_e + status.weights[path[i]];
        w_e = next_e;
    }
}
template<bool track_choices>
void path_reduction::find_MIS_on_deg_1_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    w_i = status.weights[path[0]];
    w_e = 0;
    NodeWeight next_w_e;
    for (size_t i = 1; i < path.size(); ++i) {
        if (track_choices)
            br_alg->bool_buffer[i - 1] = w_i > w_e;
        next_w_e = std::max(w_e, w_i);
        w_i = w_e + status.weights[path[i]];
        w_e = next_w_e;
    }
}
void path_reduction::find_max_path(NodeID n, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    path.clear();
    path.push_back(n);

    NodeID last;
    NodeID next;
    for (NodeID current : status.graph[n]) {
        last = n;
        std::reverse(path.begin(), path.end());
        path.push_back(current);
        while (next_node_on_path(current, last, n, next)) {
            path.push_back(next);
            last = current;
            current = next;
        }
        if (path.front() == path.back())
            break;
    }
}
void path_reduction::find_max_deg_1_path(NodeID n, sized_vector<NodeID> &path) {
    auto &status = br_alg->status;

    NodeID last = n;
    NodeID current = status.graph[n][0];
    NodeID next;

    path.clear();
    path.push_back(last);
    path.push_back(current);
    while (next_node_on_path(current, last, n, next)) {
        path.push_back(next);
        last = current;
        current = next;
    }
}
bool path_reduction::next_node_on_path(NodeID current, NodeID last, NodeID first, NodeID &next) {
    auto &status = br_alg->status;
    if (br_alg->deg(current) != 2 || current == first)
        return false;

    auto &neighbors = status.graph[current];
    next = neighbors[0] != last ? neighbors[0] : neighbors[1];
    return true;
}
template<bool add_global>
void path_reduction::enqueue_node(NodeID n) {
    if (br_alg->status.node_status[n] != IS_status::not_set)
        return;
    if (br_alg->deg(n) == 1)
        br_alg->buffers[0].push_back(n);
    /*else if (br_alg->deg(n) == 2)
        br_alg->buffers[1].push_back(n);*/
    else if (add_global)
        br_alg->add_next_level_node(n);
}
void path_reduction::apply(sized_vector<NodeID> &path) {
    auto &status = br_alg->status;
    NodeID v = path.front(), w = path.back();
    NodeWeight w_e = status.node_status[v] != IS_status::excluded ? status.weights[v] : 0;
    NodeWeight w_i = status.node_status[v] != IS_status::excluded ? 0 : status.weights[path[1]];
    find_MIS_on_path<true>(w_i, w_e, path);

    auto &choices = br_alg->bool_buffer;
    size_t i = path.size() - 2;
    bool include_v_i = status.node_status[w] == IS_status::included ? 0 : choices[i];
    for (; i >= 1; --i) {
        NodeID n = path[i];
        if (include_v_i) {
            status.node_status[n] = IS_status::included;
            status.is_weight += status.weights[n];
            include_v_i = false;
        } else {
            status.node_status[n] = IS_status::excluded;
            include_v_i = choices[i - 1];
        }
    }
}


template<reduction_type type, int new_nodes>
reduction_ptr make_iterative_struction(const ReductionConfig &config, size_t n) {
    const auto s = config.struction_type;
    if (s == Struction_Type::ORIGINAL) return reduction_ptr(new iterative_struction<original_struction<false>,type,new_nodes>(n));
    if (s == Struction_Type::MODIFIED) return reduction_ptr(new iterative_struction<original_struction<true>,type,new_nodes>(n));
    if (s == Struction_Type::EXTENDED_REDUCED) return reduction_ptr(new iterative_struction<extended_struction<true>,type,new_nodes>(n));
    return reduction_ptr(new iterative_struction<extended_struction<false>,type,new_nodes>(n));
};

reduction_ptr make_decreasing_struction(const ReductionConfig &config, size_t n) {
    return make_iterative_struction<reduction_type::struction_decrease, -1>(config, n);
};

reduction_ptr make_plateau_struction(const ReductionConfig &config, size_t n) {
    return make_iterative_struction<reduction_type::struction_plateau, 0>(config, n);
};

reduction_ptr make_increasing_struction(const ReductionConfig &config, size_t n) {
    if (config.key_type == Key_Type::RANDOM) return reduction_ptr(new blow_up_struction<RandomKey>(config, n));
    if (config.key_type == Key_Type::DEGREE) return reduction_ptr(new blow_up_struction<DegreeKey>(config, n));
    if (config.key_type == Key_Type::INCREASE) return reduction_ptr(new blow_up_struction<IncreaseKey>(config, n));
    return reduction_ptr(new blow_up_struction<ApproximateIncreaseKey>(config, n));
};
