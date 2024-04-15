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

        if (v < status.n && status.node_status[v] == IS_status::not_set) {
                f(v);
        }
    }
}
NodeID general_reduction::get_max_weight_neighbor(NodeID v, branch_and_reduce_algorithm* br_alg) { auto& status = br_alg->status;
    NodeID max_neighbor = v;
    NodeWeight max_weight = 0;

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
    auto& status = br_alg->status;
    if (status.weights[v] >= neighbors_weight) {
        br_alg->set(v, IS_status::included);
        return true;
    }
    return false;
}
NodeWeight general_reduction::solve_induced_subgraph_from_set(graph_access& graph, ReductionConfig &config, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& nodes_vec, const fast_set& nodes_set, sized_vector<NodeID>& reverse_mapping, bool apply_solution) {
    br_alg->build_induced_subgraph(graph, nodes_vec, nodes_set, reverse_mapping);
    return solve_graph(graph, config, apply_solution);
}
NodeWeight general_reduction::solve_induced_neighborhood_subgraph(graph_access& neighborhood_graph, ReductionConfig &config, branch_and_reduce_algorithm* br_alg, NodeID v, bool apply_solution) {
    br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
    return solve_graph(neighborhood_graph, config, apply_solution);
}
NodeWeight general_reduction::solve_graph(graph_access& graph, ReductionConfig &config, bool apply_solution) {
    if (graph.number_of_nodes() == 0 ) return 0;
    if (graph.number_of_edges() == 0 ) {
        NodeWeight solution = 0;
        forall_nodes(graph, node) {
            if (graph.getNodeWeight(node)>0)
            {
                graph.setPartitionIndex(node, 1);
                solution += graph.getNodeWeight(node);
            }
        } endfor
        return solution;
    }
    branch_and_reduce_algorithm solver(graph, config, true);
    if (!solver.run_branch_reduce()) {
        std::cerr << "%br_call time out" << std::endl;
        return 0;
    }
    if (apply_solution) {
        solver.apply_branch_reduce_solution(graph);
    }
    return solver.get_current_is_weight();
}
bool general_reduction::is_reduced(NodeID v, branch_and_reduce_algorithm* br_alg) {
    return br_alg->status.node_status[v] != IS_status::not_set;
}

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_neighborhood) return false;
    br_alg->reduction_timer.restart();

    auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for_each_changed_vertex(br_alg, [this, &br_alg, &status](NodeID v) {
        NodeWeight neighbor_weights = this->get_neighborhood_weight(v, br_alg);
        try_neighborhood_reduction(v, br_alg, neighbor_weights);
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool fold1_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_fold1) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for_each_changed_vertex(br_alg, [this, &br_alg, &status](NodeID v) {
        if (br_alg->deg(v) == 0) br_alg->set(v, IS_status::included);
        if (br_alg->deg(v) == 1) 
        {
            NodeID neighbor = status.graph[v][0];
            assert(status.node_status[neighbor] == IS_status::not_set);
            if (status.weights[neighbor] > status.weights[v]) 
                   this->fold(br_alg, {v, neighbor});
            else 
                   br_alg->set(v, IS_status::included, true);
        }
    });

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
void fold1_reduction::fold(branch_and_reduce_algorithm* br_alg, const fold_nodes& nodes) {

    auto& status = br_alg->status;

	restore_vec.push_back({nodes, status.weights[nodes.deg1_node], status.graph[nodes.fold_node]});
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

bool fold2_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_fold2) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	for_each_changed_vertex(br_alg, [this, &br_alg, &status](NodeID v) {
		if (br_alg->deg(v) != 2) return;

        // set bigger and smaller neighbor
        NodeID bigger  = status.graph[v][0];
        NodeID smaller = status.graph[v][1];
        if (try_neighborhood_reduction(v, br_alg, status.weights[bigger] + status.weights[smaller])) return;
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
                if (br_alg->config.disable_triangle) return;
                break;
            }
        }
        
        if (triangle && status.weights[bigger] <= status.weights[v]) {
              br_alg->set(v, IS_status::included);
        } else if (triangle && status.weights[bigger] > status.weights[v] && status.weights[smaller] <= status.weights[v]) {
            this->fold_triangle_mid_weight(br_alg, {v, {bigger, smaller}});
        } else if (triangle) {
            this->fold_triangle_min_weight(br_alg, {v, {bigger, smaller}});
        } else if (status.weights[v] >= status.weights[bigger]) {
            if (br_alg->config.disable_v_shape_max) return;
            this->fold_v_shape_max_weight(br_alg, {v, {bigger, smaller}});
        } else if (status.weights[v] >= status.weights[smaller]) 
        {
            if (br_alg->config.disable_v_shape_mid) return;
            this->fold_v_shape_mid_weight(br_alg, {v, {bigger, smaller}});
        } else {
            assert(status.weights[v] < status.weights[smaller] && "v is not the smallest");
            if (br_alg->config.disable_v_shape_min) return;
            this->fold_v_shape_min_weight(br_alg, {v, {bigger, smaller}});
        }

	});

    if (v_shape_min_count >= status.modified_stack.capacity() - status.n) { // initial cpapcity is 2n (one n for v_shape min)
        assert(status.modified_stack.capacity() >= status.modified_stack.size() && "stack size error");
        assert(status.folded_stack.capacity() >= status.folded_stack.size() && "stack size error");
        std::cout << "resize modified stack " << status.modified_stack.capacity(); 
        status.modified_stack.resize(status.modified_stack.capacity() + v_shape_min_count);
        status.folded_stack.resize(status.folded_stack.capacity() + v_shape_min_count);
        std::cout << " to " << status.modified_stack.capacity();
        std::cout << " v_shape_min_count " << v_shape_min_count;
        std::cout << " used size: " << status.modified_stack.size();
        std::cout << " n " << status.n << std::endl;
    }
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
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

bool single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_basic_se) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

    NodeWeight partial_neighbor_sum = 0;
	for_each_changed_vertex(br_alg, [&](NodeID v) {
        this->get_neighborhood_set(v, br_alg, neighbors);

        for (NodeID neighbor : status.graph[v]) {
            if (status.weights[v] <= status.weights[neighbor]) { // otherwise not applicable to this edge
                // compute w(N(neighbor)\N(v))
                NodeWeight partial_neighbor_sum = 0;
                for (NodeID second_neighbor : status.graph[neighbor]) {
                    if (status.node_status[second_neighbor] == IS_status::not_set && neighbors.add(second_neighbor)) {
                        partial_neighbor_sum += status.weights[second_neighbor];
                        neighbors.remove(second_neighbor);
                    }
                }
             
                // note: weight of v is in partial_neighbor_sum included
                if (partial_neighbor_sum <= status.weights[neighbor]) { 
                    br_alg->set(v, IS_status::excluded);
                    break;
                }
            }
        }
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool extended_single_edge_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_extended_se) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;
    NodeWeight max_neighbor_weight = 0;
	auto& neighbors = br_alg->set_1;
	auto& checked = br_alg->set_2;

	for_each_changed_vertex(br_alg, [&](NodeID v) {
            checked.clear();

            NodeWeight neighbors_weight = this->get_neighborhood_weight(v, br_alg);
            if (this->try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;
            this->get_neighborhood_set(v, br_alg, neighbors);
            NodeID max_neighbor = this->get_max_weight_neighbor(v, br_alg);
			NodeWeight max_neighbor_weight = status.weights[max_neighbor];
            
            while (status.weights[v] >= neighbors_weight - max_neighbor_weight) {

            
                for (NodeID neighbor : status.graph[max_neighbor]) {
                    if (neighbor == v)  continue; 
                    if (status.node_status[neighbor] == IS_status::not_set) {
                        // exclude neighborhood intersection and update neighborhood
                        if (neighbors.get(neighbor)) { 
                            br_alg->set(neighbor, IS_status::excluded);
                            neighbors.remove(neighbor);
                            neighbors_weight -= status.weights[neighbor];
                        }
                    }
                }

                if (this->try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;

                // check if different edge satisfies reduction
                checked.add(max_neighbor);
                max_neighbor_weight = 0;
                for (NodeID neighbor : status.graph[v]) {
                    if (checked.get(neighbor)) continue;
                    if (max_neighbor_weight < status.weights[neighbor]) { 
                        max_neighbor_weight = status.weights[neighbor];
                        max_neighbor = neighbor;
                    }
                }
            }
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool domination_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    // special case of single edge reduction (only used in original versions)
    if(br_alg->config.disable_domination) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
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

            is_subset = true;

            for (NodeID neighbor2 : status.graph[neighbor]) {
                if (!neighbors.get(neighbor2)) {
                    is_subset = false;
                    break;
                }
            }

            if (is_subset && status.weights[neighbor] >= status.weights[v]) {
                br_alg->set(v, IS_status::excluded);
                break;
            }
        }
    });


    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
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
    sized_vector<NodeID> cut_v_included_i(config.cut_vertex_max_component_size);
    sized_vector<NodeID> cut_v_included_e(config.cut_vertex_max_component_size);
    sized_vector<NodeID> cut_v_excluded_i(config.cut_vertex_max_component_size);
    sized_vector<NodeID> cut_v_excluded_e(config.cut_vertex_max_component_size);
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
        assert(cut_component.size() <= config.cut_vertex_max_component_size && "ERROR: cut_vertex_reduction::reduce: component size too large");
        if (cut_component.size() <= 1) continue; //fold1
        // std::cout << "cut_vertex: " << cut_v << " component size: " << cut_component.size()<< std::endl;
        // check if component is correct
        // #ifdef DEBUG
        // fast_set cut_component_set_test(status.n);
        // cut_component_set_test.add(cut_v);
        // for (NodeID node : cut_component) { cut_component_set_test.add(node); }
        // for (NodeID node : cut_component) { 
        //     for (NodeID neighbor : status.graph[node]) {
        //         assert(cut_component_set_test.get(neighbor) && "ERROR: cut_vertex_reduction::reduce: cut_component not correct");
        //     }
        // }
        // #endif

        cut_component_set.clear();
        get_neighborhood_set(cut_v, br_alg, cut_v_neighbor_set);
        for (auto neighbor : cut_component) {
            cut_component_set.add(neighbor);
            visited[neighbor] = true;
        }
        // graph including neighborhood of cut_v
        config.time_limit = std::min(config.reduction_time_limit*0.1, 0.5*br_alg->config.reduction_time_limit - br_alg->t.elapsed());
        reverse_mapping.assign(status.n, status.n);
        br_alg->ch.disable_cout();
        NodeWeight large_cutMWIS_weight = solve_induced_subgraph_from_set(cut_graph, config, br_alg, cut_component, cut_component_set, reverse_mapping, true);
        if (large_cutMWIS_weight == 0) continue;
        br_alg->ch.enable_cout();

        // save solution for later reduce
        cut_v_excluded_e.clear();
        cut_v_excluded_i.clear();
        cut_v_excluded_e.push_back(cut_v);
        for (NodeID node : cut_component)
        {
            assert(reverse_mapping[node] != status.n && "ERROR: cut_vertex_reduction::reduce: node not in reverse_mapping");
            if (cut_graph.getPartitionIndex(reverse_mapping[node]) == 1) cut_v_excluded_i.push_back(node);
            else cut_v_excluded_e.push_back(node);
        }


        // set weight of neighborhood of cut_v to 0, ie solve G-N(cut_v)
        for (NodeID neighbor : status.graph[cut_v])
        {
            if (reverse_mapping[neighbor] == status.n) continue;
            cut_graph.setNodeWeight(reverse_mapping[neighbor], 0);
            cut_graph.setPartitionIndex(reverse_mapping[neighbor], 0);
        } 
        config.time_limit = std::min(config.reduction_time_limit*0.1, br_alg->config.reduction_time_limit - br_alg->t.elapsed());
        br_alg->ch.disable_cout();
        NodeWeight small_cutMWIS_weight = solve_graph(cut_graph, config, true);
        br_alg->ch.enable_cout();
        if (small_cutMWIS_weight == 0) continue;
        if (status.weights[cut_v] + small_cutMWIS_weight <= large_cutMWIS_weight) // cut_vertex is excluded -> directly apply
        {
            for (NodeID node : cut_v_excluded_i)
            {
                br_alg->set(node, IS_status::included);
            }

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

        // #ifdef DEBUG // check vectors are correct
        //print info
        // std::cout << "v_cut: " << cut_v << " neighbors: "  ; for (NodeID node : status.graph[cut_v]) { std::cout << node << " "; } std::cout << std::endl;
        // std::cout << "cut_comp: " ; for (NodeID node : cut_component) { std::cout << node << " "; } std::cout <<std::endl;
        // std::cout << "cut_v_excluded_i: " << cut_v_excluded_i.size() << " : "; for (NodeID node : cut_v_excluded_i) { std::cout << node << " "; } std::cout << std::endl;
        // std::cout << "cut_v_excluded_e: " << cut_v_excluded_e.size() << " : "; for (NodeID node : cut_v_excluded_e) { std::cout << node << " "; } std::cout << std::endl;
        // std::cout << "cut_v_included_i: " << cut_v_included_i.size() << " : "; for (NodeID node : cut_v_included_i) { std::cout << node << " "; } std::cout << std::endl;
        // std::cout << "cut_v_included_e: " << cut_v_included_e.size() << " : "; for (NodeID node : cut_v_included_e) { std::cout << node << " "; } std::cout << std::endl;
        // std::cout << "large_cutMWIS_weight: " << large_cutMWIS_weight ;
        // std::cout << " small_cutMWIS_weight: " << small_cutMWIS_weight ;
        // std::cout << " cut_v weight: " << status.weights[cut_v] << std::endl;

        // fast_set cut_v_neighbors(status.n);
        // for (NodeID node : status.graph[cut_v]) { cut_v_neighbors.add(node); }
        // for (NodeID node : cut_v_included_i) { assert(!cut_v_neighbors.get(node) && "ERROR: cut_v_neighbors must be excluded in this case"); }
        // for (NodeID node : cut_v_excluded_i) { assert(cut_v != node && "ERROR: cut_v must be excluded in this case"); }
        // #endif
            fold_data data = {cut_v, status.weights[cut_v], large_cutMWIS_weight, small_cutMWIS_weight, cut_component};
            fold(br_alg, data, cut_v_included_i, cut_v_included_e, cut_v_excluded_i, cut_v_excluded_e);
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
                // std::cout << "visited: " << local_n << " overall unvisited: " << status.n - std::count(visited.begin(), visited.end(), true) << std::endl;
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
            } else if (smallComponent.size() <= config.cut_vertex_max_component_size) {
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
        if (component.size() > br_alg->config.cut_vertex_max_component_size) {
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

bool clique_neighborhood_reduction_fast::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_clique_neighborhood_fast) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
    auto& neighbors = br_alg->buffers[0];
    auto& neighborhood = br_alg->set_1;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
        if (this->try_neighborhood_reduction(v, br_alg, neighbor_weights)) return;
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

	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool clique_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_clique_neighborhood) return false;
    br_alg->reduction_timer.restart();
	this->br_alg = br_alg;
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        // if (br_alg->deg(v) > 10) return;
        if (partition_into_cliques(v)) {
			br_alg->set(v, IS_status::included);
		}
	});

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}
bool clique_neighborhood_reduction::partition_into_cliques(NodeID v) {
	auto& weights= br_alg->status.weights;
	auto& neighbors_vec = br_alg->buffers[0];
	auto& clique_neighbors_set = br_alg->set_1;
    target_weight = weights[v];

    neighbor_weights = get_neighborhood_weight(v, br_alg);
	if (neighbor_weights <= target_weight)  return true; 
    get_neighborhood_vector(v, br_alg, neighbors_vec);

	// partition neigbors of v into cliques
	NodeID max_neighbor;
	NodeWeight max_neighbor_weight;

	while (neighbors_vec.size() >= 2 && br_alg->config.reduction_time_limit > br_alg->t.elapsed()) {
        max_neighbor_weight = 0;
        clique_neighbors_set.clear();
        for (auto neighbor : neighbors_vec) {
            if (weights[neighbor] > max_neighbor_weight) {
                max_neighbor = neighbor;
                max_neighbor_weight = weights[neighbor];
            }
            clique_neighbors_set.add(neighbor);
        }
		clique_neighbors_set.remove(max_neighbor);

		if (expand_clique(max_neighbor, neighbors_vec, clique_neighbors_set))
			return true;
	}

	return false;
}
bool clique_neighborhood_reduction::expand_clique(NodeID max_neighbor, sized_vector<NodeID>& neighbors_vec, fast_set& clique_neighbors_set) {
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
        if (local_max == status.n) break;
		if (intersection_empty)	break;

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
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	auto& set_1 = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	auto& isolated = br_alg->buffers[1];
	std::vector<NodeID> non_isolated;
	size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID node) {
        get_neighborhood_set(node, br_alg, set_1);
        get_neighborhood_vector(node, br_alg, neighbors);
        set_1.add(node);

        // check if clique
        non_isolated.clear();
        isolated.clear();
        isolated.push_back(node);

        size_t max_isolated_idx = 0;
        weighted_node max_isolated{ node, status.weights[node] };
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

            if (count != neighbors.size()) return;
        }

        // one of "isolated" members has highest weight of clique: Add to IS
        // also handles completely isolated cliques
        if (max_isolated.weight >= max_non_isolated.weight) {
            br_alg->set(max_isolated.node, IS_status::included);
            return;
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
    });


    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
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
    br_alg->reduction_timer.restart();
	auto& weights = br_alg->status.weights;
	auto& graph = br_alg->status.graph;
	auto& remaining_n = br_alg->status.remaining_nodes;
	auto& funnel_set = br_alg->set_1;
	auto& neighbors = br_alg->buffers[0];
	size_t oldn = remaining_n;
    NodeID funnel_neighbor = br_alg->status.n;

    for_each_changed_vertex(br_alg, [&](NodeID node) {
        get_neighborhood_vector(node, br_alg, neighbors);
        funnel_set.clear();
        if (std::any_of(neighbors.begin(), neighbors.end(), [&](NodeID neighbor) { return weights[neighbor] >= weights[node]; })) {
            get_neighborhood_set(node, br_alg, funnel_set);
            funnel_set.add(node);

            if (is_funnel(node, funnel_neighbor, br_alg, funnel_set, neighbors)) {
                fold({node, funnel_neighbor}, funnel_set, br_alg);
            }
        }
    });

    reduced_nodes += (oldn - remaining_n);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != remaining_n;
}
bool funnel_reduction::is_funnel(NodeID node, NodeID& funnel_neighbor, branch_and_reduce_algorithm* br_alg, fast_set& funnel_set, sized_vector<NodeID>& funnel_nodes) {
    auto& status = br_alg->status;
    auto& weights = status.weights;
    auto& graph = status.graph;
    auto& neighbors = br_alg->buffers[1];
    funnel_neighbor = status.n;
    for (NodeID v : funnel_nodes) {
        assert(v != node && "ERROR: funnel_reduction::is_funnel: node must not be in funnel_nodes");
        assert(!is_reduced(v, br_alg) && "ERROR: funnel_reduction::is_funnel: node must be unset");
        if (weights[v] >= weights[node]) {
            funnel_neighbor = v;

            bool skip = false;
            for (auto neighbor :graph[node])
            {
                if (is_reduced(neighbor, br_alg)) continue;
                if (neighbor == funnel_neighbor) continue;
                if (weights[neighbor] > weights[node]) { // check with common neighbors only?
                    skip = true;
                    break;
                }
            }
            if(skip) continue;

            funnel_nodes.remove(funnel_neighbor);
            funnel_set.remove(funnel_neighbor);
            if (is_clique(br_alg, funnel_set, funnel_nodes)) return true; 
            else {
                funnel_nodes.push_back(funnel_neighbor);
                funnel_set.add(funnel_neighbor);
                funnel_neighbor = status.n;
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
    auto& funnel_neighbors = br_alg->set_2; 
    NodeID node = data.node;
    NodeID funnel_neighbor = data.funnel_neighbor;


    funnel_neighbors.clear();
    for (auto neighbor : graph[funnel_neighbor]) {
        if (!is_reduced(neighbor, br_alg)) funnel_neighbors.add(neighbor);
    }

    // common neighbors are not part of the IS in either case 
    std::vector<NodeID> remaining_neighbors;
    remaining_neighbors.reserve(graph[node].size());
    for (auto neighbor : graph[node]) {
        if (is_reduced(neighbor, br_alg)) continue;
        assert(status.node_status[neighbor] == IS_status::not_set && "ERROR: funnel_reduction::fold: common neighbors must be unset");
        if (neighbor == funnel_neighbor) continue;
        if (funnel_neighbors.get(neighbor)) 
            br_alg->set(neighbor, IS_status::excluded);
        else
            remaining_neighbors.push_back(neighbor); 
    }
    remaining_neighbors.resize(remaining_neighbors.size());

    return;
    status.folded_stack.push_back(get_reduction_type());
    status.reduction_offset += weights[node];

    // add additional edges connecting the remaining neighbors with funnel neighors outside the funnel set  
    std::vector<NodeID> outside_funnel_neighbors;
    outside_funnel_neighbors.reserve(status.graph[funnel_neighbor].size());
    for (auto neighbor : status.graph[funnel_neighbor])
    {
        if (is_reduced(neighbor, br_alg)) continue;
        if (!funnel_set.get(neighbor)) 
            outside_funnel_neighbors.push_back(neighbor);
    }

    restore_vec.push_back({data.node, data.funnel_neighbor, outside_funnel_neighbors, remaining_neighbors, {}});
    for (size_t idx = 0; idx < remaining_neighbors.size(); idx++) {
        restore_vec.back().node_vecs.push_back({});
        NodeID neighbor = remaining_neighbors[idx];
        assert(status.node_status[neighbor] == IS_status::not_set && "ERROR: funnel_reduction::fold: remaining neighbors must be unset");
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
void funnel_reduction::restore(branch_and_reduce_algorithm* br_alg) {
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
        status.weights[neighbor] -= status.weights[data.funnel_neighbor];
        status.weights[neighbor] += status.weights[data.node];
    }
    status.reduction_offset -= status.weights[data.node];

	restore_vec.pop_back();
}
void funnel_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& node_status = br_alg->status.node_status;
	auto& is_weight = br_alg->status.is_weight;
	auto& weights = br_alg->status.weights;
	auto& data = restore_vec.back();
	auto& remaining = data.remaining_neighbors;
	auto& outside = data.outside_funnel_neighbors;

    // print remaining neighor node status:
    // std::cout << "remaining neighbors: "; for (auto neighbor : remaining) { std::cout << neighbor << " (" << node_status[neighbor] << ") "; } std::cout << std::endl;
    bool include_node = true;
    if (std::any_of(remaining.begin(), remaining.end(), [&](NodeID neighbor) { return node_status[neighbor] == IS_status::included; })) {
        include_node = false;
    }
    if (include_node)
    { // if no origianl funnel_neighbor neighbors are included, we can still include funnel_neighbor
        if (std::none_of(outside.begin(), outside.end(), [&](NodeID neighbor) { return node_status[neighbor] == IS_status::included; })) {
            include_node = false;
        }
    }

    restore(br_alg);

    if (include_node) {
        node_status[data.node] = IS_status::included;
        node_status[data.funnel_neighbor] = IS_status::excluded;
        is_weight += weights[data.node];
    } else {
        node_status[data.node] = IS_status::excluded;
        node_status[data.funnel_neighbor] = IS_status::included;
        is_weight += weights[data.funnel_neighbor];
    }
    // std::cout << "node: " << data.node << " (" << node_status[data.node] << ") funnel_neighbor: " << data.funnel_neighbor << " (" << node_status[data.funnel_neighbor] << ")\n";
}

bool twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_twin) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;
	size_t oldn = status.remaining_nodes;
	NodeID twin;

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
        if (this->try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;
        get_neighborhood_vector(v, br_alg, neighbors);

        twin_candidates_set.clear();
        bool candidates_empty = true;

        for (NodeID neighbor : status.graph[neighbors[0]]) {
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
                if (twin_candidates_set.get(candidate)) {
                    tmp_set.add(candidate);
                    candidates_empty = false;
                    twin = candidate;
                }
            }

            std::swap(twin_candidates_set, tmp_set);
        }

        if (candidates_empty) return;

        if (status.weights[v] + status.weights[twin] >= neighbors_weight) {
            br_alg->set(v, IS_status::included);
            br_alg->set(twin, IS_status::included);
        } else {
            fold(br_alg, v, twin);
        }
    });


    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
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

bool heavy_set_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    auto config = br_alg->config;
    if (config.disable_heavy_set) return false;
    br_alg->reduction_timer.restart();

	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;
	size_t oldw = status.reduction_offset;

	auto& build_graph_neighbors = br_alg->buffers[0];
	auto& reverse_mapping = br_alg->buffers[1];

	//sets to build different subgraphs
	auto& build_graph_neighbors_set = br_alg->set_1;
	auto& heavy_vertex_neighbors_set = br_alg->set_2;
    fast_set second_heavy_vertex_neighbors_set(status.n);
	graph_access subgraph;

    br_alg->ch.disable_cout();
    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
        if (this->try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;
        NodeID heavy_vertex = get_max_weight_neighbor(v, br_alg);
        if (heavy_vertex == v) return;
		::NodeWeight heavy_vertex_weight = status.weights[heavy_vertex];

		for (NodeID neighbor : status.graph[v]) {
            if (br_alg->deg(neighbor) > config.heavy_set) return; //subgraph too large
		}

        get_neighborhood_set(heavy_vertex, br_alg, heavy_vertex_neighbors_set);
        get_neighborhood_set(heavy_vertex, br_alg, build_graph_neighbors_set);
        get_neighborhood_vector(heavy_vertex, br_alg, build_graph_neighbors);

		//find second heavy vertex (not adjacent to heavy vertex)
        NodeID second_heavy_vertex = v;
		::NodeWeight second_heavy_vertex_weight = 0;
		for (NodeID neighbor : status.graph[v]) {
	        if (status.node_status[neighbor] == IS_status::not_set) {
                if (neighbor == heavy_vertex) continue; // look for different nodes
                if (heavy_vertex_neighbors_set.get(neighbor)) continue; // look for non adjacent nodes
                if (br_alg->deg(neighbor) + br_alg->deg(heavy_vertex) > config.disable_heavy_set) continue; //subgraph too large

    		    if (status.weights[neighbor] > second_heavy_vertex_weight) {
			    	second_heavy_vertex_weight = status.weights[neighbor];
                    second_heavy_vertex = neighbor;
                }
            }
        }

        if (second_heavy_vertex_weight == 0) return;
		second_heavy_vertex_neighbors_set.clear();

		for (NodeID neighbor : status.graph[second_heavy_vertex])
		{
		 	if (status.node_status[neighbor] != IS_status::not_set) continue;
		    second_heavy_vertex_neighbors_set.add(neighbor);
		    build_graph_neighbors.push_back(neighbor);
		    build_graph_neighbors_set.add(neighbor);
        }

        //build neighborhood graph 
		config.time_limit = 8.0 / 10.0;
        bool only_check_single = false;

        //check different cases on the neighborhood graph:
		// case 1)
        // compute MWIS in N(heavy_vertex) \cup N(second_heavy_vertex):
        ::NodeWeight MWIS_weight_case1 = this->solve_induced_subgraph_from_set(subgraph, config, br_alg, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);

        if (second_heavy_vertex_weight >= MWIS_weight_case1) {
            br_alg->set(heavy_vertex, IS_status::included);
            br_alg->set(second_heavy_vertex, IS_status::included);
			return;
        } else if (heavy_vertex_weight + second_heavy_vertex_weight < MWIS_weight_case1)
            only_check_single = true;


		// case 2)
        // compute MWIS in N(heavy_vertex): 
        // first set graph to N(heavy_vertex)
		build_graph_neighbors.clear();
		std::for_each(status.graph[heavy_vertex].begin(), status.graph[heavy_vertex].end(), [&](NodeID neighbor) {
			if (status.node_status[neighbor] == IS_status::not_set) {
				build_graph_neighbors.push_back(neighbor);
			}
		});
        ::NodeWeight MWIS_weight_case2 = this->solve_induced_subgraph_from_set(subgraph, config, br_alg, build_graph_neighbors, heavy_vertex_neighbors_set, reverse_mapping);

        if (heavy_vertex_weight >= MWIS_weight_case2) {
            br_alg->set(heavy_vertex, IS_status::included);
            only_check_single = true;
        }

        if (!only_check_single) {
            // case 3)
            // compute MWIS in N(heavy_vertex)\N(second_heavy_vertex): 

		    build_graph_neighbors.clear();
		    build_graph_neighbors_set.clear();
		    std::for_each(status.graph[heavy_vertex].begin(), status.graph[heavy_vertex].end(), [&](NodeID neighbor) {
			    if (status.node_status[neighbor] == IS_status::not_set && !second_heavy_vertex_neighbors_set.get(neighbor)) {
			    	build_graph_neighbors.push_back(neighbor);
			    	build_graph_neighbors_set.add(neighbor);
			    }
		    });
            ::NodeWeight MWIS_weight_case3 = this->solve_induced_subgraph_from_set(subgraph, config, br_alg, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);
                   
		    // case 4)
            // compute MWIS in N(second_heavy_vertex)\N(heavy_vertex): 
            // first set graph 

			build_graph_neighbors.clear();
			build_graph_neighbors_set.clear();
			std::for_each(status.graph[second_heavy_vertex].begin(), status.graph[second_heavy_vertex].end(), [&](NodeID neighbor) {
				if (status.node_status[neighbor] == IS_status::not_set && !heavy_vertex_neighbors_set.get(neighbor)) {
					build_graph_neighbors.push_back(neighbor);
					build_graph_neighbors_set.add(neighbor);
				}
			});
            ::NodeWeight MWIS_weight_case4 = this->solve_induced_subgraph_from_set(subgraph, config, br_alg, build_graph_neighbors, build_graph_neighbors_set, reverse_mapping);

            if (heavy_vertex_weight + second_heavy_vertex_weight >= MWIS_weight_case1 &&
                heavy_vertex_weight >= MWIS_weight_case3 &&
                second_heavy_vertex_weight >= MWIS_weight_case4) {
                br_alg->set(heavy_vertex, IS_status::included);
                br_alg->set(second_heavy_vertex, IS_status::included);
			    return;
            }
        }

		// case 5)
        // compute MWIS in N(second_heavy_vertex):
        // first set graph 
        if (only_check_single) {
        	build_graph_neighbors.clear();
			std::for_each(status.graph[second_heavy_vertex].begin(), status.graph[second_heavy_vertex].end(), [&](NodeID neighbor) {
				if (status.node_status[neighbor] == IS_status::not_set) {
					build_graph_neighbors.push_back(neighbor);
				}
			});
            ::NodeWeight MWIS_weight_case5 = this->solve_induced_subgraph_from_set(subgraph, config, br_alg, build_graph_neighbors, second_heavy_vertex_neighbors_set, reverse_mapping);

            if (second_heavy_vertex_weight >= MWIS_weight_case5) {
                br_alg->set(second_heavy_vertex, IS_status::included);
            }
        }
    });
    br_alg->ch.enable_cout();

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool generalized_neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_generalized_neighborhood) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	auto config = br_alg->config;
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;

	br_alg->ch.disable_cout();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        NodeWeight neighbors_weight = get_neighborhood_weight(v, br_alg);
        if (try_neighborhood_reduction(v, br_alg, neighbors_weight)) return;

        NodeID max_neighbor = get_max_weight_neighbor(v, br_alg);
        NodeWeight max_neighbor_weight = status.weights[max_neighbor]; 

        //MWIS in N(v) >= max_neighbor_weight > w(v)
        if (status.weights[v] < max_neighbor_weight) return;

        if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) 
            return;
        
        //std::cerr << "%gnr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

        // compute MWIS in N(v)
        config.time_limit = status.graph[v].size() / 10.0;
        br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
        branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

        if (!neighborhood_br_alg.run_branch_reduce()) {
            std::cerr << "%generalized_neighborhood_reduction br_call time out" << std::endl;
            return;
        }

        if (status.weights[v] >= neighborhood_br_alg.get_current_is_weight())
            br_alg->set(v, IS_status::included);
    });

	br_alg->ch.enable_cout();

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
	return oldn != status.remaining_nodes;
}

bool generalized_fold_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_generalized_fold) return false;
    br_alg->reduction_timer.restart();
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
	auto& MWIS_set = br_alg->set_2;
	auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;

	graph_access neighborhood_graph;
	auto config = br_alg->config;

	br_alg->ch.disable_cout();

    for_each_changed_vertex(br_alg, [&](NodeID v) {
        get_neighborhood_set(v, br_alg, neighbors_set);
        get_neighborhood_vector(v, br_alg, neighbors);

        NodeID max_neighbor = get_max_weight_neighbor(v, br_alg);
        NodeWeight max_neighbor_weight = status.weights[max_neighbor];

        if (status.graph[v].size() > branch_and_reduce_algorithm::REDU_RECURSION_LIMIT) {
            return;
        }

        //std::cerr << "%gfr try neighborhood of size " << status.graph[v].size() << " start node: " << v << std::endl;

        // compute MWIS in N(v)
        config.time_limit = status.graph[v].size() / 10.0;

        br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);
        #ifdef DEBUG
        solution_check<graph_access> sc(neighborhood_graph);
        sc.check_graph();
        #endif
        branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

        if (!neighborhood_br_alg.run_branch_reduce()) {
            std::cerr << "%generalized_fold_reduction br_call time out" << std::endl;
            return;
        }

        NodeWeight MWIS_weight = neighborhood_br_alg.get_current_is_weight();
        NodeWeight min_MWIS_neighbor_weight = std::numeric_limits<NodeWeight>::max();

        if (status.weights[v] >= MWIS_weight) {
            // same as in generalized_neighborhood_reduction
            br_alg->set(v, IS_status::included);
            return;
        }

        neighborhood_br_alg.apply_branch_reduce_solution(neighborhood_graph);
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
            return;
        }

        bool check_failed = false;

        // check that no other IS in N(v) exists with weight greater than v
        for (const NodeID neighbor : status.graph[v]) {
            if (!MWIS_set.get(neighbor))
                continue;

            neighbors.remove(std::find(neighbors.begin(), neighbors.end(), neighbor));
            neighbors_set.remove(neighbor);

            br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);
            branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

            if (!neighborhood_br_alg.run_branch_reduce()) {
                std::cerr << "%generalized_fold_reduction br_call loop time out" << std::endl;
                check_failed = true;
            }
            else if (neighborhood_br_alg.get_current_is_weight() >= status.weights[v]) {
                check_failed = true;
            }

            neighbors.push_back(neighbor);
            neighbors_set.add(neighbor);

            if (check_failed)
                break;
        }

        if (!check_failed) {
            fold(br_alg, v, MWIS_set, MWIS_weight);
            return;
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

                br_alg->build_induced_subgraph(neighborhood_graph, neighbors, neighbors_set, reverse_mapping);

                config.time_limit = neighbors.size() / 10.0;
                branch_and_reduce_algorithm neighborhood_br_alg(neighborhood_graph, config, true);

                if (!neighborhood_br_alg.run_branch_reduce()) {
                    std::cerr << "%generalized_fold_reduction br_call loop time out" << std::endl;
                    remove_node = false;
                }
                else {
                    // if the weight of every MWIS in N(v) which contains "node" is smaller than w(v) then we can remove "node"
                    remove_node = neighborhood_br_alg.get_current_is_weight() + status.weights[node] <= status.weights[v];
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
    });

	br_alg->ch.enable_cout();

    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
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
    bool applied = false;
    for_each_changed_vertex(br_alg, [&](NodeID v) {
        if (br_alg->deg(v) > br_alg->config.struction_degree || !s.reduce(br_alg, v, s.removed_vertices(br_alg, v) + vertex_increase)) return;

        status.folded_stack.push_back(get_reduction_type());
        applied = true;
    });

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
        enqueue_node<false>(v);
    });

    while (reduce_degree_one_node() || reduce_path());
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
