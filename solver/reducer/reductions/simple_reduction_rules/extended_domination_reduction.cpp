#include "extended_domination_reduction.h"

#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool extended_domination_reduction::reduce(reduce_algorithm* br_alg) {
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
inline bool extended_domination_reduction::reduce_vertex(reduce_algorithm* br_alg, NodeID v) {
    if(br_alg->config.disable_domination) return false;
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;
	#ifndef generate_training_data
        NodeWeight neighbor_weights = this->get_neighborhood_weight(v, br_alg);
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
void extended_domination_reduction::fold(reduce_algorithm* br_alg, NodeID v, NodeID neighbor) {

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
void extended_domination_reduction::restore(reduce_algorithm* br_alg) {
    auto& status = br_alg->status;
    auto& data = restore_vec.back();

    status.weights[data.v] += data.neighborWeight;
    status.graph.add_edge_undirected(data.v, data.neighbor);
    restore_vec.pop_back();
}
void extended_domination_reduction::apply(reduce_algorithm* br_alg) {
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
int extended_domination_reduction::generate_data(reduce_algorithm* br_alg, NodeID v, std::vector<NodeID>& label) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->set_1;
	size_t oldn = status.remaining_nodes;

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
                label.push_back(v);
	            return true;
            }
        }
    }
    
	return false;
}