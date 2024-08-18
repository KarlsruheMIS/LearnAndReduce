
#include "high_degree_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "graph_extractor.h"
#include "strongly_connected_components.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool high_degree_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_high_degree) return false;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif

	auto& status = br_alg->status;
    auto& config = br_alg->config;
	auto& excluded_nodes = br_alg->buffers[0];
    size_t oldn = status.remaining_nodes;
    NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
    graph_access globalG;
    auto& reverse_mapping = br_alg->buffers[1];
    reverse_mapping.resize(status.remaining_nodes);
	br_alg->build_graph_access(globalG, reverse_mapping);

	std::vector<int> comp_map(status.remaining_nodes, 0);
	size_t comp_count = strongly_connected_components().strong_components(globalG, comp_map);


	for (size_t node = 0; node < status.remaining_nodes; node++)
	{
		globalG.setPartitionIndex(node, comp_map[node]);
	}

	std::vector<size_t> comp_size(comp_count, 0);
	for (auto comp_id : comp_map)
	{
		comp_size[comp_id]++;
	}

	std::vector<size_t> comp_idx(comp_count);
	std::iota(comp_idx.begin(), comp_idx.end(), 0);

	std::sort(comp_idx.begin(), comp_idx.end(), [&comp_size](const size_t lhs, const size_t rhs)
			  { return comp_size[lhs] < comp_size[rhs]; });

	graph_extractor extractor;
    NodeWeight MWIS_weight = 0;

	for (size_t i : comp_idx)
	{
		if (br_alg->t.elapsed() > config.time_limit)
		{
            #ifdef REDUCTION_INFO
                reduced_nodes += (oldn - status.remaining_nodes);
                reduction_time += br_alg->reduction_timer.elapsed();
            #endif
            return oldn != status.remaining_nodes;
		}

		br_alg->local_mapping.clear();
		graph_access G;
		extractor.extract_block(globalG, G, i, br_alg->local_mapping);
		br_alg->local_graph = &G;

        excluded_nodes.clear();
        if (G.number_of_nodes() < config.subgraph_node_limit)
        {
            bool solved_exact = solve_graph(MWIS_weight, G, br_alg->config, true); 
            if (solved_exact) {
                forall_nodes(G, node) {
                    if (G.getPartitionIndex(node) == 1) {
                        br_alg->set(reverse_mapping[br_alg->local_mapping[node]], IS_status::included);
                    }
                } endfor
                continue;
            }
        }
        // compute local search bound
        fold_nodes MWIS_nodes;
        excluded_nodes.clear();
        NodeWeight best_weight_so_far = br_alg->run_ils(config, G, br_alg->buffers[0], config.max_swaps);
        bool best_is_ils_solution = true;

	    std::vector<NodeID> node_order(G.number_of_nodes());
	    std::iota(node_order.begin(), node_order.end(), 0);

		// include highest degree vertex
		NodeID v = std::max_element(node_order.begin(), node_order.end(), [&](const NodeID lhs, const NodeID rhs)
			  { return br_alg->deg(lhs) < br_alg->deg(rhs) || ( br_alg->deg(lhs) == br_alg->deg(rhs) && status.weights[lhs] > status.weights[rhs] ); })[0];
            
        // store component nodes
        std::vector<bool> component_nodes(status.n, false);
        for (NodeID node = 0; node < G.number_of_nodes(); node++)
            component_nodes[reverse_mapping[br_alg->local_mapping[node]]] = true;

        while (G.number_of_nodes() - G.getNodeDegree(v) - excluded_nodes.size() < 2*config.subgraph_node_limit)
        {

            if (status.node_status[reverse_mapping[br_alg->local_mapping[v]]] == IS_status::not_set)
            {
                MWIS_weight = 0;
                if (!reduce_in_component(br_alg, v, G, MWIS_weight)) break;

                if (MWIS_weight > best_weight_so_far || (best_is_ils_solution && MWIS_weight == best_weight_so_far))
                {
                    best_weight_so_far = MWIS_weight;
                    best_is_ils_solution = false;
                    MWIS_nodes.MWIS.clear();

                    // get MWIS nodes
                    MWIS_nodes.main = reverse_mapping[br_alg->local_mapping[v]];
                    forall_nodes(G, node) {
                        if (G.getPartitionIndex(node) == 1) {
                            MWIS_nodes.MWIS.push_back(reverse_mapping[br_alg->local_mapping[node]]);
                        }
                    } endfor

                }

                excluded_nodes.push_back(v);

                // exclude nodes by setting weight to 0
                for (NodeID node : excluded_nodes)
                    G.setNodeWeight(v, 0);
            }

            // Remove the last element from node_order
            auto include_node_index = std::find(node_order.begin(), node_order.end(), v);
            std::iter_swap(include_node_index, node_order.end() - 1);
            node_order.pop_back();
		    v = std::max_element(node_order.begin(), node_order.end(), [&](const NodeID lhs, const NodeID rhs)
			        {return br_alg->deg(lhs) < br_alg->deg(rhs) || ( br_alg->deg(lhs) == br_alg->deg(rhs) && status.weights[lhs] > status.weights[rhs] ); })[0];

        }
        if (!best_is_ils_solution)
        {
            std::vector<NodeID> global_exclude_nodes;
            global_exclude_nodes.reserve(excluded_nodes.size());
            for (NodeID node : excluded_nodes)
            {
                if (reverse_mapping[br_alg->local_mapping[node]] != MWIS_nodes.main)
                    global_exclude_nodes.push_back(reverse_mapping[br_alg->local_mapping[node]]);
            }
            MWIS_nodes.is_component = comp_count > 1;
            fold(br_alg, MWIS_nodes, best_weight_so_far, global_exclude_nodes, component_nodes); 
        } else {
            for (auto node : excluded_nodes)
                br_alg->set(reverse_mapping[br_alg->local_mapping[node]], IS_status::excluded);
        }
    }

    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
	return oldn != status.remaining_nodes;
}

bool high_degree_reduction::reduce_in_component(branch_and_reduce_algorithm* br_alg, NodeID v, graph_access& G, NodeWeight& MWIS_weight) {
	auto& status = br_alg->status;
    auto& neighbors = br_alg->buffers[0];
	auto& neighbors_set = br_alg->set_1;
    auto& reverse_mapping = br_alg->buffers[1];
	size_t oldn = status.remaining_nodes;

    // include vertex v in MWIS by setting all neighbor weights to 0
    forall_out_edges(G, e, v) {
        NodeID neighbor = G.getEdgeTarget(e);
        G.setNodeWeight(neighbor, 0);
    } endfor

    // compute MWIS in G-N[v]
    bool solved_exact = solve_graph(MWIS_weight, G, br_alg->config, true); 

    //  restore include vertex v in MWIS
    forall_out_edges(G, e, v) {
        NodeID neighbor = G.getEdgeTarget(e);
        G.setNodeWeight(neighbor, status.weights[reverse_mapping[br_alg->local_mapping[neighbor]]]);
    } endfor

	return solved_exact;
}

void high_degree_reduction::fold(branch_and_reduce_algorithm* br_alg, fold_nodes& nodes, NodeWeight MWIS_weight, std::vector<NodeID>& excluded_nodes, std::vector<bool>& component_nodes) {
	auto& status = br_alg->status;

	restore_vec.emplace_back();
	restore_data& data = restore_vec.back();
	data.main_weight = status.weights[nodes.main];
	data.MWIS_weight = MWIS_weight;
    data.component_nodes = component_nodes;
    data.nodes = nodes;

	br_alg->set(nodes.main, IS_status::folded, true);

	status.folded_stack.push_back(get_reduction_type());
    for (NodeID node : excluded_nodes)
	    br_alg->set(node, IS_status::excluded);

}
void high_degree_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();
    br_alg->unset(data.nodes.main);

	restore_vec.pop_back();
}
void high_degree_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& global_status = br_alg->global_status;
	auto data = restore_vec.back();
	auto MWIS_weight = data.MWIS_weight;
	restore(br_alg);

    NodeWeight current_componnent_is_weight = 0;
    if (!data.nodes.is_component)
    {
        current_componnent_is_weight = status.is_weight;
    }
    // get component nodes for current solution weight
    for (NodeID node = 0; node < status.n; node++)
    {
        if (data.component_nodes[node])
        {
            if (status.node_status[node] == IS_status::included)
                current_componnent_is_weight += status.weights[node];
        }
    }
    
    // std::cout << "current weight: " << current_componnent_is_weight << " MWIS weight: " << MWIS_weight << std::endl;    
	if (MWIS_weight > current_componnent_is_weight) {
		for (auto node : data.nodes.MWIS) {
			status.node_status[node] = IS_status::included;
            for (NodeID neighbor : status.graph[node])
            {
                status.node_status[neighbor] = IS_status::excluded;
            }
		}

		status.is_weight += MWIS_weight - current_componnent_is_weight;
		global_status.is_weight += MWIS_weight - current_componnent_is_weight;
	}
    else
		status.node_status[data.nodes.main] = IS_status::excluded;
}

bool high_degree_reduction::solve_graph(NodeWeight& solution, graph_access& graph, ReductionConfig &config, bool apply_solution) {
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
    
    // setup config for solving graph in high degree reduction
    auto c = config;
    c.disable_heuristic_exclude = true;
    c.disable_heuristic_include = true;
    c.use_partition_cover = false;
    // c.time_limit = graph.number_of_nodes() / 10.0;
    c.time_limit = config.reduction_time_limit*0.05;
    c.max_swaps = 1000;

    branch_and_reduce_algorithm solver(graph, c, true);
	solver.ch.disable_cout();

    if (!solver.run_branch_reduce(std::numeric_limits<NodeWeight>::max())) {
	    solver.ch.enable_cout();
        std::cerr << "timeout during solve_graph in high degree reduction" << std::endl;
        return false;
    }
    if (apply_solution) {
        solver.apply_branch_reduce_solution(graph);
    }
	solver.ch.enable_cout();
    solution = solver.get_current_is_weight();
    return true;
}