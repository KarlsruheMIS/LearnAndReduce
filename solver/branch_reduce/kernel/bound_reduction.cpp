
#include "bound_reduction.h"
#include "branch_and_reduce_algorithm.h"
#include "graph_extractor.h"
#include "strongly_connected_components.h"
#include "partition_cover.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool bound_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_bound_reduction) return false;
    auto& status = br_alg->status;
    size_t oldn = br_alg->status.remaining_nodes;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif

    auto& config = br_alg->config;
    NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
    partition_solution_weights.resize(k);
    partition_mappings.resize(k);
    auto& lb_solution_set = br_alg->set_1;
    auto& neighbors_set = br_alg->set_2;
    auto& ub_solution_set = br_alg->double_set;
    std::vector<NodeID> vertices_to_include;
    graph_access globalG;
    auto& reverse_mapping = br_alg->buffers[0];
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

	for (size_t j : comp_idx)
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
		extractor.extract_block(globalG, G, j, br_alg->local_mapping);
		// br_alg->local_graph = &G;

        // compute global lower bound
        NodeWeight lb = br_alg->run_ils(config, G, br_alg->buffers[1], config.max_swaps);
        lb_solution_set.clear();
        forall_nodes(G, node)
        {
            if (G.getPartitionIndex(node) == 1) {
                lb_solution_set.add(node);
            }
        } endfor

        // compute partition for upper bound
	    partition_cover ub_solver(k);
	    ub_solver.create_partition(G, config);
	    // set original weights to solve
	    forall_nodes(G, node) {
	    	G.setNodeWeight(node, status.weights[reverse_mapping[br_alg->local_mapping[node]]]);
	    } endfor
        build_component_graphs(G, config);

        // solve partitions for upper bound
        NodeWeight ub = 0;
        bool skip = false;
        ub_solution_set.clear();
        for (size_t i = 0; i < k; i++)
        {
	        partition_solution_weights[i] = solve_partition(i, config, true);
            for (NodeID node = 0; node < partition_graphs[i]->number_of_nodes(); node++) {
                if (partition_graphs[i]->getPartitionIndex(node) == 1) {
                    NodeID component_node = partition_mappings[i][node];
                    ub_solution_set.add(component_node);
                }
            }
            if (partition_solution_weights[i] == std::numeric_limits<NodeWeight>::max())
            {
                clear_partition_graphs();
                skip = true;
                break;
            }
            ub += partition_solution_weights[i];
        }
        if (skip) continue;
        assert(ub >= lb);
        printf("lb: %lu, ub: %lu\n", lb, ub);
        if (ub == lb) 
        {
            // apply lb solution
	        forall_nodes(G, node) {
	        	if (lb_solution_set.get(node))
            		br_alg->set(reverse_mapping[br_alg->local_mapping[node]], IS_status::included);
	        } endfor
            clear_partition_graphs();
            continue; 
        }

        for (size_t i = 0; i < k; i++)
        {
            auto& graph = partition_graphs[i];
            std::vector<NodeID> vertices_included;
            std::vector<NodeID> vertices_excluded;
            for (NodeID node = 0; node < graph->number_of_nodes(); node++)
            {
                if (ub_solution_set.get(partition_mappings[i][node]))
                    vertices_included.push_back(node);
                else
                    vertices_excluded.push_back(node);
            }
            std::sort(vertices_included.begin(), vertices_included.end(), [&graph](const NodeID lhs, const NodeID rhs)
                      { return graph->getNodeWeight(lhs) > graph->getNodeWeight(rhs); });
            std::sort(vertices_excluded.begin(), vertices_excluded.end(), [&graph](const NodeID lhs, const NodeID rhs)
                      { return graph->getNodeWeight(lhs) < graph->getNodeWeight(rhs); });

            // for (auto node : vertices_included)
            // {
            //     NodeID global_node = reverse_mapping[br_alg->local_mapping[partition_mappings[i][node]]];
            //     // start excluding one node
            //     if (ub - graph->getNodeWeight(node) <= lb) {
            //         NodeWeight new_part_ub = recompute_partition_weight_excluding_node(i, node, config); 
            //         if (new_part_ub == std::numeric_limits<NodeWeight>::max()) continue;
            //         NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
            //         printf("excluding lb: %lu, new_ub: %lu\n", lb, new_ub);
            //         if (new_ub < lb || (new_ub == lb && lb_solution_set.get(partition_mappings[i][node])))
            //         {
            //             br_alg->set(global_node, IS_status::included);
            //             continue;
            //         }
            //         // try excluding two nodes from the upper bound solution:
            //         // if successfull, we can exclude their common neighbors
            //         std::vector<NodeID> nodes = {node};
            //         for (NodeID second_node : vertices_included)
            //         {
            //             if (second_node == node) continue;
            //             nodes.push_back(second_node);
            //             NodeWeight new_part_ub = recompute_partition_weight_excluding_nodes(i, nodes, config); 
            //             if (new_part_ub == std::numeric_limits<NodeWeight>::max()) {
            //                 nodes.pop_back();
            //                 continue;
            //             }
            //             nodes.pop_back();
            //             NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
            //             printf("excluding lb: %lu, new_ub: %lu\n", lb, new_ub);
            //             if (new_ub < lb || (new_ub == lb && lb_solution_set.get(partition_mappings[i][node])))
            //             {
            //                 // exclude common neighbors
            //                 neighbors_set.clear();
            //                 for (EdgeID e = graph->get_first_edge(node); e < graph->get_first_invalid_edge(node); e++)
            //                     neighbors_set.add(graph->getEdgeTarget(e));

            //                 for (EdgeID e = graph->get_first_edge(second_node); e < graph->get_first_invalid_edge(second_node); e++)
            //                 {
            //                     NodeID neighbor = graph->getEdgeTarget(e);
            //                     if (neighbors_set.get(neighbor))
            //                         br_alg->set(reverse_mapping[br_alg->local_mapping[partition_mappings[i][neighbor]]], IS_status::excluded);
            //                 }
            //             }
            //         }
            //     } else break;
            // }

            for (NodeID node : vertices_excluded)
            {
                NodeID global_node = reverse_mapping[br_alg->local_mapping[partition_mappings[i][node]]];
                assert(status.node_status[global_node] == IS_status::not_set);
                    
                // recompute partition weight for nodes beeing included 
                NodeWeight potential_decrease = 0;
                for (EdgeID e = graph->get_first_edge(node); e < graph->get_first_invalid_edge(node); e++)
                {
                    NodeID neighbor = graph->getEdgeTarget(e);
                    if (ub_solution_set.get(partition_mappings[i][neighbor]))
                        potential_decrease += graph->getNodeWeight(neighbor);
                }
                // printf("lb: %d potential decrease: %lu weight %ld ub: %ld\n",lb,  potential_decrease, partition_graphs[i]->getNodeWeight(node),ub);
                if (ub + graph->getNodeWeight(node) - potential_decrease > lb) continue;
                
                NodeWeight new_part_ub = recompute_partition_weight_including_node(i, node, config); 

                if (new_part_ub == std::numeric_limits<NodeWeight>::max()) continue;
                NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
                printf("including lb: %lu, new_ub: %lu\n", lb, new_ub);
                if (new_ub < lb || (new_ub == lb && !lb_solution_set.get(partition_mappings[i][node])))
                {
                    br_alg->set(global_node, IS_status::excluded);
                    // config.disable_bound_reduction = true;
                    // return true;
                    partition_graphs[i]->setNodeWeight(node, 0);
                }
            }
        }
    }

    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - br_alg->status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
    clear_partition_graphs();
    br_alg->config.disable_bound_reduction = true;
	return oldn != br_alg->status.remaining_nodes;
}

inline bool bound_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    return false;
}

void bound_reduction::build_component_graphs(graph_access &G, ReductionConfig &config)
{
    graph_extractor extractor;
    if (partition_graphs.size() == 0) {
        for (size_t i = 0; i < k; i++)
        {
            graph_access *graph = new graph_access();
            partition_mappings[i].clear();
            extractor.extract_block(G, *graph, i, partition_mappings[i]);
            partition_graphs.push_back(graph);
        }
    } else {
        for (size_t i = 0; i < k; i++)
        {
            partition_mappings[i].clear();
            extractor.extract_block(G, *partition_graphs[i], i, partition_mappings[i]);
        }
    }
}

NodeWeight bound_reduction::solve_partition(size_t i, ReductionConfig &config, bool apply_solution)
{   
    NodeWeight solution_weight = 0;
    double time_limit = config.time_limit; 
    config.time_limit = 20;
    branch_and_reduce_algorithm solver(*partition_graphs[i], config, true);
    solver.ch.disable_cout();
    if (solver.run_branch_reduce()) {
            solution_weight = solver.get_current_is_weight();
            if (apply_solution) solver.apply_branch_reduce_solution(*partition_graphs[i]);
            solver.ch.enable_cout();
    } else {
            solver.ch.enable_cout();
            std::cout << "component not solved optimally" << std::endl;
            solution_weight = std::numeric_limits<NodeWeight>::max();
    }
    config.time_limit = time_limit;
    return solution_weight;
}

NodeWeight bound_reduction::recompute_partition_weight_including_node(size_t i, NodeID node, ReductionConfig& config)
{
    // NodeWeight new_weight = 1;
    // NodeWeight original_node_weight = partition_graphs[i]->getNodeWeight(node);
    std::vector<NodeWeight> original_neighbor_weights;
    for (EdgeID e = partition_graphs[i]->get_first_edge(node); e < partition_graphs[i]->get_first_invalid_edge(node); e++)
    {
        NodeID neighbor = partition_graphs[i]->getEdgeTarget(e);
        // new_weight += partition_graphs[i]->getNodeWeight(neighbor);
        original_neighbor_weights.push_back(partition_graphs[i]->getNodeWeight(neighbor));
        partition_graphs[i]->setNodeWeight(neighbor,0);
    }
    // partition_graphs[i]->setNodeWeight(node, new_weight);
    assert(partition_graphs[i]->getPartitionIndex(node) == 0);
    NodeWeight solution_weight = solve_partition(i, config);
    NodeID weight_index = 0;
    for (EdgeID e = partition_graphs[i]->get_first_edge(node); e < partition_graphs[i]->get_first_invalid_edge(node); e++)
    {
        NodeID neighbor = partition_graphs[i]->getEdgeTarget(e);
        partition_graphs[i]->setNodeWeight(neighbor,original_neighbor_weights[weight_index++]);
    }
 
    // partition_graphs[i]->setNodeWeight(node, original_node_weight);
    return solution_weight; // -new_weight + original_node_weight;
}

NodeWeight bound_reduction::recompute_partition_weight_excluding_node(size_t i, NodeID node, ReductionConfig& config)
{
    NodeWeight original_node_weight = partition_graphs[i]->getNodeWeight(node);
    partition_graphs[i]->setNodeWeight(node, 0);
    assert(partition_graphs[i]->getPartitionIndex(node) == 1);
    NodeWeight solution_weight = solve_partition(i, config);
    partition_graphs[i]->setNodeWeight(node, original_node_weight);
    return solution_weight;
}

NodeWeight bound_reduction::recompute_partition_weight_excluding_nodes(size_t i, std::vector<NodeID>& nodes, ReductionConfig& config)
{
    std::vector<NodeWeight> original_node_weights;
    for (auto node : nodes)
    {
        original_node_weights.push_back(partition_graphs[i]->getNodeWeight(node));
        partition_graphs[i]->setNodeWeight(node, 0);
    }
    NodeWeight solution_weight = solve_partition(i, config);
    for (size_t i = 0; i < nodes.size(); i++)
    {
        partition_graphs[i]->setNodeWeight(nodes[i], original_node_weights[i]);
    }
    return solution_weight;
}

void bound_reduction::clear_partition_graphs()
{
    for (size_t i = 0; i < k; i++)
    {
        delete partition_graphs[i];
        partition_graphs[i] = nullptr;
    }
    partition_graphs.clear();
}