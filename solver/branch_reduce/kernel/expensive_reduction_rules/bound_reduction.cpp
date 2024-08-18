
#include "bound_reduction.h"
#include "branch_and_reduce_algorithm.h"
#include "graph_extractor.h"
#include "strongly_connected_components.h"
// #include "partition_cover.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool bound_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
    if (br_alg->config.disable_bound_reduction) return false;
    auto& status = br_alg->status;
    size_t oldn = br_alg->status.remaining_nodes;
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
    #endif

    auto& config = br_alg->config;
// 
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

        reduce_component(br_alg, G);
    }

    #ifdef REDUCTION_INFO
        reduced_nodes += (oldn - br_alg->status.remaining_nodes);
        reduction_time += br_alg->reduction_timer.elapsed();
    #endif
    // if (partition_graphs.size() > 0) clear_partition_graphs();
    br_alg->config.disable_bound_reduction = true;
	return oldn != br_alg->status.remaining_nodes;
}
NodeWeight bound_reduction::get_lower_bound(branch_and_reduce_algorithm* br_alg, graph_access& G, fast_set& lb_solution)
{
    auto& lb_solution_set = br_alg->set_1;
    auto& config = br_alg->config;
    // compute global lower bound
    NodeWeight lb = br_alg->run_ils(config, G, br_alg->buffers[1], config.max_swaps);
    lb_solution_set.clear();
    forall_nodes(G, node)
    {
        if (G.getPartitionIndex(node) == 1) {
            lb_solution_set.add(node);
        }
    } endfor
    return lb;
}

void bound_reduction::reduce_component(branch_and_reduce_algorithm* br_alg, graph_access& G) {
    NodeWeight no_limit = std::numeric_limits<NodeWeight>::max();
    auto& status = br_alg->status;
    auto& config = br_alg->config;
    auto& reverse_mapping = br_alg->buffers[0];
    auto& lb_solution_set = br_alg->set_1;
    auto& ub_solution_set = br_alg->double_set;
    std::vector<NodeID> vertices_to_include;
    std::vector<NodeID> reverse_clique_mapping;

    // compute global lower bound
    NodeWeight lb = get_lower_bound(br_alg, G, lb_solution_set); 

    int n_cliques = 0;
    NodeWeight ub = get_clique_cover_bound(br_alg, n_cliques, ub_solution_set, vertices_to_include, reverse_clique_mapping);
    // NodeWeight ub =  get_partition_bound(br_alg, G, ub_solution_set);

    assert(lb <= ub);
    // bool use_partition = true;
    // bool use_partition = false;
    // bool use_partition = n_cliques > 20;
    // if (use_partition)
    // {
    //     NodeWeight partition_ub = get_partition_bound(br_alg, G, ub_solution_set);
    //     assert(lb <= partition_ub);
    //     if (partition_ub >= ub) use_partition = false;
    //     else ub = partition_ub;
    // }

    printf("lb: %lu, ub: %lu\n", lb, ub);
    assert(lb <= ub);
    if (ub == lb) 
    {
        // apply lb solution
	    forall_nodes(G, node) {
	    	if (lb_solution_set.get(node))
            {
                assert(status.node_status[reverse_mapping[br_alg->local_mapping[node]]] == IS_status::not_set);
        		br_alg->set(reverse_mapping[br_alg->local_mapping[node]], IS_status::included);
            }
	    } endfor
        // if (use_partition) clear_partition_graphs();
        return; 
    }
    // if (use_partition) reduce_by_partition_bound(br_alg, G, lb, ub); 
    else reduce_by_clique_cover_bound(br_alg, G, lb, ub, vertices_to_include, reverse_clique_mapping, n_cliques);
}

// void bound_reduction::reduce_by_partition_bound(branch_and_reduce_algorithm* br_alg, graph_access& G, NodeWeight lb, NodeWeight ub)
// {
//     auto& config = br_alg->config;
//     auto& reverse_mapping = br_alg->buffers[0];
//     auto& lb_solution_set = br_alg->set_1;
//     auto& ub_solution_set = br_alg->double_set;
//     auto& neighbors_set = br_alg->set_2;

//     for (size_t i = 0; i < k; i++)
//     {
//         auto& graph = partition_graphs[i];
//         std::vector<NodeID> vertices_included;
//         std::vector<NodeID> vertices_excluded;
//         for (NodeID node = 0; node < graph->number_of_nodes(); node++)
//         {
//             if (ub_solution_set.get(partition_mappings[i][node]))
//                 vertices_included.push_back(node);
//             else
//                 vertices_excluded.push_back(node);
//         }
//         std::vector<NodeWeight> potential_decrease(graph->number_of_nodes(), 0);
//         for (NodeID node : vertices_excluded)
//         {
//             for (EdgeID e = graph->get_first_edge(node); e < graph->get_first_invalid_edge(node); e++)
//             {
//                 NodeID neighbor = graph->getEdgeTarget(e);
//                 if (ub_solution_set.get(partition_mappings[i][neighbor]))
//                     potential_decrease[node] += graph->getNodeWeight(neighbor);
//             }
//             assert(potential_decrease[node] >= graph->getNodeWeight(node));
//             potential_decrease[node] -= graph->getNodeWeight(node);
//         }

//         std::sort(vertices_included.begin(), vertices_included.end(), [&graph](const NodeID lhs, const NodeID rhs)
//                   { return graph->getNodeWeight(lhs) > graph->getNodeWeight(rhs); });
//         std::sort(vertices_excluded.begin(), vertices_excluded.end(), [&graph, &potential_decrease](const NodeID lhs, const NodeID rhs)
//                   { return  potential_decrease[lhs] > potential_decrease[rhs]; });

//         for (NodeID node : vertices_excluded)
//         {
//             NodeID global_node = reverse_mapping[br_alg->local_mapping[partition_mappings[i][node]]];
//             if (br_alg->status.node_status[global_node] != IS_status::not_set) continue;
//             if (ub - 0.5*potential_decrease[node] > lb) break;
            
//             NodeWeight new_part_ub = recompute_partition_weight_including_node(i, node, config); 

//             if (new_part_ub == std::numeric_limits<NodeWeight>::max()) continue;
//             NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
//             if (new_ub < lb || (new_ub == lb && !lb_solution_set.get(partition_mappings[i][node])))
//             {
//                 printf("excluding lb: %lu, old ub: %d, new_ub: %lu\n", lb, ub, new_ub);
//                 br_alg->set(global_node, IS_status::excluded);
//             }
//         }


//         for (auto node : vertices_included)
//         {
//             NodeID global_node = reverse_mapping[br_alg->local_mapping[partition_mappings[i][node]]];
//             if (br_alg->status.node_status[global_node] != IS_status::not_set) continue;
//             // start excluding one node
//             if (ub - graph->getNodeWeight(node) <= lb) {
//                 NodeWeight new_part_ub = recompute_partition_weight_excluding_node(i, node, config); 
//                 if (new_part_ub == std::numeric_limits<NodeWeight>::max()) continue;
//                 NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
//                     // printf("excluding lb: %lu, old ub: %d new_ub: %lu\n", lb, ub, new_ub);
//                 if (new_ub < lb || (new_ub == lb && lb_solution_set.get(partition_mappings[i][node])))
//                 {
//                     printf("including lb: %lu, old ub: %d new_ub: %lu\n", lb, ub, new_ub);
//                     br_alg->set(global_node, IS_status::included);
//                     continue;
//                 }
//                 // try excluding two nodes from the upper bound solution:
//                 // if successfull, we can exclude their common neighbors
//                 // std::vector<NodeID> nodes = {node};
//                 // for (NodeID second_node : vertices_included)
//                 // {
//                 //     if (second_node == node) continue;
//                 //     nodes.push_back(second_node);
//                 //     NodeWeight new_part_ub = recompute_partition_weight_excluding_nodes(i, nodes, config); 
//                 //     if (new_part_ub == std::numeric_limits<NodeWeight>::max()) {
//                 //         nodes.pop_back();
//                 //         continue;
//                 //     }
//                 //     nodes.pop_back();
//                 //     NodeWeight new_ub = ub - partition_solution_weights[i] + new_part_ub;
//                 //     if (new_ub < lb || (new_ub == lb && lb_solution_set.get(partition_mappings[i][node])))
//                 //     {
//                 //         printf("excluding lb: %lu, new_ub: %lu\n", lb, new_ub);
//                 //         // exclude common neighbors
//                 //         neighbors_set.clear();
//                 //         for (EdgeID e = graph->get_first_edge(node); e < graph->get_first_invalid_edge(node); e++)
//                 //             neighbors_set.add(graph->getEdgeTarget(e));

//                 //         for (EdgeID e = graph->get_first_edge(second_node); e < graph->get_first_invalid_edge(second_node); e++)
//                 //         {
//                 //             NodeID neighbor = graph->getEdgeTarget(e);
//                 //             NodeID global_neighbor = reverse_mapping[br_alg->local_mapping[partition_mappings[i][neighbor]]]; 
//                 //             if (neighbors_set.get(neighbor) && !is_reduced(global_neighbor, br_alg))
//                 //                 br_alg->set(global_neighbor, IS_status::excluded);
//                 //         }
//                 //     }
//                 // }
//             } else break;
//         }

//     }
// }

void bound_reduction::reduce_by_clique_cover_bound(branch_and_reduce_algorithm* br_alg, graph_access& G, NodeWeight lb, NodeWeight ub, std::vector<NodeID>& vertices_to_include, std::vector<NodeID>& reverse_mapping, int n_cliques)
{
    auto& config = br_alg->config;
    auto& status = br_alg->status;
    auto& buffers = br_alg->buffers;
    auto& node_mapping = buffers[1];
	auto& clique_weight = buffers[3];
	auto& clique_sizes = buffers[4];
	auto& clique_second_weight = buffers[6];
	auto& clique = buffers[2];

    // try excluding cover vertices -> clique cover with second highest weight 
    std::cout << "lb " << lb << " ub: " << ub  << " clique cover " << std::endl;
    for (size_t i = 0; i < n_cliques; i++)
    {
        if (lb >= ub - clique_weight[clique[i]] + clique_second_weight[clique[i]])
        {
            std::cout << "lb " << lb << " new_ub: " << ub - clique_weight[clique[i]] + clique_second_weight[clique[i]] << std::endl;
            br_alg->set(reverse_mapping[clique[i]], IS_status::included);
        }
    }
}

// NodeWeight bound_reduction::get_partition_bound(branch_and_reduce_algorithm* br_alg, graph_access& G, fast_set& ub_solution_set)
// {
//     auto& config = br_alg->config;
//     auto& status = br_alg->status;
//     auto& lb_solution_set = br_alg->set_1;
//     partition_solution_weights.resize(k);
//     partition_mappings.resize(k);
//     // compute partition for upper bound
// 	partition_cover ub_solver(k);
// 	ub_solver.create_partition(G, config);

//     build_component_graphs(G, config);

//     NodeWeight ub = 0;
//     ub_solution_set.clear();
//     // solve partitions for upper bound
//     for (size_t i = 0; i < k; i++)
//     {
// 	    partition_solution_weights[i] = solve_partition(i, config, true);
//         for (NodeID node = 0; node < partition_graphs[i]->number_of_nodes(); node++) 
//         {
//             if (partition_graphs[i]->getPartitionIndex(node) == 1) {
//                 NodeID component_node = partition_mappings[i][node];
//                 ub_solution_set.add(component_node);
//             }
//         }
//         if (partition_solution_weights[i] == std::numeric_limits<NodeWeight>::max())
//         {
//             clear_partition_graphs();
//             return std::numeric_limits<NodeWeight>::max();
//         }
//         ub += partition_solution_weights[i];
//     }
//     return ub;
// }

NodeWeight bound_reduction::get_clique_cover_bound(branch_and_reduce_algorithm* br_alg, int& n_cliques, fast_set& solution, std::vector<NodeID>& vertices_to_include, std::vector<NodeID>& reverse_mapping)
{
    auto& status = br_alg->status;
    auto& buffers = br_alg->buffers;
	auto &clique = buffers[2];
	auto &clique_weight = buffers[3];
    vertices_to_include.clear();
    reverse_mapping.resize(br_alg->status.remaining_nodes,0);

    NodeWeight ub = compute_cover_pruning_bound(br_alg, n_cliques, reverse_mapping);

	auto end_iter = std::unique(clique.begin(), clique.end());
    vertices_to_include.reserve(n_cliques);
	for (auto iter = clique.begin(); iter != end_iter; iter++)
    {
        vertices_to_include.push_back(reverse_mapping[*iter]);
		solution.add(reverse_mapping[*iter]);
    }

    return ub;
}

NodeWeight bound_reduction::compute_cover_pruning_bound(branch_and_reduce_algorithm* br_alg, int& n_cliques, std::vector<NodeID>& reverse_mapping)
{
    auto& buffers = br_alg->buffers;
    auto& status = br_alg->status;
	// Gather remaining nodes
	auto &nodes = buffers[7];
	nodes.clear();

	for (size_t node = 0; node < status.n; node++)
	{
		if (status.node_status[node] == IS_status::not_set)
			nodes.emplace_back(node);
	}

	// Sort by descending weight
	// Break ties by degree
	std::sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs)
			  { return (status.weights[lhs] > status.weights[rhs] || (status.weights[lhs] == status.weights[rhs] && br_alg->deg(lhs) > br_alg->deg(rhs))); });

	// Compute node mapping
	NodeID current_node = 0;
	auto &node_mapping = buffers[1];
	node_mapping.resize(status.n);

	for (NodeID node : nodes)
    {
        reverse_mapping[current_node] = node;
		node_mapping[node] = current_node++;
    }

	// Init cliques
	auto &clique = buffers[2];
	auto &clique_weight = buffers[3];
	auto &clique_sizes = buffers[4];
	auto &covered = buffers[5];
	auto &clique_second_weight = buffers[6];

	clique.resize(nodes.size());
	clique_weight.resize(nodes.size());
	clique_second_weight.resize(nodes.size());
	clique_sizes.resize(nodes.size());
	covered.resize(nodes.size());

	std::fill(clique_weight.begin(), clique_weight.end(), 0);
	std::fill(clique_second_weight.begin(), clique_second_weight.end(), 0);
	std::fill(clique_sizes.begin(), clique_sizes.end(), 0);
	std::fill(covered.begin(), covered.end(), false);

	for (NodeID node : nodes)
		clique[node_mapping[node]] = node_mapping[node];

	auto &neigh_clique_sizes = buffers[6];
	neigh_clique_sizes.resize(nodes.size());

	for (NodeID node : nodes)
	{
		NodeID v = node_mapping[node];
		// Find heaviest neighboring clique
		NodeID heaviest_clique = v;
		std::fill(neigh_clique_sizes.begin(), neigh_clique_sizes.end(), 0);

		for (NodeID neighbor : status.graph[node])
		{
			if (status.node_status[neighbor] == IS_status::not_set)
			{
				NodeID w = node_mapping[neighbor];

				if (covered[w])
				{
					NodeID c = clique[w];
					neigh_clique_sizes[c]++;

					if (neigh_clique_sizes[c] == clique_sizes[c] && clique_weight[c] > clique_weight[heaviest_clique])
						heaviest_clique = c;
				}
			}
		}

		// Update clique weights/sizes
		clique[v] = heaviest_clique;
        if (clique_weight[heaviest_clique] <= status.weights[node])
        {
		    clique_second_weight[heaviest_clique] = clique_weight[heaviest_clique];
		    clique_weight[heaviest_clique] =  status.weights[node];
        }
		clique_sizes[heaviest_clique]++;

		// Node is covered
		covered[v] = true;
	}

	// Keep unique cliques
	std::sort(clique.begin(), clique.end());
	auto end_iter = std::unique(clique.begin(), clique.end());

	// Add clique weights to get upper bound
	NodeWeight upper_bound = 0;
	for (auto iter = clique.begin(); iter != end_iter; iter++)
		upper_bound += clique_weight[*iter];

	n_cliques = std::distance(clique.begin(), end_iter);
	return upper_bound;
}

// void bound_reduction::build_component_graphs(graph_access &G, ReductionConfig &config)
// {
//     graph_extractor extractor;
//     if (partition_graphs.size() == 0) {
//         for (size_t i = 0; i < k; i++)
//         {
//             graph_access *graph = new graph_access();
//             partition_mappings[i].clear();
//             extractor.extract_block(G, *graph, i, partition_mappings[i]);
//             partition_graphs.push_back(graph);
//         }
//     } else {
//         for (size_t i = 0; i < k; i++)
//         {
//             partition_mappings[i].clear();
//             extractor.extract_block(G, *partition_graphs[i], i, partition_mappings[i]);
//         }
//     }
// }

// NodeWeight bound_reduction::solve_partition(size_t i, ReductionConfig &config, bool apply_solution)
// {   
//     NodeWeight solution_weight = 0;
//     auto c = config;
//     double time_limit = c.time_limit; 
//     c.time_limit = 20;
//     c.use_partition_cover = false;
//     branch_and_reduce_algorithm solver(*partition_graphs[i], c, true);
//     solver.ch.disable_cout();
//     if (solver.run_branch_reduce()) {
//             solution_weight = solver.get_current_is_weight();
//             if (apply_solution) solver.apply_branch_reduce_solution(*partition_graphs[i]);
//             solver.ch.enable_cout();
//     } else {
//             solver.ch.enable_cout();
//             std::cout << "component not solved optimally" << std::endl;
//             solution_weight = std::numeric_limits<NodeWeight>::max();
//     }
//     return solution_weight;
// }

// NodeWeight bound_reduction::recompute_partition_weight_including_node(size_t i, NodeID node, ReductionConfig& config)
// {
//     std::vector<NodeWeight> original_neighbor_weights;
//     for (EdgeID e = partition_graphs[i]->get_first_edge(node); e < partition_graphs[i]->get_first_invalid_edge(node); e++)
//     {
//         NodeID neighbor = partition_graphs[i]->getEdgeTarget(e);
//         original_neighbor_weights.push_back(partition_graphs[i]->getNodeWeight(neighbor));
//         partition_graphs[i]->setNodeWeight(neighbor,0);
//     }
//     assert(partition_graphs[i]->getPartitionIndex(node) == 0);
//     NodeWeight solution_weight = solve_partition(i, config);
//     NodeID weight_index = 0;
//     for (EdgeID e = partition_graphs[i]->get_first_edge(node); e < partition_graphs[i]->get_first_invalid_edge(node); e++)
//     {
//         NodeID neighbor = partition_graphs[i]->getEdgeTarget(e);
//         partition_graphs[i]->setNodeWeight(neighbor,original_neighbor_weights[weight_index++]);
//     }
 
//     return solution_weight; 
// }

// NodeWeight bound_reduction::recompute_partition_weight_excluding_node(size_t i, NodeID node, ReductionConfig& config)
// {
//     NodeWeight original_node_weight = partition_graphs[i]->getNodeWeight(node);
//     partition_graphs[i]->setNodeWeight(node, 0);
//     assert(partition_graphs[i]->getPartitionIndex(node) == 1);
//     NodeWeight solution_weight = solve_partition(i, config);
//     partition_graphs[i]->setNodeWeight(node, original_node_weight);
//     return solution_weight;
// }

// NodeWeight bound_reduction::recompute_partition_weight_excluding_nodes(size_t i, std::vector<NodeID>& nodes, ReductionConfig& config)
// {
//     std::vector<NodeWeight> original_node_weights;
//     for (auto node : nodes)
//     {
//         original_node_weights.push_back(partition_graphs[i]->getNodeWeight(node));
//         partition_graphs[i]->setNodeWeight(node, 0);
//     }
//     NodeWeight solution_weight = solve_partition(i, config);
//     for (size_t node = 0; node < nodes.size(); node++)
//     {
//         partition_graphs[i]->setNodeWeight(nodes[node], original_node_weights[node]);
//     }
//     return solution_weight;
// }

// void bound_reduction::clear_partition_graphs()
// {
//     for (size_t i = 0; i < k; i++)
//     {
//         delete partition_graphs[i];
//         partition_graphs[i] = nullptr;
//     }
//     partition_graphs.clear();
// }

inline bool bound_reduction::reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) {
    return false;
}
