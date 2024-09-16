/******************************************************************************
 * Copyright (C) 2015-2017 Darren Strash <strash@kit.edu>
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

#include "branch_and_reduce_algorithm.h"
#include "reductions.h"
#include "ils/ils.h"
#include "hils/hils.h"
#include "ils/local_search.h"
#include "strongly_connected_components.h"
#include "graph_extractor.h"
#include "ReduceAndPeel.h"
#include "solution_check.h"
#include "struction_log.h"
#include "LRConv.h"
#include "struction_reductions.h"
// #include "partition_cover.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <climits>

constexpr NodeID branch_and_reduce_algorithm::BRANCHING_TOKEN;
constexpr NodeID branch_and_reduce_algorithm::MODIFIED_TOKEN;

branch_and_reduce_algorithm::branch_and_reduce_algorithm(graph_access &G, const ReductionConfig &config, bool called_from_fold)
	: config(config), global_status(G), set_1(global_status.n), set_2(global_status.n), double_set(global_status.n * 2),
	  buffers(4), zero_vec(global_status.n, 0), gnn(called_from_fold ? 0 : G.number_of_nodes())
{
	buffers[0].reserve(global_status.n);
	buffers[1].reserve(global_status.n);
	buffers[2].reserve(global_status.n);
	buffers[3].reserve(global_status.n);
	bool_buffer.reserve(global_status.n);
	weight_buffer.reserve(global_status.n);

	if (!config.disable_neighborhood)
		global_status.transformations.emplace_back(new neighborhood_reduction(global_status.n));
	if (!config.disable_fold1)
		global_status.transformations.emplace_back(new fold1_reduction(global_status.n));
	if (!config.disable_fold2)
		global_status.transformations.emplace_back(new fold2_reduction(global_status.n));
	if (!config.disable_clique)
		global_status.transformations.emplace_back(new clique_reduction(global_status.n));
	if (!config.disable_domination)
		global_status.transformations.emplace_back(new domination_reduction(global_status.n));
	if (!config.disable_basic_se)
		global_status.transformations.emplace_back(new single_edge_reduction(global_status.n));
	if (!config.disable_extended_se)
		global_status.transformations.emplace_back(new extended_single_edge_reduction(global_status.n));
	if (!config.disable_twin)
		global_status.transformations.emplace_back(new twin_reduction(global_status.n));

	global_status.num_reductions = global_status.transformations.size();

	if (!called_from_fold)
	{
		if (!config.disable_funnel)
			global_status.transformations.emplace_back(new funnel_reduction(global_status.n));
		if (!config.disable_funnel_fold)
			global_status.transformations.emplace_back(new funnel_fold_reduction(global_status.n));
		if (!config.disable_clique_neighborhood_fast)
			global_status.transformations.emplace_back(new clique_neighborhood_reduction_fast(global_status.n));
		if (!config.disable_extended_domination)
			global_status.transformations.emplace_back(new extended_domination_reduction(global_status.n));
		if (!config.disable_unconfined)
			global_status.transformations.emplace_back(new unconfined_csr_reduction(global_status.n));

		if (config.reduction_style == ReductionConfig::Reduction_Style::EARLY_STRUCTION)
		{
			if (!config.disable_decreasing_struction)
				global_status.transformations.emplace_back(make_decreasing_struction(config, global_status.n));
			if (!config.disable_plateau_struction)
				global_status.transformations.push_back(make_plateau_struction(config, global_status.n));

			global_status.num_reductions = global_status.transformations.size();
			if (!config.disable_blow_up)
				global_status.transformations.push_back(make_increasing_struction(config, global_status.n));

			if (!config.disable_critical_set)
				global_status.transformations.emplace_back(new critical_set_reduction(global_status.n));
			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
			// if (!config.disable_bound_reduction)
			// 	global_status.transformations.emplace_back(new bound_reduction(global_status.n));
			// if (!config.disable_high_degree)
			// 	global_status.transformations.emplace_back(new high_degree_reduction(global_status.n));
			if (!config.disable_cut_vertex)
				global_status.transformations.emplace_back(new cut_vertex_reduction(global_status.n));
		}
		else if (config.reduction_style == ReductionConfig::Reduction_Style::EARLY_CS)
		{
			if (!config.disable_critical_set)
				global_status.transformations.emplace_back(new critical_set_reduction(global_status.n));
			if (!config.disable_decreasing_struction)
				global_status.transformations.emplace_back(make_decreasing_struction(config, global_status.n));
			if (!config.disable_plateau_struction)
				global_status.transformations.push_back(make_plateau_struction(config, global_status.n));

			global_status.num_reductions = global_status.transformations.size();
			if (!config.disable_blow_up)
				global_status.transformations.push_back(make_increasing_struction(config, global_status.n));

			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
			// if (!config.disable_bound_reduction)
			// 	global_status.transformations.emplace_back(new bound_reduction(global_status.n));
			// if (!config.disable_high_degree)
			// 	global_status.transformations.emplace_back(new high_degree_reduction(global_status.n));
			if (!config.disable_cut_vertex)
				global_status.transformations.emplace_back(new cut_vertex_reduction(global_status.n));
		}
		else
		{ // FULL and tests
			if (!config.disable_clique_neighborhood_fast)
				global_status.transformations.emplace_back(new clique_neighborhood_reduction_fast(global_status.n));
			if (!config.disable_clique_neighborhood)
				global_status.transformations.emplace_back(new clique_neighborhood_reduction(global_status.n));
			if (!config.disable_decreasing_struction)
				global_status.transformations.emplace_back(make_decreasing_struction(config, global_status.n));
			if (!config.disable_plateau_struction)
				global_status.transformations.push_back(make_plateau_struction(config, global_status.n));
			if (!config.disable_critical_set)
				global_status.transformations.emplace_back(new critical_set_reduction(global_status.n));
			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
			// if (!config.disable_bound_reduction)
			// 	global_status.transformations.emplace_back(new bound_reduction(global_status.n));
			// if (!config.disable_high_degree)
			// 	global_status.transformations.emplace_back(new high_degree_reduction(global_status.n));
			if (!config.disable_cut_vertex)
				global_status.transformations.emplace_back(new cut_vertex_reduction(global_status.n));

			global_status.num_reductions = global_status.transformations.size();
			if (!config.disable_blow_up)
				global_status.transformations.push_back(make_increasing_struction(config, global_status.n));
		}

		if (!config.disable_heuristic_exclude)
			global_status.transformations.emplace_back(new heuristic_exclude_reduction(global_status.n));
		if (!config.disable_heuristic_include)
			global_status.transformations.emplace_back(new heuristic_include_reduction(global_status.n));
	}

	global_transformation_map.resize(REDUCTION_NUM);
	for (size_t i = 0; i < global_status.transformations.size(); i++)
	{
		global_transformation_map[global_status.transformations[i]->get_reduction_type()] = i;
	}

	set_local_reductions = [this, called_from_fold, &config]()
	{
		if (this->config.reduction_style == ReductionConfig::Reduction_Style::DENSE)
		{
			status.transformations = make_reduction_vector<
				neighborhood_reduction, fold2_reduction, clique_reduction,
				domination_reduction, twin_reduction, clique_neighborhood_reduction_fast>(status.n);
		}
		else
		{
			status.transformations.clear();
			if (!config.disable_neighborhood)
				status.transformations.emplace_back(new neighborhood_reduction(global_status.n));
			if (!config.disable_fold1)
				status.transformations.emplace_back(new fold1_reduction(global_status.n));
			if (!config.disable_fold2)
				status.transformations.emplace_back(new fold2_reduction(global_status.n));
			if (!config.disable_clique)
				status.transformations.emplace_back(new clique_reduction(global_status.n));
			if (!config.disable_domination)
				status.transformations.emplace_back(new domination_reduction(global_status.n));
			if (!config.disable_basic_se)
				status.transformations.emplace_back(new single_edge_reduction(global_status.n));
			if (!config.disable_extended_se)
				status.transformations.emplace_back(new extended_single_edge_reduction(global_status.n));
			if (!config.disable_twin)
				status.transformations.emplace_back(new twin_reduction(global_status.n));
			if (!config.disable_extended_domination)
				status.transformations.emplace_back(new extended_domination_reduction(global_status.n));
			if (!config.disable_decreasing_struction)
				status.transformations.emplace_back(make_decreasing_struction(config, global_status.n));
		}
		status.num_reductions = status.transformations.size();

		local_transformation_map.resize(REDUCTION_NUM);
		for (size_t i = 0; i < status.transformations.size(); i++)
		{
			local_transformation_map[status.transformations[i]->get_reduction_type()] = i;
		}
	};

	for (auto &re : status.transformations)
	{
		if (re->get_model_path() != "")
			gnn.change_parameters(re->get_model_path());
	}
}

branch_and_reduce_algorithm::~branch_and_reduce_algorithm() {}

void branch_and_reduce_algorithm::push_nodes(size_t nodes)
{
	resize(status.n + nodes);
	status.push_nodes(nodes);
}
void branch_and_reduce_algorithm::pop_nodes(size_t nodes)
{
	resize(status.n - nodes);
	status.pop_nodes(nodes);
}
void branch_and_reduce_algorithm::resize(size_t size)
{
	set_1.resize(size);
	set_2.resize(size);
	double_set.resize(2 * size);
	zero_vec.resize(size, 0);
	for (auto &transformation : status.transformations)
	{
		transformation->marker.resize(size);
	}
}

size_t branch_and_reduce_algorithm::deg(NodeID node) const
{
	return status.graph[node].size();
}

template <bool update_dependency_checking>
void branch_and_reduce_algorithm::set(NodeID node, IS_status mis_status, bool push_modified)
{
	assert(status.node_status[node] == IS_status::not_set && "Node status set");
	assert(status.remaining_nodes > 0 && "No nodes remaining to set");
	status.node_status[node] = mis_status;

	status.remaining_nodes--;
	status.graph.hide_node(node);

	if (push_modified)
		status.modified_stack.push_back(node);

	if (mis_status == IS_status::included)
	{
		status.is_weight += status.weights[node];

		for (auto neighbor : status.graph[node])
		{
			assert(status.node_status[neighbor] == IS_status::not_set && "neighbor already set");
			status.node_status[neighbor] = IS_status::excluded;
			status.remaining_nodes--;
			status.graph.hide_node(neighbor);

			status.modified_stack.push_back(neighbor);
			// if (update_dependency_checking)
			add_next_level_neighborhood(neighbor);
		}
	}
	else /*if (update_dependency_checking)*/
	{
		add_next_level_neighborhood(node);
	}
}
template void branch_and_reduce_algorithm::set<true>(NodeID node, IS_status mis_status, bool push_modified);
template void branch_and_reduce_algorithm::set<false>(NodeID node, IS_status mis_status, bool push_modified);

void branch_and_reduce_algorithm::unset(NodeID node, bool restore)
{
	if (status.node_status[node] == IS_status::not_set)
		return;
	if (status.node_status[node] == IS_status::included)
	{
		status.is_weight -= status.weights[node];
	}

	status.node_status[node] = IS_status::not_set;
	status.remaining_nodes++;

	if (restore)
		status.graph.restore_node(node);
}

void branch_and_reduce_algorithm::fill_global_greedy()
{
	auto &nodes = buffers[0];
	nodes.resize(global_status.n);
	nodes.clear();

	for (size_t node = 0; node < global_status.n; node++)
	{
		if (global_status.node_status[node] == IS_status::not_set)
			nodes.push_back(node);
	}

	// sort in descending order by node weights
	std::sort(nodes.begin(), nodes.end(), [this](NodeID lhs, NodeID rhs)
			  { return global_status.weights[lhs] > global_status.weights[rhs]; });

	for (NodeID node : nodes)
	{
		bool free_node = true;

		for (NodeID neighbor : global_status.graph[node])
		{
			if (global_status.node_status[neighbor] == IS_status::included)
			{
				free_node = false;
				break;
			}
		}

		if (free_node)
		{
			global_status.node_status[node] = IS_status::included;
			global_status.is_weight += global_status.weights[node];
		}
		else
		{
			global_status.node_status[node] = IS_status::excluded;
		}
	}
}

void branch_and_reduce_algorithm::greedy_initial_is(graph_access &G, std::vector<NodeID> &tmp_buffer)
{
	auto nodes = tmp_buffer;
	nodes.resize(G.number_of_nodes());

	for (size_t i = 0; i < nodes.size(); i++)
	{
		nodes[i] = i;
	}

	// sort in descending order by node weights
	std::sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs)
			  { return G.getNodeWeight(lhs) > G.getNodeWeight(rhs); });

	for (NodeID node : nodes)
	{
		bool free_node = true;

		forall_out_edges(G, edge, node)
		{
			NodeID neighbor = G.getEdgeTarget(edge);
			if (G.getPartitionIndex(neighbor) == 1)
			{
				free_node = false;
				break;
			}
		}
		endfor

			if (free_node) G.setPartitionIndex(node, 1);
	}
}

void branch_and_reduce_algorithm::compute_ils_pruning_bound()
{
	is_ils_best_solution = true;
	best_solution_status = status;

	ch.disable_cout();
	auto config_cpy = config;
	config_cpy.time_limit = config_cpy.time_limit * status.n / total_ils_node_count / 100;
	best_weight = status.reduction_offset + status.is_weight + run_ils(config_cpy, *local_graph, buffers[0], config.max_swaps);
	best_is_time = t.elapsed();
	ch.enable_cout();

	struction_log::instance()->set_best(get_current_is_weight() + best_weight, best_is_time);
	std::cout << get_current_is_weight() + best_weight << ", " << best_is_time << std::endl;
}

// NodeWeight branch_and_reduce_algorithm::compute_partition_pruning_bound(int k)
// {
// 	graph_access G;
// 	std::vector<NodeID> reverse_mapping(status.remaining_nodes,0);
// 	build_graph_access(G, reverse_mapping);
// 	partition_cover ub_solver(k);
// 	ub_solver.create_partition(G, config);
// 	// set original weights to solve
// 	forall_nodes(G, node)
// 	{
// 		G.setNodeWeight(node, status.weights[reverse_mapping[node]]);
// 	} endfor
// 	NodeWeight upper_bound = ub_solver.solve_partition(G, config);
// 	return upper_bound;
// }

NodeWeight branch_and_reduce_algorithm::compute_cover_pruning_bound(int &n_cliques)
{
	// Gather remaining nodes
	auto &nodes = buffers[0];
	nodes.clear();

	for (size_t node = 0; node < status.n; node++)
	{
		if (status.node_status[node] == IS_status::not_set)
			nodes.emplace_back(node);
	}

	// Sort by descending weight
	// Break ties by degree
	std::sort(nodes.begin(), nodes.end(), [&](NodeID lhs, NodeID rhs)
			  { return (status.weights[lhs] > status.weights[rhs] || (status.weights[lhs] == status.weights[rhs] && deg(lhs) > deg(rhs))); });

	// Compute node mapping
	NodeID current_node = 0;
	auto &node_mapping = buffers[1];
	node_mapping.resize(status.n);

	for (NodeID node : nodes)
		node_mapping[node] = current_node++;

	// Init cliques
	auto &clique = buffers[2];
	auto &clique_sizes = buffers[3];
	auto &covered = buffers[4];
	auto &clique_weight = weight_buffer;

	clique.resize(nodes.size());
	clique_weight.resize(nodes.size());
	clique_sizes.resize(nodes.size());
	covered.resize(nodes.size());

	std::fill(clique_weight.begin(), clique_weight.end(), 0);
	std::fill(clique_sizes.begin(), clique_sizes.end(), 0);
	std::fill(covered.begin(), covered.end(), false);

	for (NodeID node : nodes)
		clique[node_mapping[node]] = node_mapping[node];

	auto &neigh_clique_sizes = buffers[5];
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
		clique_weight[heaviest_clique] = std::max(clique_weight[heaviest_clique], status.weights[node]);
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

size_t branch_and_reduce_algorithm::run_ils(const ReductionConfig &config, graph_access &G, std::vector<NodeID> &tmp_buffer, size_t max_swaps)
{
	if (!config.perform_hils)
	{
		greedy_initial_is(G, tmp_buffer);
		ils local_search(config);
		local_search.perform_ils(G, max_swaps);
	}
	else
	{
		if (config.reduce_and_peel)
		{
			ReduceAndPeel reducer(G);
			reducer.reduceAndPeel();
		}
		else
			greedy_initial_is(G, tmp_buffer);
		hils local_search(config);
		local_search.perform_ils(G, max_swaps);
	}

	size_t solution_weight = 0;

	forall_nodes(G, node)
	{
		if (G.getPartitionIndex(node) == 1)
		{
			solution_weight += G.getNodeWeight(node);
		}
	}
	endfor

		return solution_weight;
}

void branch_and_reduce_algorithm::init_transformation_step(reduction_ptr &reduction)
{
	if (config.initial_filter && !reduction->has_run && reduction->has_filtered_marker)
	{
		reduction->marker.current.clear();
		if (reduction->get_model_path() != "" && status.remaining_nodes > 1 && config.gnn_filter)
		{

			timer t;
			gnn.change_parameters(reduction->get_model_path());
			double t_parse = t.elapsed();

			t.restart();
			const float *y = gnn.predict(this);

			reduction->marker.current.clear();
			if (reduction->get_reduction_type() == heuristic_exclude || reduction->get_reduction_type() == heuristic_include)
			{
				for (NodeID u = 0; u < this->status.graph.size(); u++)
				{
					if (status.node_status[u] != IS_status::not_set)
						continue;

					if (y[u] > 0.0f)
						reduction->marker.current.push_back(u);
				}
				printf("%s added %ld vertices from %ld\n", reduction->get_model_path().c_str(), reduction->marker.current.size(), status.remaining_nodes);
			}
			else
			{
				int c = 0;
				for (NodeID u = 0; u < this->status.graph.size(); u++)
				{
					if (status.node_status[u] != IS_status::not_set)
						continue;
					c++;
					if (y[u] > 0.0f)
						reduction->marker.current.push_back(u);
				}

				printf("%s %lf parse, added %ld/%d in %lf seconds\n", reduction->get_model_path().c_str(), t_parse, reduction->marker.current.size(), c, t.elapsed());
			}
		}
		else
		{
			for (NodeID node = 0; node < status.n; node++)
			{
				if (status.node_status[node] != branch_and_reduce_algorithm::IS_status::not_set)
					continue;
				if (reduction->is_suited(node, this))
				{
					reduction->marker.current.push_back(node);
				}
			}
		}

		reduction->has_run = true;
		reduction->has_filtered_marker = false;
		reduction->marker.clear_next();
		return;
	}
	if (!reduction->has_run)
	{
		reduction->marker.fill_current_ascending(status.n);
		reduction->marker.clear_next();
		reduction->has_run = true;
	}
	else
	{
		reduction->marker.get_next();
	}
}

void branch_and_reduce_algorithm::add_next_level_node(NodeID node)
{
	for (size_t i = 0; i < status.transformations.size(); i++)
	{
		auto &reduction = status.transformations[i];
		if (reduction->has_run)
			reduction->marker.add(node);
	}
}

void branch_and_reduce_algorithm::add_next_level_neighborhood(NodeID node)
{
	// node has been excluded in mis -> neighboring vertices are interesting for next round of reduction
	for (auto neighbor : status.graph[node])
	{
		add_next_level_node(neighbor);
	}
}

void branch_and_reduce_algorithm::add_next_level_neighborhood(const std::vector<NodeID> &nodes)
{
	for (auto node : nodes)
	{
		add_next_level_neighborhood(node);
	}
}

void branch_and_reduce_algorithm::reduce_graph()
{
	initial_reduce();
	status = std::move(global_status);
}

void branch_and_reduce_algorithm::initial_reduce()
{
	std::swap(global_transformation_map, local_transformation_map);
	status = std::move(global_status);

	bool further_impovement = status.remaining_nodes > 0;
	min_kernel = status.remaining_nodes;

	while (further_impovement && status.remaining_nodes > 0 && t.elapsed() <= config.time_limit)
	{
		reduce_graph_internal_before_blow_up();

		further_impovement = false;

		if (!config.disable_blow_up && status.transformations.size() != status.num_reductions && !heuristically_reducing)
			cyclic_blow_up();

		if (status.remaining_nodes == 0)
		{
			if (config.print_reduction_info)
				print_reduction_info();
			break;
		}

		size_t oldn = status.remaining_nodes;
		reduce_graph_internal_after_blow_up();
		further_impovement = oldn != status.remaining_nodes;
		if (!further_impovement && heuristically_reducing)
		{ // testing all again after heuristic reductions finished
			heuristically_reducing = false;
			further_impovement = true;
		}
		if (config.print_reduction_info)
			print_reduction_info();
	}
	status.modified_stack.push_back(BRANCHING_TOKEN);

	global_status = std::move(status);
	std::swap(global_transformation_map, local_transformation_map);
}

void branch_and_reduce_algorithm::reduce_graph_internal_before_blow_up()
{
	size_t active_reduction_index = 0;
	while (active_reduction_index < status.transformations.size() && t.elapsed() <= config.time_limit)
	{
		if (status.transformations[active_reduction_index]->get_reduction_type() == struction_blow)
			break;

		auto &reduction = status.transformations[active_reduction_index];

		init_transformation_step(reduction);
		bool progress = reduction->reduce(this);
		/* if (progress && config.print_reduction_info) */
		/* { */
		/* 	print_reduction_progress(); */
		/* } */
		active_reduction_index = progress ? 0 : active_reduction_index + 1;
		if (status.remaining_nodes == 0)
			break;
	}

	timeout |= t.elapsed() > config.time_limit;
}

void branch_and_reduce_algorithm::reduce_graph_internal_after_blow_up()
{
	// get struction_blow index
	auto it = std::find_if(status.transformations.begin(), status.transformations.end(), [&](auto &reduction)
						   { return reduction->get_reduction_type() == struction_blow; });
	if (it == status.transformations.end())
		return;

	size_t active_reduction_index = std::distance(status.transformations.begin(), it) + 1;

	while (active_reduction_index < status.transformations.size() && t.elapsed() <= config.time_limit)
	{
		auto &reduction = status.transformations[active_reduction_index];

		init_transformation_step(reduction);
		bool progress = reduction->reduce(this);
		if (progress)
			return;
		active_reduction_index++;
		if (status.remaining_nodes == 0)
			break;
	}

	timeout |= t.elapsed() > config.time_limit;
	// std::cout << "\ncurrent weight: " << status.is_weight << "  weight offset: " << status.reduction_offset << "  remaining nodes: " << status.remaining_nodes << std::endl;
}

bool branch_and_reduce_algorithm::blow_up_graph_internal()
{
	bool blown_up = false;
	size_t active_blow_up_index = status.num_reductions;
	auto &blow_up = status.transformations[active_blow_up_index];
	init_transformation_step(blow_up);
	blown_up |= blow_up->reduce(this);
	timeout |= t.elapsed() > config.time_limit;
	return blown_up;
}

void branch_and_reduce_algorithm::cyclic_blow_up()
{
	blowing_up = true;
	int phase_count = 0;

#ifdef OUTPUT_GRAPH_CONVERGENCE
	std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif
	const double alpha = config.global_blow_up_factor;
	const int X = config.max_unimproving_phases;
	while (status.remaining_nodes < alpha * min_kernel && phase_count < X)
	{
		if (timeout || !blow_up_graph_internal())
			break;
		// REDUCE
		reduce_graph_internal_before_blow_up();

		if (status.remaining_nodes < min_kernel)
		{
			min_kernel = status.remaining_nodes;
#ifdef OUTPUT_GRAPH_CONVERGENCE
			std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif
			phase_count = 0;
		}
		else if (config.backtrack_style == ReductionConfig::IMMEDIATE_EXCLUDE ||
				 config.backtrack_style == ReductionConfig::IMMEDIATE_TIE_BREAKING)
		{
			++phase_count;
			undo_blow_up();
		}
	}
	if (config.backtrack_style != ReductionConfig::NO_BACKTRACK)
		undo_blow_up();
#ifdef OUTPUT_GRAPH_CONVERGENCE
	std::cout << t.elapsed() << "," << min_kernel << std::endl;
#endif

	blowing_up = false;
	// std::cout << "\ncurrent weight: " << status.is_weight << "  weight offset: " << status.reduction_offset << "  remaining nodes: " << status.remaining_nodes << std::endl;
}

bool branch_and_reduce_algorithm::branch_reduce_recursive()
{
	recursive_mapping.resize(status.n, 0);
	recursive_comp_map.resize(status.n, 0);
	build_graph_access(recursive_graph, recursive_mapping);
	size_t comp_count = strongly_connected_components().strong_components(recursive_graph, recursive_comp_map);

	if (comp_count == 1)
		return false;

	graph_extractor extractor;
	const auto time_limit = config.time_limit;

	status.modified_stack.pop_back();

	for (size_t node = 0; node < status.remaining_nodes; node++)
	{
		recursive_graph.setPartitionIndex(node, recursive_comp_map[node] + 2);
	}

	for (size_t i = 0; i < comp_count; i++)
	{
		if (t.elapsed() > time_limit)
		{
			timeout = true;
			break;
		}

		recursive_local_mapping.clear();
		graph_access G;
		extractor.extract_block(recursive_graph, G, i + 2, recursive_local_mapping);

		config.time_limit = time_limit - t.elapsed();
		config.disable_heuristic_exclude = true;
		config.disable_heuristic_include = true;

		branch_and_reduce_algorithm br_alg(G, config, true);
		ch.disable_cout();
		timeout = !br_alg.run_branch_reduce();
		ch.enable_cout();

		if (timeout)
			break;

		br_alg.apply_branch_reduce_solution(G);

		forall_nodes(G, node)
		{
			recursive_graph.setPartitionIndex(recursive_local_mapping[node], G.getPartitionIndex(node));
		}
		endfor
	}

	config.time_limit = time_limit;

	if (timeout)
	{
		return false;
	}

	forall_nodes(recursive_graph, node)
	{
		if (recursive_graph.getPartitionIndex(node) == 1)
			set<false>(recursive_mapping[node], IS_status::included);
	}
	endfor

		return true;
}

void branch_and_reduce_algorithm::branch_reduce_single_component()
{
	best_weight = 0;

	if (status.n == 0)
	{
		return;
	}
	if (status.n == 1)
	{
		set<false>(0, IS_status::included);
		std::cout << (get_current_is_weight() + status.is_weight + status.reduction_offset) << "," << t.elapsed() << std::endl;
		return;
	}

	resize(status.n);

	std::vector<NodeID> node_order(status.n);
	std::iota(node_order.begin(), node_order.end(), 0);
	std::sort(node_order.begin(), node_order.end(), [&](const NodeID lhs, const NodeID rhs)
			  { return deg(lhs) > deg(rhs) || (deg(lhs) == deg(rhs) && status.weights[lhs] > status.weights[rhs]); });

	if (status.n > ILS_SIZE_LIMIT)
		compute_ils_pruning_bound();

	if (best_weight > weight_bound)
	{
		// std::cerr << "pruning for solving subgraphs in reducion" << std::endl;
		return;
	}

	size_t i = 0;
	while (i < status.n)
	{
		if (t.elapsed() > config.time_limit && best_weight != 0)
		{
			timeout = true;
			break;
		}

		NodeID branch_node = node_order[i];
		if (i == status.n - 1)
		{
			if (status.node_status[branch_node] == IS_status::not_set)
				set<false>(branch_node, IS_status::included);

			update_best_solution();
			reverse_branching();
			i = status.branching_stack.back().pos;
			continue;
		}
		bool improvement_possible = false;
		if (status.remaining_nodes > 0)
		{
			int n_cliques = 0;
			NodeWeight better_bound = compute_cover_pruning_bound(n_cliques) + status.reduction_offset + status.is_weight;

			improvement_possible = better_bound > best_weight;
		}
		if (!improvement_possible)
		{
			// no improvement_possible --> backtrack
			if (status.branching_stack.size() > 0 && status.branching_stack.back().node == branch_node)
			{
				status.branching_stack.pop_back();
				unset(branch_node);

				if (i == 0)
					break;
			}

			update_best_solution();
			reverse_branching();
			i = status.branching_stack.back().pos;
			continue;
		}

		// improvement_possible
		if (status.node_status[branch_node] == IS_status::not_set)
		{
			// first time reaching this node in this branch --> include it to IS
			set(branch_node, IS_status::included, false);
			status.branching_stack.emplace_back(branch_node, i);
		}
		else if (status.branching_stack.size() > 0 && status.branching_stack.back().node == branch_node)
		{
			if (status.node_status[branch_node] == IS_status::included)
			{
				// second time reaching this node in this branch --> exclude it from IS; was included before
				status.node_status[branch_node] = IS_status::excluded;
				status.is_weight -= status.weights[branch_node];
				add_next_level_neighborhood(branch_node);
			}
			else
			{
				// third time reaching this node in this branch --> backtack
				status.branching_stack.pop_back();
				unset(branch_node);

				if (i == 0)
					break;

				update_best_solution();
				reverse_branching();
				i = status.branching_stack.back().pos;
				continue;
			}
		}
		else
		{
			i++;
			continue;
		}

		size_t old_n = status.n;
		reduce_graph_internal_before_blow_up();

		assert(status.remaining_nodes <= status.n && "Graph size and remaining nodes mismatch");
		status.modified_stack.push_back(BRANCHING_TOKEN);
		if (old_n < status.n)
		{
			node_order.resize(status.n);
			std::iota(node_order.begin() + old_n, node_order.begin() + status.n, old_n);
			std::sort(node_order.begin() + old_n, node_order.begin() + status.n, [&](const NodeID lhs, const NodeID rhs)
					  { return deg(lhs) > deg(rhs) || (deg(lhs) == deg(rhs) && status.weights[lhs] > status.weights[rhs]); });

			resize(status.n);
		}

		if (status.remaining_nodes > SPLIT_CC_LIMIT && branch_reduce_recursive())
		{
			status.modified_stack.push_back(BRANCHING_TOKEN);
			update_best_solution();
			reverse_branching();
			i = status.branching_stack.back().pos;
		}
		else
		{
			i++;
		}
	}

	restore_best_local_solution();
}

graph_access &branch_and_reduce_algorithm::kernelize()
{
	initial_reduce();
	build_global_graph_access();
	return global_graph;
}

bool branch_and_reduce_algorithm::run_branch_reduce(NodeWeight bound)
{
	weight_bound = bound;
	return run_branch_reduce();
}
bool branch_and_reduce_algorithm::run_branch_reduce()
{
	t.restart();
	if (!config.disable_bound_reduction)
	{
		buffers.resize(7);
		buffers[4].reserve(global_status.n);
		buffers[5].reserve(global_status.n);
		buffers[6].reserve(global_status.n);
	}
	initial_reduce();
	min_kernel = status.remaining_nodes;
	kernelization_time = t.elapsed();
	kernelization_offset = global_status.is_weight + global_status.reduction_offset;
	struction_log::instance()->set_best(kernelization_offset, kernelization_time);

	if (config.console_log)
	{
		std::cout << "%reduction_nodes " << global_status.remaining_nodes << "\n";
		std::cout << "%reduction_offset " << kernelization_offset << "\n";
		std::cout << "%reduction_time " << kernelization_time << "\n";
	}

	if (global_status.remaining_nodes == 0)
	{
		if (config.console_log)
			std::cout << "solution weight and time: " << get_current_is_weight() << "," << t.elapsed() << std::endl;
		max_min_kernel_comp = 0;

		restore_best_global_solution();
		if (config.console_log)
			std::cout << "solution weight and time: " << get_current_is_weight() << "," << t.elapsed() << std::endl;
		update_best_global_solution();
		return true;
	}

	build_global_graph_access();

	std::vector<int> comp_map(global_status.remaining_nodes, 0);
	size_t comp_count = strongly_connected_components().strong_components(global_graph, comp_map);

	// std::cout << "%components " << comp_count << "\n";

	for (size_t node = 0; node < global_status.remaining_nodes; node++)
	{
		global_graph.setPartitionIndex(node, comp_map[node]);
	}

	std::vector<size_t> comp_size(comp_count, 0);
	for (auto comp_id : comp_map)
	{
		comp_size[comp_id]++;
	}

	total_ils_node_count = 0;
	std::vector<size_t> comp_idx(comp_count, 0);
	for (size_t i = 0; i < comp_count; i++)
	{
		comp_idx[i] = i;
		if (comp_size[i] > ILS_SIZE_LIMIT)
		{
			total_ils_node_count += comp_size[i];
		}
	}

	std::sort(comp_idx.begin(), comp_idx.end(), [&comp_size](const size_t lhs, const size_t rhs)
			  { return comp_size[lhs] < comp_size[rhs]; });

	// std::cout << "%max_component " << comp_size[comp_idx.back()] << "\n";
	max_min_kernel_comp = comp_size[comp_idx.back()];

	graph_extractor extractor;

	buffers.resize(7);
	buffers[4].reserve(global_status.n);
	buffers[5].reserve(global_status.n);
	buffers[6].reserve(global_status.n);

	for (size_t i : comp_idx)
	{
		if (t.elapsed() > config.time_limit)
		{
			timeout = true;
			break;
		}

		// std::cout << "%connected component " << i << ":  " << comp_size[i] << std::endl;

		local_mapping.clear();
		graph_access G;
		extractor.extract_block(global_graph, G, i, local_mapping);
		local_graph = &G;

		status = graph_status(*local_graph);

		set_local_reductions();

		branch_reduce_single_component();

		if (best_weight > weight_bound)
			return false;

		for (size_t node = 0; node < local_mapping.size(); node++)
		{
			global_status.node_status[global_mapping[local_mapping[node]]] = status.node_status[node];
		}

		global_status.is_weight += best_weight;
		// best_weight = global_status.is_weight;
		assert(status.n <= global_status.remaining_nodes && "Graph size and remaining nodes mismatch");
		global_status.remaining_nodes -= status.n;
		local_graph = nullptr;
	}

	if (timeout)
	{
		fill_global_greedy();
	}
	update_best_global_solution();
	restore_best_global_solution();
	return !timeout;
}

void branch_and_reduce_algorithm::update_best_global_solution()
{
	NodeWeight current_weight = global_status.is_weight + global_status.reduction_offset;
	if (current_weight > best_weight)
	{
		struction_log::instance()->set_best(current_weight, t.elapsed());
		best_solution_status = global_status;
		best_weight = current_weight;
		best_is_time = t.elapsed();
		is_ils_best_solution = false;

		std::cout << (current_weight) << "," << t.elapsed() << std::endl;
	}
}
void branch_and_reduce_algorithm::update_best_solution()
{
	NodeWeight current_weight = status.is_weight + status.reduction_offset;
	if (current_weight > best_weight)
	{
		struction_log::instance()->set_best(get_current_is_weight() + current_weight, t.elapsed());
		best_solution_status = status;
		best_weight = current_weight;
		best_is_time = t.elapsed();
		is_ils_best_solution = false;

		std::cout << (get_current_is_weight() + current_weight) << "," << t.elapsed() << std::endl;
	}
}

void branch_and_reduce_algorithm::undo_blow_up()
{
	while (status.remaining_nodes > min_kernel)
	{
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();
		if (node == MODIFIED_TOKEN)
		{ // for reductions that do not reduce but only modify the graph
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->restore(this);
			continue;
		}
		if (status.node_status[node] == IS_status::folded)
		{
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->restore(this);
		}
		else
		{
			unset(node);
		}
	}
}

void branch_and_reduce_algorithm::reverse_branching()
{
	// discard topmost branching token
	if (status.modified_stack.size() > 0)
	{
		status.modified_stack.pop_back();
	}
	else
	{
		return;
	}

	while (status.modified_stack.size() > 0 && !is_token(status.modified_stack.back()))
	{
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();

		if (node == MODIFIED_TOKEN)
		{ // for reductions that do not reduce but only modify the graph
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->restore(this);
			continue;
		}
		else if (status.node_status[node] == IS_status::folded)
		{
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->restore(this);
			continue;
		}

		unset(node);
	}
}

void branch_and_reduce_algorithm::restore_best_local_solution()
{
	status = best_solution_status;
	if (is_ils_best_solution)
	{
		assert(local_graph->number_of_nodes() == status.n && "graph size and status.n mismatch");
		for (size_t node = 0; node < status.n; node++)
		{
			if (local_graph->getPartitionIndex(node) == 1)
			{
				status.node_status[node] = IS_status::included;
			}
			else
			{
				status.node_status[node] = IS_status::excluded;
			}
		}
		return;
	}

	status.modified_stack.pop_back();

	while (status.modified_stack.size() > 0)
	{
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();

		if (node == BRANCHING_TOKEN)
		{
			status.graph.restore_node(status.branching_stack.back().node);
			status.branching_stack.pop_back();
			continue;
		}

		if (node == MODIFIED_TOKEN)
		{ // for reductions that do not reduce but only modify the graph
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->apply(this);
		}
		else if (status.node_status[node] == IS_status::folded)
		{
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[local_transformation_map[type]]->apply(this);
		}
		else
		{
			status.graph.restore_node(node);
		}
	}
}

void branch_and_reduce_algorithm::restore_best_global_solution()
{
	status = std::move(global_status);
	status.modified_stack.pop_back();

	while (status.modified_stack.size() > 0)
	{
		NodeID node = status.modified_stack.back();
		status.modified_stack.pop_back();
		if (node == MODIFIED_TOKEN)
		{ // for reductions that do not reduce but only modify the graph
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[global_transformation_map[type]]->apply(this);
		}
		else if (status.node_status[node] == IS_status::folded)
		{
			auto type = status.folded_stack.back();
			status.folded_stack.pop_back();
			status.transformations[global_transformation_map[type]]->apply(this);
		}
		else
		{
			status.graph.restore_node(node);
		}
	}
}

NodeWeight branch_and_reduce_algorithm::get_current_is_weight() const
{
	return global_status.is_weight + global_status.reduction_offset;
}

NodeWeight branch_and_reduce_algorithm::get_is_weight() const
{
	return status.is_weight + status.reduction_offset;
}

NodeID branch_and_reduce_algorithm::get_heuristically_reduced_vertices() const
{
	return global_status.heuristically_reduced_n;
}

void branch_and_reduce_algorithm::build_global_graph_access()
{
	global_mapping.resize(global_status.remaining_nodes, 0);
	std::swap(status, global_status);
	build_graph_access(global_graph, global_mapping);
	std::swap(status, global_status);
}

void branch_and_reduce_algorithm::build_graph_access(graph_access &G, std::vector<NodeID> &reverse_mapping) const
{
	std::vector<NodeID> mapping(status.graph.size(), UINT_MAX);
	size_t edge_count = 0;

	// Get number of edges and reorder nodes
	size_t node_counter = 0;
	for (NodeID node = 0; node < status.graph.size(); ++node)
	{
		if (status.node_status[node] == IS_status::not_set)
		{
			edge_count += status.graph[node].size();

			mapping[node] = node_counter;
			reverse_mapping[node_counter] = node;
			node_counter++;
		}
	}

	// Create the adjacency array
	std::vector<EdgeID> xadj(status.remaining_nodes + 1);
	std::vector<NodeID> adjncy(edge_count + 1);
	size_t adjncy_counter = 0;

	for (size_t i = 0; i < status.remaining_nodes; ++i)
	{
		xadj[i] = adjncy_counter;

		for (auto neighbor : status.graph[reverse_mapping[i]])
		{
			if (mapping[neighbor] == i)
				continue;
			if (mapping[neighbor] == UINT_MAX)
				continue;
			adjncy[adjncy_counter++] = mapping[neighbor];
		}

		std::sort(std::begin(adjncy) + xadj[i], std::begin(adjncy) + adjncy_counter);
	}
	xadj[status.remaining_nodes] = adjncy_counter;

	// Build the graph
	G.build_from_metis(status.remaining_nodes, &xadj[0], &adjncy[0]);

	forall_nodes(G, node)
	{
		G.setNodeWeight(node, status.weights[reverse_mapping[node]]);
	}
	endfor
}

void branch_and_reduce_algorithm::build_induced_neighborhood_subgraph(graph_access &G, NodeID source_node)
{
	buffers[0].clear();
	set_1.clear();

	for (NodeID neighbor : status.graph[source_node])
	{
		buffers[0].push_back(neighbor);
		set_1.add(neighbor);
	}

	build_induced_subgraph(G, buffers[0], set_1, buffers[1]);
}

void branch_and_reduce_algorithm::build_induced_subgraph(graph_access &G, const std::vector<NodeID> &nodes, const fast_set &nodes_set, std::vector<NodeID> &reverse_mapping)
{
	assert(nodes.size() > 1 && "Induced subgraph must have at least two nodes");
	size_t edge_count = 0;

	for (size_t i = 0; i < nodes.size(); i++)
	{
		NodeID node = nodes[i];
		for (auto neighbor : status.graph[node])
		{
			if (nodes_set.get(neighbor))
				edge_count++;
		}

		reverse_mapping[node] = i;
	}

	G.start_construction(nodes.size(), edge_count);

	for (NodeID node : nodes)
	{
		NodeID new_node = G.new_node();
		G.setNodeWeight(new_node, status.weights[node]);

		for (auto neighbor : status.graph[node])
		{
			if (nodes_set.get(neighbor))
			{
				G.new_edge(new_node, reverse_mapping[neighbor]);
			}
		}
	}

	G.finish_construction();
}

void branch_and_reduce_algorithm::reverse_reduction(graph_access &G, graph_access &reduced_G, std::vector<NodeID> &reverse_mapping)
{
	// build_global_graph_access();

	forall_nodes(reduced_G, node)
	{
		if (reduced_G.getPartitionIndex(node) == 1)
		{
			status.node_status[reverse_mapping[node]] = IS_status::included;
		}
		else
		{
			status.node_status[reverse_mapping[node]] = IS_status::excluded;
		}
	}
	endfor

	global_status = std::move(status);
	restore_best_global_solution();
	apply_branch_reduce_solution(G);
}

void branch_and_reduce_algorithm::apply_branch_reduce_solution(graph_access &G)
{
	forall_nodes(G, node)
	{
		if (status.node_status[node] == IS_status::included)
		{
			G.setPartitionIndex(node, 1);
		}
		else
		{
			G.setPartitionIndex(node, 0);
		}
	}
	endfor
}


/****************************************************
 * added for generating training_data
 * ************************************************/

void branch_and_reduce_algorithm::generate_initial_reduce_data(graph_access &G, std::vector<std::vector<int>> &reduction_data)
{
	std::swap(global_transformation_map, local_transformation_map);
	status = std::move(global_status);
	std::vector<reduction_type> global = {critical_set, struction_plateau, struction_blow, unconfined_csr, cut_vertex};
	std::vector<NodeID> label;
	label.reserve(status.n);

	for (size_t i = 0; i < status.transformations.size(); i++)
	{
		auto t = status.transformations[i]->get_reduction_type();
		if (std::find(global.begin(), global.end(), t) != global.end())
		{
			auto &&reduction = status.transformations[i];
			if (status.transformations[i]->get_reduction_type() == unconfined_csr)
				reduction->marker.fill_current_ascending(status.n);
			label.clear();
			int l = (int)reduction->generate_global_data(this, label);
			for (NodeID node : label)
				reduction_data[i][node] = l;
		}
	}

	for (NodeID node = 0; node < status.n; node++)
	{
		for (size_t i = 0; i < status.transformations.size(); i++)
		{
			auto t = status.transformations[i]->get_reduction_type();
			if (std::find(global.begin(), global.end(), t) != global.end())
				continue;
			auto &&reduction = status.transformations[i];
			label.clear();
			int l = reduction->generate_data(this, node, label);
			if (label.size() == 0)
				label.push_back(node);
			for (NodeID node : label)
				reduction_data[i][node] = l;
			assert(status.remaining_nodes == status.n && "graph should be original size");
		}
	}

	global_status = std::move(status);
	std::swap(global_transformation_map, local_transformation_map);
}

// function to get transformations
void branch_and_reduce_algorithm::get_transformation_names(std::vector<std::string> &names)
{
	names.clear();
	for (size_t i = 0; i < global_status.transformations.size(); i++)
	{
		names.push_back(global_status.transformations[i]->get_reduction_name());
	}
}

void branch_and_reduce_algorithm::print_reduction_progress()
{
	// print progress bar for remaining nodes in graph
	size_t i = 0;
	double factor = 1.0 / status.n;
	std::cout << "\r [";
	for (; i < (status.n - status.remaining_nodes) * factor * 50; i++)
		std::cout << "=";
	for (; i < status.n * factor * 50; i++)
		std::cout << " ";
	std::cout << "] " << (status.n - status.remaining_nodes) * factor * 50 << "% \n";
	std::cout << "\033[A" << std::flush;
}

void branch_and_reduce_algorithm::print_reduction_info()
{
	std::cout << "\n\n\t\t Reduction Info" << std::endl;
	std::cout << "====================================================" << std::endl;
	for (size_t i = 0; i < status.transformations.size(); i++)
	{
		auto &reduction = status.transformations[i];
		std::cout << i << " ";
		reduction->print_reduction_type();
		std::cout << " reduced " << reduction->reduced_nodes << "\tnodes in " << reduction->reduction_time << " s" << std::endl;
	}
}
