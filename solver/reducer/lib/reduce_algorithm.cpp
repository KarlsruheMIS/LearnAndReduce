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

#include "reduce_algorithm.h"

#include "definitions.h"
#include "reductions.h"
#include "solution_check.h"
#include "struction_log.h"
#include "LRConv.h"
#include "struction_reductions.h"

#include "tiny_solver.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <climits>

// node limit for subgraph solver
#define MAX_NODES 128

constexpr NodeID reduce_algorithm::BRANCHING_TOKEN;
constexpr NodeID reduce_algorithm::MODIFIED_TOKEN;

reduce_algorithm::reduce_algorithm(graph_access &G, const ReductionConfig &config, bool called_from_fold)
	: config(config), global_status(G), set_1(global_status.n), set_2(global_status.n), double_set(global_status.n * 2),
	  buffers(4), zero_vec(global_status.n, 0), gnn(called_from_fold ? 0 : G.number_of_nodes())
{
	subgraph_solver = tiny_solver_init(G.number_of_nodes());
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
	if (!config.disable_extended_twin)
		global_status.transformations.emplace_back(new extended_twin_reduction(global_status.n));

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
			if (!config.disable_extended_domination_reverse)
				global_status.transformations.emplace_back(new extended_domination_reverse_reduction(global_status.n));
			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
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

			if (!config.disable_extended_domination_reverse)
				global_status.transformations.emplace_back(new extended_domination_reverse_reduction(global_status.n));
			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
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
			if (!config.disable_extended_domination_reverse)
				global_status.transformations.emplace_back(new extended_domination_reverse_reduction(global_status.n));
			if (!config.disable_generalized_fold)
				global_status.transformations.emplace_back(new generalized_fold_reduction(global_status.n));
			if (!config.disable_heavy_set)
				global_status.transformations.emplace_back(new heavy_set_reduction(global_status.n));
			if (!config.disable_heavy_set3)
				global_status.transformations.emplace_back(new heavy_set3_reduction(global_status.n));
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

reduce_algorithm::~reduce_algorithm()
{
	tiny_solver_free(subgraph_solver);
}

void reduce_algorithm::push_nodes(size_t nodes)
{
	resize(status.n + nodes);
	status.push_nodes(nodes);
}
void reduce_algorithm::pop_nodes(size_t nodes)
{
	resize(status.n - nodes);
	status.pop_nodes(nodes);
}
void reduce_algorithm::resize(size_t size)
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

size_t reduce_algorithm::deg(NodeID node) const
{
	return status.graph[node].size();
}

template <bool update_dependency_checking>
void reduce_algorithm::set(NodeID node, IS_status mis_status, bool push_modified)
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
template void reduce_algorithm::set<true>(NodeID node, IS_status mis_status, bool push_modified);
template void reduce_algorithm::set<false>(NodeID node, IS_status mis_status, bool push_modified);

void reduce_algorithm::unset(NodeID node, bool restore)
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

void reduce_algorithm::init_transformation_step(reduction_ptr &reduction)
{
	if (!reduction->has_run && reduction->has_filtered_marker)
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
				if (config.verbose)
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

				if (config.verbose)
					printf("%s %lf parse, added %ld/%d in %lf seconds\n", reduction->get_model_path().c_str(), t_parse, reduction->marker.current.size(), c, t.elapsed());
			}
		}
		else
		{
			for (NodeID node = 0; node < status.n; node++)
			{
				if (status.node_status[node] != reduce_algorithm::IS_status::not_set)
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

void reduce_algorithm::add_next_level_node(NodeID node)
{
	for (size_t i = 0; i < status.transformations.size(); i++)
	{
		auto &reduction = status.transformations[i];
		if (reduction->has_run)
			reduction->marker.add(node);
	}
}

void reduce_algorithm::add_next_level_neighborhood(NodeID node)
{
	// node has been excluded in mis -> neighboring vertices are interesting for next round of reduction
	for (auto neighbor : status.graph[node])
	{
		add_next_level_node(neighbor);
	}
}

void reduce_algorithm::add_next_level_neighborhood(const std::vector<NodeID> &nodes)
{
	for (auto node : nodes)
	{
		add_next_level_neighborhood(node);
	}
}

void reduce_algorithm::reduce_graph()
{
	initial_reduce();
	status = std::move(global_status);
}

void reduce_algorithm::initial_reduce()
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

void reduce_algorithm::reduce_graph_internal_before_blow_up()
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

void reduce_algorithm::reduce_graph_internal_after_blow_up()
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

bool reduce_algorithm::blow_up_graph_internal()
{
	bool blown_up = false;
	size_t active_blow_up_index = status.num_reductions;
	auto &blow_up = status.transformations[active_blow_up_index];
	init_transformation_step(blow_up);
	blown_up |= blow_up->reduce(this);
	timeout |= t.elapsed() > config.time_limit;
	return blown_up;
}

void reduce_algorithm::cyclic_blow_up()
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

graph_access &reduce_algorithm::kernelize()
{
	initial_reduce();
	build_global_graph_access();
	return global_graph;
}

void reduce_algorithm::undo_blow_up()
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

void reduce_algorithm::restore_best_local_solution()
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
		assert(status.modified_stack.back() != BRANCHING_TOKEN);
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

void reduce_algorithm::restore_best_global_solution()
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

NodeWeight reduce_algorithm::get_current_is_weight() const
{
	return global_status.is_weight + global_status.reduction_offset;
}

NodeWeight reduce_algorithm::get_is_weight() const
{
	return status.is_weight + status.reduction_offset;
}

NodeID reduce_algorithm::get_heuristically_reduced_vertices() const
{
	return global_status.heuristically_reduced_n;
}

csr_graph reduce_algorithm::build_global_graph_csr()
{
	auto &reverse_map = buffers[0];
	auto &map = buffers[1];
	auto &V = buffers[2];
	auto &E = buffers[3];
	auto &W = weight_buffer;

	reverse_map.resize(status.remaining_nodes, status.n);
	map.resize(status.n, status.n);

	NodeID N = status.remaining_nodes;
	EdgeID M = 0;
	NodeID u = 0;
	for (NodeID v = 0; v < status.n; v++)
	{
		if (!status.node_status[v] == IS_status::not_set)
			continue;
		map[v] = u;
		reverse_map[u] = v;
		u++;
		M += status.graph[v].size();
	}
	V.resize(N + 1);
	E.resize(M);
	W.resize(N);

	// build CSR
	std::vector<NodeID> tmp;
	tmp.reserve(N);
	EdgeID ei = 0;
	for (NodeID u : reverse_map)
	{
		assert(map[u] < status.n && "map is not set");
		V[map[u]] = ei;
		W[map[u]] = status.weights[u];

		// get the sorted, unreduced neighborhood
		tmp.clear();
		for (NodeID neighbor : status.graph[u])
		{
			if (status.node_status[neighbor] == IS_status::not_set)
				tmp.push_back(map[neighbor]);
		}

		std::sort(tmp.begin(), tmp.end());
		for (NodeID v : tmp)
			E[ei++] = v;
	}
	V[N] = ei;
	return (csr_graph){N, V.data(), E.data(), W.data()};
}

void reduce_algorithm::tiny_solver_solve_neighbourhood(const NodeWeight Wl, const NodeID v)
{
	tiny_solver_clear(subgraph_solver);
	auto &neighbor_set = set_1;
	set_1.clear();
	if (deg(v) > MAX_NODES)
	{
		subgraph_solver->node_limit_exceeded = 1;
		return;
	}
	// set mappings for the solver
	NodeID n = 0;
	for (NodeID u : status.graph[v])
	{
		assert(status.node_status[u] == IS_status::not_set && "node already set");
		neighbor_set.add(u);
		subgraph_solver->forward_map[u] = n;
		subgraph_solver->reverse_map[n] = u;
		subgraph_solver->subgraph_W[n] = status.weights[u];
		n++;
	}
	subgraph_solver->subgraph_N = n;

	// set the adjacency matrix
	for (NodeID u : status.graph[v])
	{
		NodeID u_subgraph = subgraph_solver->forward_map[u];
		subgraph_solver->subgraph[u_subgraph][u_subgraph] = 1;
		for (NodeID v : status.graph[u])
		{
			if (u <= v)
				continue;
			NodeID v_subgraph = subgraph_solver->forward_map[v];
			if (!neighbor_set.get(v))
				continue;
			subgraph_solver->subgraph[u_subgraph][v_subgraph] = 1;
			subgraph_solver->subgraph[v_subgraph][u_subgraph] = 1;
		}
	}

	tiny_solver_solve(subgraph_solver, config.subgraph_time_limit, Wl);
}

void reduce_algorithm::tiny_solver_solve_subgraph(const NodeWeight Wl, const std::vector<NodeID> &nodes, const fast_set &nodes_set)
{
	if (nodes.size() == 0)
		return;
	tiny_solver_clear(subgraph_solver);
	assert(std::all_of(nodes.begin(), nodes.end(), [&](NodeID v)
					   { return status.node_status[v] == IS_status::not_set; }) &&
		   "all nodes need to be unset");

	if (nodes.size() > MAX_NODES)
	{
		subgraph_solver->node_limit_exceeded = 1;
		return;
	}
	NodeID n = 0;
	// set mappings for the solver
	for (NodeID u : nodes)
	{
		assert(nodes_set.get(u) && "node not in set");
		subgraph_solver->forward_map[u] = n;
		subgraph_solver->reverse_map[n] = u;
		subgraph_solver->subgraph_W[n] = status.weights[u];
		n++;
	}
	subgraph_solver->subgraph_N = n;

	// set the adjacency matrix
	for (NodeID u : nodes)
	{
		NodeID u_subgraph = subgraph_solver->forward_map[u];
		subgraph_solver->subgraph[u_subgraph][u_subgraph] = 1;
		for (NodeID v : status.graph[u])
		{
			if (u <= v)
				continue;
			NodeID v_subgraph = subgraph_solver->forward_map[v];
			if (!nodes_set.get(v))
				continue;
			subgraph_solver->subgraph[u_subgraph][v_subgraph] = 1;
			subgraph_solver->subgraph[v_subgraph][u_subgraph] = 1;
		}
	}
	tiny_solver_solve(subgraph_solver, config.subgraph_time_limit, Wl);
}
void reduce_algorithm::build_global_graph_access()
{
	global_mapping.resize(global_status.remaining_nodes, 0);
	std::swap(status, global_status);
	build_graph_access(global_graph, global_mapping);
	std::swap(status, global_status);
}
void reduce_algorithm::build_graph_access(graph_access &G, std::vector<NodeID> &reverse_mapping)
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

void reduce_algorithm::reverse_reduction(graph_access &G, graph_access &reduced_G, std::vector<NodeID> &reverse_mapping)
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

void reduce_algorithm::apply_branch_reduce_solution(graph_access &G)
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

void reduce_algorithm::generate_initial_reduce_data(graph_access &G, std::vector<std::vector<int>> &reduction_data)
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
			while (status.modified_stack.size() > 0)
			{
				NodeID node = status.modified_stack.back();
				status.modified_stack.pop_back();
				unset(node);
			}
			assert(status.remaining_nodes == status.n && "graph should be original size");
		}
	}

	global_status = std::move(status);
	std::swap(global_transformation_map, local_transformation_map);
}

void reduce_algorithm::get_transformation_names(std::vector<std::string> &names)
{
	names.clear();
	for (size_t i = 0; i < global_status.transformations.size(); i++)
	{
		names.push_back(global_status.transformations[i]->get_reduction_name());
	}
}

/****************************************************
 * added for printing
 * ************************************************/

void reduce_algorithm::print_reduction_progress()
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

void reduce_algorithm::print_reduction_info()
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
void reduce_algorithm::print_dyn_graph()
{
	print_dyn_graph(global_status);
}
void reduce_algorithm::print_dyn_graph(graph_status &s)
{
	std::cout << "__________________________________" << std::endl;
	std::cout << "Dynamic graph:" << std::endl;
	for (size_t node = 0; node < s.graph.size(); node++)
	{
		switch (s.node_status[node])
		{
		case IS_status::included:
			std::cout << " \t \t\t included w(" << s.weights[node] << ")\t";
			break;
		case IS_status::excluded:
			std::cout << " \t \t\t excluded w(" << s.weights[node] << ")\t";
			break;
		case IS_status::folded:
			std::cout << " \t \t\t folded w(" << s.weights[node] << ")\t";
			break;
		case IS_status::not_set:
			std::cout << "not_set\t\t w(" << s.weights[node] << ")\t";
			break;
		}
		s.graph.print_neighbors(node);
	}
	std::cout << "__________________________________" << std::endl;
}
void reduce_algorithm::print_graph(graph_access &G)
{
	std::cout << "__________________________________" << std::endl;
	std::cout << "Graph:" << std::endl;
	forall_nodes(G, node)
	{
		std::cout << "Node: " << node << " w(" << G.getNodeWeight(node) << ")\t index(" << G.getPartitionIndex(node) << ")\t :";
		forall_out_edges(G, e, node)
		{
			std::cout << G.getEdgeTarget(e) << " ";
		}
		endfor
				std::cout
			<< std::endl;
	}
	endfor
}
void reduce_algorithm::print_subgraph(graph_status &G, std::vector<NodeID> &nodes)
{
	std::cout << "__________________________________" << std::endl;
	std::cout << "Subgraph:" << std::endl;
	for (auto node : nodes)
	{
		std::cout << "Node: " << node << " w(" << status.weights[node] << ")\t status(" << status.node_status[node] << ")\t :";
		for (auto neighbor : status.graph[node])
		{
			std::cout << neighbor << " ";
		}
		std::cout << "\n"
				  << std::endl;
	}
}