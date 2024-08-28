/******************************************************************************
* branch_and_reduce_algorithm.h
*
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

#pragma once

// system includes
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <iostream>
#include <sstream>

// local includes
#include "fast_set.h"
#include "config.h"
#include "timer.h"
#include "cout_handler.h"
#include "definitions.h"
#include "data_structure/graph_access.h"
#include "dynamic_graph.h"
#include "mwis_finder.h"
#include "key_functions.h"
#include "LRConv.h"

// reduction rules
#include "general_reduction.h"
#include "neighborhood_reduction.h"
#include "fold1_reduction.h"
#include "fold2_reduction.h"
#include "clique_reduction.h"
#include "twin_reduction.h"
#include "domination_reduction.h"
#include "extended_domination_reduction.h"
#include "funnel_reduction.h"
#include "funnel_fold_reduction.h"
#include "single_edge_reduction.h"
#include "extended_single_edge_reduction.h"
#include "clique_neighborhood_reduction.h"
#include "clique_neighborhood_fast_reduction.h"
#include "generalized_fold_reduction.h"
#include "critical_set_reduction.h"
#include "cut_vertex_reduction.h"
#include "heavy_set3_reduction.h"
#include "heavy_set_reduction.h"
// #include "high_degree_reduction.h"
// #include "bound_reduction.h"
#include "struction_reductions.h"
#include "heuristic_exclude_reduction.h"
#include "heuristic_include_reduction.h"
#include "reductions.h"


class branch_and_reduce_algorithm {
public:
	enum IS_status { not_set, included, excluded, folded };

	size_t min_kernel;
	size_t max_min_kernel_comp = 0;

    double kernelization_time;
    NodeWeight kernelization_offset = 0;
    double best_is_time;
    bool timeout = false;
	cout_handler ch;

	//generate training data
	// std::vector<bool> is_included_vertex;
	// std::vector<bool> is_excluded_vertex;
    void get_transformation_names(std::vector<std::string> &names);

private:
    friend general_reduction;
    friend neighborhood_reduction;
	friend fold1_reduction;
	friend fold2_reduction;
    friend clique_neighborhood_reduction;
    friend clique_neighborhood_reduction_fast;
    friend critical_set_reduction;
    friend heavy_set_reduction;
    friend heavy_set3_reduction;
    friend clique_reduction;
	friend funnel_reduction;
	friend funnel_fold_reduction;
	friend single_edge_reduction;
	friend extended_single_edge_reduction;
    friend domination_reduction;
	friend extended_domination_reduction;	
    friend twin_reduction;
	// friend high_degree_reduction;
	// friend bound_reduction;
    friend generalized_fold_reduction;
    friend cut_vertex_reduction;
	friend heuristic_include_reduction;
	friend heuristic_exclude_reduction;

    // friend path_reduction;
    template<typename struction_type, reduction_type type, int new_nodes>
    friend class iterative_struction;
    template<typename pick_function, typename struction_type, reduction_type type>
    friend class blow_up_struction;

    friend mwis_finder;
    friend extended_struction<false>;
    friend extended_struction<true>;
    friend original_struction<false>;
    friend original_struction<true>;

    friend ApproximateIncreaseKey;
    friend DegreeKey;
    friend IncreaseKey;

	friend LRConv;

    struct node_pos {
		NodeID node;
		size_t pos;

		node_pos(NodeID node = 0, size_t pos = 0) : node(node), pos(pos) {}
	};

	struct graph_status {
		size_t n = 0;
		size_t m = 0;
		size_t remaining_nodes = 0;
		NodeWeight is_weight = 0;
		NodeWeight reduction_offset = 0;
		dynamic_graph graph;
		std::vector<NodeWeight> weights;
		std::vector<IS_status> node_status;
		std::vector<std::vector<bool>> reduction_node_status;

        std::vector<reduction_ptr> transformations; //reductions + blow_ups.
        size_t num_reductions;
        // int blow_up_index;
		NodeID heuristically_reduced_n = 0;

		std::vector<reduction_type> folded_stack;
		std::vector<node_pos> branching_stack;
		std::vector<NodeID> modified_stack;

		graph_status() = default;

		graph_status(graph_access& G) :
                n(G.number_of_nodes()), m(G.number_of_edges()), remaining_nodes(n), graph(G), weights(n, 0), node_status(n, IS_status::not_set)
		{
			modified_stack.reserve(2*n + 1);
			branching_stack.reserve(n);
			folded_stack.reserve(2*n);

			forall_nodes(G, node) {
				weights[node] = G.getNodeWeight(node);
			} endfor
		}

		void resize(size_t size) {
		    weights.resize(size, 0);
		    node_status.resize(size, IS_status::not_set);
			if (reduction_node_status.size() > 0 )
				reduction_node_status.resize(size, std::vector<bool>(num_reductions, true));
            n = size;
		}

		void push_nodes(size_t nodes) {
            remaining_nodes += nodes;
            graph.push_nodes(nodes);
            resize(n + nodes);
		}

		void pop_nodes(size_t nodes) {
            for (NodeID n = this->n - nodes; n < this->n; ++n)
                remaining_nodes -= node_status[n] == IS_status ::not_set;
            graph.pop_nodes(nodes);
            resize(n - nodes);
		}
	};

    static constexpr NodeID BRANCHING_TOKEN = std::numeric_limits<NodeID>::max();
    static constexpr NodeID MODIFIED_TOKEN = BRANCHING_TOKEN - 1; 

	static bool is_token(NodeID node) {
		return node == BRANCHING_TOKEN;
	}

	// lower graph size limit for when to use ils pruning
	static constexpr size_t ILS_SIZE_LIMIT = 50;

	// min number of remaining nodes to split up connected components
	static constexpr size_t SPLIT_CC_LIMIT = 100;

	// max number of neighbors to recurse on induced neighborhood graph in reductions
	static constexpr size_t REDU_RECURSION_LIMIT = 150;

	// max number of nodes to use all reductions during branch reduce
	static constexpr size_t FULL_REDUCTIONS_RECURSION_LIMIT = 50;


    graph_access global_graph;
	ReductionConfig config;
	graph_status best_solution_status;
	NodeWeight best_weight = 0;
	NodeWeight weight_bound = std::numeric_limits<NodeWeight>::max();
	timer t;
	timer reduction_timer;
	bool is_ils_best_solution = false;

	bool blowing_up = false;
	bool heuristically_reducing = false;

	graph_status global_status;
	std::vector<NodeID> global_mapping;
	std::vector<size_t> global_transformation_map;
	std::vector<reduction_type> global_transformations;
	std::vector<reduction_type> expensive_transformations;

	size_t total_ils_node_count;

	graph_status status;
    graph_access* local_graph;
	std::vector<NodeID> local_mapping;
	std::vector<size_t> local_transformation_map;
	std::vector<reduction_ptr> local_reductions;

    graph_access recursive_graph;
	std::vector<NodeID> recursive_mapping;
	std::vector<int> recursive_comp_map;
	std::vector<NodeID> recursive_local_mapping;

	std::function<void()> set_local_reductions;

	// presized objects for temporary use
	fast_set set_1;
	fast_set set_2;
	fast_set double_set;
    std::vector<std::vector<NodeID>> buffers;
    std::vector<NodeWeight> weight_buffer;
    std::vector<bool> bool_buffer;
    std::vector<NodeID> zero_vec;

	LRConv gnn;

    void resize(size_t size);

    template<bool update_dependency_checking = true>
	void set(NodeID node, IS_status status, bool push_modified = true);
	void unset(NodeID node, bool restore = true);

	void fill_global_greedy();
	void compute_ils_pruning_bound();
	NodeWeight compute_cover_pruning_bound(int& n_cliques);
    // NodeWeight compute_partition_pruning_bound(int k);

	void init_transformation_step(reduction_ptr &reduction);
	void init_global_transformation_step(reduction_ptr & reduction);
	void add_next_level_node(NodeID node);
	void add_next_level_neighborhood(NodeID node);
	void add_next_level_neighborhood(const std::vector<NodeID>& nodes);

    // void initial_filter_vertices_for_reduction();

	void reduce_graph_internal_before_blow_up();
    void reduce_graph_internal_after_blow_up();
    void reduce_graph_by_vertex_internal(bool full);
    bool blow_up_graph_internal();
    void cyclic_blow_up();
	bool branch_reduce_recursive();
	void branch_reduce_single_component();
    void initial_reduce();


    void update_best_solution();
    void update_best_global_solution();
    void undo_blow_up();
    void reverse_branching();
	void restore_best_local_solution();
	void restore_best_global_solution();

    void apply_local_reduction(size_t type);
    void restore_local_reduction(size_t type);
    void restore_reduction(size_t type);
    void apply_reduction(size_t type);

    void build_global_graph_access();
	void build_induced_neighborhood_subgraph(graph_access& G, NodeID source_node);
	void build_induced_subgraph(graph_access& G, const std::vector<NodeID>& nodes, const fast_set& nodes_set, std::vector<NodeID>& reverse_mapping);

    void push_nodes(size_t nodes);
    void pop_nodes(size_t nodes);
public:
	branch_and_reduce_algorithm(graph_access& G, const ReductionConfig& config, bool called_from_fold = false);
	~branch_and_reduce_algorithm();


    graph_access &kernelize();
    size_t deg(NodeID node) const;
	void reduce_graph();
    bool run_branch_reduce(NodeWeight weight_bound);
    bool run_branch_reduce();

    static size_t run_ils(const ReductionConfig& config, graph_access& G, std::vector<NodeID>& tmp_buffer, size_t max_swaps);
	static void greedy_initial_is(graph_access& G, std::vector<NodeID>& tmp_buffer);

	NodeWeight get_current_is_weight() const;
	NodeWeight get_is_weight() const;
    NodeID get_heuristically_reduced_vertices() const;
    void reverse_reduction(graph_access & G, graph_access & reduced_G, std::vector<NodeID> & reverse_mapping);
	void apply_branch_reduce_solution(graph_access & G);

	void build_graph_access(graph_access & G, std::vector<NodeID>& reverse_mapping) const;

	// added for mmwis
 	void force_into_independent_set(std::vector<NodeID> const &nodes);
    void exclude_nodes(std::vector<NodeID> const &nodes);
    void update_independent_set(std::vector<bool> & independent_set);
    void set_node_status(std::vector<bool> & independent_set , graph_access & G, graph_access & reduced, std::vector<NodeID> & reverse_mapping);
    // added for training data
    // void get_training_data_for_graph_size(graph_access &graph, NodeID n, std::vector<std::vector<bool>> &reduction_data, size_t i);
    void get_training_data_for_graph_size(graph_access &graph, NodeID n, std::vector<std::vector<bool>> &reduction_data, std::vector<bool> &include_data, std::vector<bool> &exclude_data, size_t i);
    void get_exact_training_data_for_graph_size(graph_access &graph, NodeID n);
    void pick_nodes_by_BFS(NodeID n, std::vector<NodeID> &nodes_vec, fast_set &nodes_set);
	void pick_nodes_by_BFS_sample(NodeID n, std::vector<NodeID> &nodes_vec, fast_set &nodes_set);
    void pick_nodes_by_nodeID(NodeID n, std::vector<NodeID> &nodes_vec, fast_set &nodes_set);
    void generate_initial_reduce_data(std::vector<std::vector<bool>> &reduction_data, size_t i);
	void generate_initial_reduce_data(graph_access &G, std::vector<std::vector<bool>> &reduction_data); 

    // printing
    void print_reduction_info();
    void print_dyn_graph(){
		print_dyn_graph(global_status);
	}
	void print_dyn_graph(graph_status & s)
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

	void print_graph(graph_access&  G)
	{
		std::cout << "__________________________________" << std::endl;
		std::cout << "Graph:" << std::endl;
		forall_nodes(G, node)
		{
			std::cout << "Node: " << node << " w(" << G.getNodeWeight(node) << ")\t index(" << G.getPartitionIndex(node) << ")\t :";
			forall_out_edges(G, e, node)
			{
				std::cout << G.getEdgeTarget(e) << " ";
			} endfor	
			std::cout << std::endl;
		} endfor
	}

	// template<typename T>
	// void print_subgraph(graph_status& G, T & nodes)
	void print_subgraph(graph_status& G, std::vector<NodeID>& nodes)
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
			std::cout << "\n"<< std::endl;
		}
	}
    void print_reduction_progress();
};
