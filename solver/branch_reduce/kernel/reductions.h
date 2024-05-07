/******************************************************************************
 * reductions.h
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

#ifndef _REDUCTIONS_H
#define _REDUCTIONS_H

// local includes
#include "definitions.h"
#include "fast_set.h"
#include "sized_vector.h"
#include "dynamic_graph.h"
#include "extended_struction.h"
#include "original_struction.h"
#include "vertex_marker.h"
#include "key_functions.h"
#include "maxNodeHeap.h"
#include "reduction_config.h"

// system includes
#include <vector>
#include <memory>
#include <array>
#include <data_structure/priority_queues/MaxHeap.h>

class branch_and_reduce_algorithm;
using Struction_Type = ReductionConfig::Struction_Type;
using Key_Type = ReductionConfig::Key_Type;

// initial reductions need to be at the beginning
enum reduction_type
{
    fold1,
    neighborhood,
    fold2,
    clique,
    clique_neighborhood_fast,
    clique_neighborhood,
    domination,
    twin,
    critical_set,
    generalized_fold,
    single_edge,
    extended_single_edge,
    heavy_vertex,
    heavy_set,
    heavy_set3,
    cut_vertex,
    component,
    funnel,
    funnel_fold,
    path,
    struction_decrease,
    struction_plateau,
    struction_blow
};
constexpr size_t REDUCTION_NUM = 23;

struct general_reduction
{
    general_reduction(size_t n) : marker(n) {}
    virtual ~general_reduction() {}
    virtual general_reduction *clone() const = 0;

    virtual reduction_type get_reduction_type() const = 0;
    virtual void print_reduction_type() {
        std::cout << get_reduction_name() << ": \t" ;
    };
    virtual std::string get_reduction_name() = 0;
    virtual std::string get_model_path() { return ""; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) = 0;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) {}
    virtual void restore(branch_and_reduce_algorithm *br_alg) {}
    virtual void apply(branch_and_reduce_algorithm *br_alg) {}
    virtual void reset(branch_and_reduce_algorithm *br_alg, size_t comp_size) {}
    virtual bool is_suited(NodeID v, branch_and_reduce_algorithm *br_alg) { return true; }

    bool has_run = false;
    bool has_filtered_marker = false;
    NodeID reduced_nodes = 0;
    double reduction_time = 0.0;
    vertex_marker marker;

	template<typename F>
	void for_each_changed_vertex(branch_and_reduce_algorithm* br_alg, F f);
    inline NodeWeight get_neighborhood_weight(NodeID v, branch_and_reduce_algorithm* br_alg);
    inline NodeID get_max_weight_neighbor(NodeID v, branch_and_reduce_algorithm* br_alg);
    inline void get_neighborhood_set(NodeID v, branch_and_reduce_algorithm* br_alg, fast_set& neighborhood_set);
    inline void get_neighborhood_vector(NodeID v, branch_and_reduce_algorithm* br_alg, sized_vector<NodeID>& neighborhood_vec);
    inline bool try_neighborhood_reduction(NodeID v, branch_and_reduce_algorithm* br_alg, NodeWeight neighborhood_weight);
    inline bool solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution, graph_access &graph, branch_and_reduce_algorithm *br_alg, sized_vector<NodeID> &nodes_vec, const fast_set &nodes_set, sized_vector<NodeID> &reverse_mapping, bool apply_solution=false);
    inline bool solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, graph_access &neighborhood_graph, branch_and_reduce_algorithm *br_alg, NodeID v, bool apply_solution=false);
    inline bool solve_graph(NodeWeight &solution, graph_access &graph, ReductionConfig &config, NodeWeight weight_bound, bool apply_solution=false);
    inline bool is_reduced(NodeID v, branch_and_reduce_algorithm *br_alg);
};

// simple reductions:
struct neighborhood_reduction : public general_reduction
{
    neighborhood_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = true;
    }
    ~neighborhood_reduction() {}
    virtual neighborhood_reduction *clone() const final { return new neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::neighborhood; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual std::string get_reduction_name() final { return "neighborhood"; }

    bool is_suited(NodeID v, branch_and_reduce_algorithm *br_alg);
};

struct fold1_reduction : public general_reduction
{
    fold1_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = true;
    }
    ~fold1_reduction() {}
    virtual fold1_reduction *clone() const final { return new fold1_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::fold1; }
    virtual std::string get_reduction_name() final { return "fold1"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;
    bool is_suited(NodeID v, branch_and_reduce_algorithm *br_alg);

private:
    struct fold_nodes
    {
        NodeID deg1_node;
        NodeID fold_node;
    };

    struct restore_data
    {
        fold_nodes nodes;
        NodeWeight deg1_weight;
        // dynamic_graph::neighbor_list old_neighbors;
    };

    void fold(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);

    std::vector<restore_data> restore_vec;
};

struct fold2_reduction : public general_reduction
{
    fold2_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~fold2_reduction() {}
    virtual fold2_reduction *clone() const final { return new fold2_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::fold2; }
    virtual std::string get_reduction_name() final { return "fold2"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;
    bool is_suited(NodeID v, branch_and_reduce_algorithm *br_alg);

private:
    struct fold_nodes
    {
        NodeID deg2_node;
        std::vector<NodeID> neighbors;
    };
    enum fold_case
    {
        triangle_mid,
        triangle_min,
        v_shape_max,
        v_shape_mid,
        v_shape_min
    };

    struct restore_data
    {
        fold_nodes nodes;
        NodeWeight deg2_weight;
        fold_case fold_case;
        std::array<std::vector<NodeID>, 2> node_vecs;
    };

    void fold_triangle_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_triangle_min_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_max_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_min_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);

    int v_shape_min_count = 0;

    std::vector<restore_data> restore_vec;
};

// clique reductions:
struct clique_neighborhood_reduction_fast : public general_reduction
{
    clique_neighborhood_reduction_fast(size_t n) : general_reduction(n) {}
    ~clique_neighborhood_reduction_fast() {}
    virtual clique_neighborhood_reduction_fast *clone() const final { return new clique_neighborhood_reduction_fast(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood_fast; }
    virtual std::string get_reduction_name() final { return "clique_nbh_fast"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

struct clique_neighborhood_reduction : public general_reduction
{
    clique_neighborhood_reduction(size_t n) : general_reduction(n) {}
    ~clique_neighborhood_reduction() {}
    virtual clique_neighborhood_reduction *clone() const final { return new clique_neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood; }
    virtual std::string get_reduction_name() final { return "clique_nbh"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

    bool partition_into_cliques(NodeID v, branch_and_reduce_algorithm *br_alg);
    bool expand_clique(NodeID max_neighbor, sized_vector<NodeID> &neighbors_vec, fast_set &clique_neighbors_set, branch_and_reduce_algorithm *br_alg);

    NodeWeight target_weight;
    NodeWeight neighbor_weights;
};

struct clique_reduction : public general_reduction
{
    clique_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~clique_reduction() {}
    virtual clique_reduction *clone() const final { return new clique_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique; }
    virtual std::string get_reduction_name() final { return "clique"; }
    virtual std::string get_model_path() final { return "models/clique.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    struct weighted_node
    {
        NodeID node;
        NodeWeight weight;
    };

    struct restore_data
    {
        weighted_node isolated;
        std::vector<NodeID> non_isolated;

        restore_data() = default;
        restore_data(const weighted_node &isolated, std::vector<NodeID> &&non_isolated) : isolated(isolated), non_isolated(std::move(non_isolated)){};
    };

    void fold(branch_and_reduce_algorithm *br_alg, const weighted_node &isolated, std::vector<NodeID> &&non_isolated);

    std::vector<restore_data> restore_vec;
};

struct funnel_reduction : public general_reduction
{
    funnel_reduction(size_t n) : general_reduction(n) {}
    ~funnel_reduction() {}
    virtual funnel_reduction *clone() const final { return new funnel_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::funnel; }
    virtual std::string get_reduction_name() final { return "funnel"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    struct fold_data
    {
        NodeID node;
        NodeID funnel_neighbor;
    };

    struct restore_data
    {
        NodeID node;
        NodeID funnel_neighbor;
        sized_vector<NodeID> remaining_neighbors;
        std::vector<std::vector<NodeID>> node_vecs;
    };

    bool is_funnel(NodeID v, NodeID &funnel_neighbor, branch_and_reduce_algorithm *br_alg, fast_set &funnel_set, sized_vector<NodeID> &funnel_nodes);
    bool is_clique(branch_and_reduce_algorithm *br_alg, fast_set &clique_set, sized_vector<NodeID> &clique_nodes);
    void fold(const fold_data &data, fast_set &funnel_set, branch_and_reduce_algorithm *br_alg);

    std::vector<restore_data> restore_vec;
};
struct funnel_fold_reduction : public general_reduction
{
    funnel_fold_reduction(size_t n) : general_reduction(n) {}
    ~funnel_fold_reduction() {}
    virtual funnel_fold_reduction *clone() const final { return new funnel_fold_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::funnel_fold; }
    virtual std::string get_reduction_name() final { return "funnel fold"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    struct fold_data
    {
        NodeID node;
        NodeID funnel_neighbor;
    };

    struct restore_data
    {
        NodeID node;
        NodeID funnel_neighbor;
        sized_vector<NodeID> remaining_neighbors;
        std::vector<std::vector<NodeID>> node_vecs;
    };

    bool is_funnel(NodeID v, NodeID &funnel_neighbor, branch_and_reduce_algorithm *br_alg, fast_set &funnel_set, sized_vector<NodeID> &funnel_nodes);
    bool is_clique(branch_and_reduce_algorithm *br_alg, fast_set &clique_set, sized_vector<NodeID> &clique_nodes);
    void fold(const fold_data &data, fast_set &funnel_set, branch_and_reduce_algorithm *br_alg);

    std::vector<restore_data> restore_vec;
};

// larger reductions:
struct critical_set_reduction : public general_reduction
{
    critical_set_reduction(size_t n) : general_reduction(n) {}
    ~critical_set_reduction() {}
    virtual critical_set_reduction *clone() const final { return new critical_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::critical_set; }
    virtual std::string get_reduction_name() final { return "critical_set"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
};

struct component_reduction : public general_reduction
{
    component_reduction(size_t n) : general_reduction(n) {}
    ~component_reduction() {}
    virtual component_reduction *clone() const final { return new component_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::component; }
    virtual std::string get_reduction_name() final { return "component"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

private:
    bool check_components(branch_and_reduce_algorithm *br_alg, NodeID u, NodeID &cut_vertex, sized_vector<NodeID> &smallComponent);
    bool build_small_component(NodeID u, branch_and_reduce_algorithm *br_alg, sized_vector<NodeID> &component, std::vector<bool> &component_visited);
    void dfs_fill_visited(NodeID u, branch_and_reduce_algorithm *br_alg, std::vector<bool> &component_visited);
};

struct cut_vertex_reduction : public general_reduction
{
    cut_vertex_reduction(size_t n) : general_reduction(n) {}
    ~cut_vertex_reduction() {}
    virtual cut_vertex_reduction *clone() const final { return new cut_vertex_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::cut_vertex; }
    virtual std::string get_reduction_name() final { return "cut_vertex"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

    bool find_cut_vertex(branch_and_reduce_algorithm *br_alg, NodeID &cut_v, sized_vector<NodeID> &cut_component, std::vector<NodeID> &reverse_mapping, fast_set &tested);
    bool DFS(branch_and_reduce_algorithm *br_alg, NodeID u, int &step, NodeID &cut_vertex, sized_vector<NodeID> &smallComponent);
    bool get_fold_data(branch_and_reduce_algorithm *br_alg, NodeID cut_v, sized_vector<NodeID> &cut_v_included_i, sized_vector<NodeID> &cut_v_included_e, sized_vector<NodeID> &cut_v_excluded_i, sized_vector<NodeID> &cut_v_excluded_e, NodeWeight &large_cutMWIS_weight, NodeWeight &small_cutMWIS_weight);

private:
    struct fold_data
    {
        NodeID cut_vertex;
        NodeWeight cut_vertex_weight;
        NodeWeight large_cutMWIS_weight;
        NodeWeight small_cutMWIS_weight;
        sized_vector<NodeID> cut_component;
    };

    struct restore_data
    {
        fold_data data;
        sized_vector<NodeID> case_cut_v_included_nodes_to_include;
        sized_vector<NodeID> case_cut_v_included_nodes_to_exclude;
        sized_vector<NodeID> case_cut_v_excluded_nodes_to_include;
        sized_vector<NodeID> case_cut_v_excluded_nodes_to_exclude;
    };

    bool check_components(branch_and_reduce_algorithm *br_alg, NodeID u, NodeID &cut_vertex, sized_vector<NodeID> &smallComponent);
    bool build_small_component(NodeID u, branch_and_reduce_algorithm *br_alg, sized_vector<NodeID> &component, std::vector<bool> &component_visited);
    void dfs_fill_visited(NodeID u, branch_and_reduce_algorithm *br_alg, std::vector<bool> &component_visited);
    void fold(branch_and_reduce_algorithm *br_alg, fold_data &data, sized_vector<NodeID> &cut_v_included_i, sized_vector<NodeID> &cut_v_included_e, sized_vector<NodeID> &cut_v_excluded_i, sized_vector<NodeID> &cut_v_excluded_e);

    std::vector<restore_data> restore_vec;
};

struct single_edge_reduction : public general_reduction
{
    single_edge_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~single_edge_reduction() {}
    virtual single_edge_reduction *clone() const final { return new single_edge_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::single_edge; }
    virtual std::string get_reduction_name() final { return "single_edge"; }
    virtual std::string get_model_path() final { return "models/single_edge.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

struct extended_single_edge_reduction : public general_reduction
{
    extended_single_edge_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~extended_single_edge_reduction() {}
    virtual extended_single_edge_reduction *clone() const final { return new extended_single_edge_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::extended_single_edge; }
    virtual std::string get_reduction_name() final { return "extended_single_edge"; }
    virtual std::string get_model_path() final { return "models/extended_single_edge.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

struct twin_reduction : public general_reduction
{
    twin_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~twin_reduction() {}
    virtual twin_reduction *clone() const final { return new twin_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::twin; }
    virtual std::string get_reduction_name() final { return "twin"; }
    virtual std::string get_model_path() final { return "models/twin.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    struct restore_data
    {
        NodeID main;
        NodeID twin;
    };

    void fold(branch_and_reduce_algorithm *br_alg, NodeID main, NodeID twin);

    std::vector<restore_data> restore_vec;
};

struct domination_reduction : public general_reduction
{
    domination_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~domination_reduction() {}
    virtual domination_reduction *clone() const final { return new domination_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::domination; }
    virtual std::string get_reduction_name() final { return "domination"; }
    virtual std::string get_model_path() final { return "models/domination.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

struct generalized_fold_reduction : public general_reduction
{
    generalized_fold_reduction(size_t n) : general_reduction(n) {}
    ~generalized_fold_reduction() {}
    virtual generalized_fold_reduction *clone() const final { return new generalized_fold_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_fold; }
    virtual std::string get_reduction_name() final { return "generalized_fold"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    struct fold_nodes
    {
        NodeID main;
        std::vector<NodeID> MWIS;
    };

    struct restore_data
    {
        fold_nodes nodes;
        NodeWeight main_weight;
        NodeWeight MWIS_weight;
        dynamic_graph::neighbor_list main_neighbor_list;
        std::vector<std::vector<NodeID>> MWIS_node_vecs;
    };

    void fold(branch_and_reduce_algorithm *br_alg, NodeID main_node, fast_set &MWIS_set, NodeWeight MWIS_weight);

    std::vector<restore_data> restore_vec;
};

struct heavy_vertex_reduction : public general_reduction
{
    heavy_vertex_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heavy_vertex_reduction() {}
    virtual heavy_vertex_reduction *clone() const final { return new heavy_vertex_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_vertex; }
    virtual std::string get_reduction_name() final { return "heavy_vertex"; }
    virtual std::string get_model_path() final { return "models/heavy_vertex.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

struct heavy_set_reduction : public general_reduction
{
    heavy_set_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heavy_set_reduction() {}
    virtual heavy_set_reduction *clone() const final { return new heavy_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set; }
    virtual std::string get_reduction_name() final { return "heavy_set"; }
    virtual std::string get_model_path() final { return "models/heavy_set.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

private:
    enum v_combination
    {
        oo,
        uo,
        ov,
        uv
    }; // uo = u fixed as included and v is excluded (start with no increasing to all vertices)
    bool is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, branch_and_reduce_algorithm *br_alg);
    void unset_weights(graph_access &graph, sized_vector<NodeID> &nodes, sized_vector<NodeID> &reverse_mapping);
    void set_weights(graph_access &graph, sized_vector<NodeID> &nodes, sized_vector<NodeID> &reverse_mapping, std::vector<NodeWeight> &weights);
};

struct heavy_set3_reduction : public general_reduction
{
    heavy_set3_reduction(size_t n) : general_reduction(n) {}
    ~heavy_set3_reduction() {}
    virtual heavy_set3_reduction *clone() const final { return new heavy_set3_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set3; }
    virtual std::string get_reduction_name() final { return "heavy_set3"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

private:
    enum v_combination
    {
        ooo,
        uoo,
        ovo,
        oow,
        uvo,
        uow,
        ovw,
        uvw
    }; // uoo = u fixed as included and v, w (=o) are excluded (start with no increasing to all vertices)

    bool is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, fast_set &u_neighbors_set, NodeID w, branch_and_reduce_algorithm *br_alg);
    bool check_u_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_v_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_w_combination(std::vector<NodeWeight> &MWIS_weights);
    void unset_weights(graph_access &graph, sized_vector<NodeID> &nodes, sized_vector<NodeID> &reverse_mapping);
    void set_weights(graph_access &graph, sized_vector<NodeID> &nodes, sized_vector<NodeID> &reverse_mapping, std::vector<NodeWeight> &weights);
};

struct path_reduction : public general_reduction
{
    path_reduction(size_t n) : general_reduction(n) {}
    ~path_reduction() {}
    virtual path_reduction *clone() const final { return new path_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::path; }
    virtual std::string get_reduction_name() final { return "path"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

private:
    branch_and_reduce_algorithm *br_alg;
    NodeWeight w_i_i, w_i_e, w_e_i, w_e_e;                   // c_{f}_{l} ^= weight of maximum weight IS on path (in/ex)cluding first node (f) and (in/ex)cluding last node (l)
    std::vector<std::pair<NodeID, NodeWeight>> node_weights; // tmp vector to store changed node weights
    struct restore_data
    {
        std::vector<std::pair<NodeID, NodeWeight>> node_weights;
        sized_vector<NodeID> path;
        size_t start;
        size_t end;
        size_t offset;
        bool relink = false;
    };

    std::vector<restore_data> restore_vec;

    template <bool add_global = true>
    void enqueue_node(NodeID n);
    bool dequeue_node(sized_vector<NodeID> &queue, NodeID &n, size_t degree);
    bool reduce_degree_one_node();
    bool reduce_path();

    void fold_path(size_t start, size_t end, restore_data &data);
    void reassign_weight(NodeID n, NodeWeight w, restore_data &data);
    void reconnect(NodeID v, NodeID w, restore_data &data);
    void add_reduction_offset(size_t reduction_offset, restore_data &data);
    bool are_connected(NodeID v, NodeID w) const;

    void find_max_deg_1_path(NodeID n, sized_vector<NodeID> &path);
    void find_max_path(NodeID n, sized_vector<NodeID> &path);
    bool next_node_on_path(NodeID current, NodeID last, NodeID first, NodeID &next);

    void find_MIS_on_path(sized_vector<NodeID> &path);

    template <bool track_choices = false>
    void find_MIS_on_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path);
    template <bool track_choices = false>
    void find_MIS_on_deg_1_path(NodeWeight &w_i, NodeWeight &w_e, sized_vector<NodeID> &path);

    void apply(sized_vector<NodeID> &path);
};

template <typename struction_type, reduction_type type, int vertex_increase>
struct iterative_struction : public general_reduction
{
    iterative_struction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~iterative_struction() {}
    virtual iterative_struction *clone() const final { return new iterative_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual std::string get_reduction_name() final { return "struction"; }
    virtual std::string get_model_path() final { 
        if (type == struction_decrease) 
            return "models/decreasing_struction.gnn";
        else
            return ""; 
    }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final { s.restore(br_alg); };
    virtual void apply(branch_and_reduce_algorithm *br_alg) final { s.apply(br_alg); };

private:
    struction_type s;
};

using reduction_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_decrease, -1>;
using plateu_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_plateau, 0>;

template <typename key_function = ApproximateIncreaseKey, typename struction_type = extended_struction<false>, reduction_type type = reduction_type ::struction_blow /**/>
struct blow_up_struction : public general_reduction
{
    blow_up_struction(const ReductionConfig &config, size_t n) : general_reduction(n), update_set(n), distribution(0, 1),
                                                                 generator(config.seed), f(config) {}
    ~blow_up_struction() {}
    virtual blow_up_struction *clone() const final { return new blow_up_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual std::string get_reduction_name() final { return "blow_up_struction"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    // virtual bool reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final { s.apply(br_alg); };
    virtual void reset(branch_and_reduce_algorithm *br_alg, size_t comp_size) final;

private:
    size_t phase_start_kernel;

    fast_set update_set;
    std::vector<std::pair<NodeID, Gain>> update_list;
    std::vector<NodeID> added_list;

    key_function f;
    struction_type s;
    bool restored = false;
    size_t blow_ups;
    branch_and_reduce_algorithm *br_alg;

    MaxHeap<float> queue;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution;

    void init_blow_up_phase();

    bool clean_up_queue();
    bool is_done();

    Gain denoise(float key)
    {
        return static_cast<Gain>(key);
    }

    float apply_noise(Gain key)
    {
        return key + distribution(generator);
    }

    void update_queue_by_key(NodeID n, Gain key)
    {
        if (queue.contains(n))
        {
            if (blow_ups != 0 && update_set.add(n))
                update_list.emplace_back(n, denoise(queue.getKey(n)));
            queue.changeKey(n, apply_noise(key));
        }
        else
        {
            if (blow_ups != 0 && update_set.add(n))
                added_list.push_back(n);
            queue.insert(n, apply_noise(key));
        }
    }

    void update_queue(NodeID n);
};

struct reduction_ptr
{
    general_reduction *reduction = nullptr;

    reduction_ptr() = default;

    ~reduction_ptr()
    {
        release();
    }

    reduction_ptr(general_reduction *reduction) : reduction(reduction){};

    reduction_ptr(const reduction_ptr &other) : reduction(other.reduction->clone()){};

    reduction_ptr &operator=(const reduction_ptr &other)
    {
        release();
        reduction = other.reduction->clone();
        return *this;
    };

    reduction_ptr(reduction_ptr &&other) : reduction(std::move(other.reduction))
    {
        other.reduction = nullptr;
    };

    reduction_ptr &operator=(reduction_ptr &&other)
    {
        reduction = std::move(other.reduction);
        other.reduction = nullptr;
        return *this;
    };

    general_reduction *operator->() const
    {
        return reduction;
    }

    void release()
    {
        if (reduction)
        {
            delete reduction;
            reduction = nullptr;
        }
    };
};

template <class Last>
void make_reduction_vector_helper(std::vector<reduction_ptr> &vec, size_t n)
{
    vec.emplace_back(new Last(n));
};

template <class First, class Second, class... Redus>
void make_reduction_vector_helper(std::vector<reduction_ptr> &vec, size_t n)
{
    vec.emplace_back(new First(n));
    make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template <class... Redus>
std::vector<reduction_ptr> make_reduction_vector(size_t n)
{
    std::vector<reduction_ptr> vec;
    make_reduction_vector_helper<Redus...>(vec, n);
    return vec;
};

template <reduction_type type, int new_nodes>
reduction_ptr make_iterative_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_decreasing_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_plateau_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_increasing_struction(const ReductionConfig &config, size_t n);

#endif // REDUCTIONS_H
