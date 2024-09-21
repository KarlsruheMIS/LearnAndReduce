
#ifndef GENERAL_REDUCTION_H
#define GENERAL_REDUCTION_H

// local includes
#include "definitions.h"
#include "graph_access.h"
#include "reduction_config.h"
#include "vertex_marker.h"

// system includes
#include <iostream>
#include <iomanip> // For std::setw

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
    extended_domination,
    extended_domination_reverse,
    extended_twin,
    twin,
    critical_set,
    generalized_fold,
    single_edge,
    extended_single_edge,
    heavy_set,
    heavy_set3,
    cut_vertex,
    high_degree,
    funnel,
    funnel_fold,
    bound,
    unconfined_csr,
    struction_decrease,
    struction_plateau,
    struction_blow,
    heuristic_include,
    heuristic_exclude
};

class reduce_algorithm;

struct general_reduction
{
    general_reduction(size_t n) : marker(n) {}
    virtual ~general_reduction() {}
    virtual general_reduction *clone() const = 0;

    virtual reduction_type get_reduction_type() const = 0;
    virtual void print_reduction_type()
    {
        std::cout << std::left << std::setw(24) << get_reduction_name() << " \t";
    };
    virtual std::string get_reduction_name() = 0;
    virtual std::string get_model_path() { return ""; }
    virtual bool reduce(reduce_algorithm *br_alg) = 0;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) { return false; }
    virtual void restore(reduce_algorithm *br_alg) {}
    virtual void apply(reduce_algorithm *br_alg) {}
    virtual void reset(reduce_algorithm *br_alg, size_t comp_size) {}

    bool has_run = false;
    bool has_filtered_marker = false;
    NodeID reduced_nodes = 0;
    double reduction_time = 0.0;
    vertex_marker marker;

    NodeWeight get_neighborhood_weight(NodeID v, reduce_algorithm *br_alg);
    NodeID get_max_weight_neighbor(NodeID v, reduce_algorithm *br_alg);
    void get_neighborhood_set(NodeID v, reduce_algorithm *br_alg, fast_set &neighborhood_set);
    void get_neighborhood_vector(NodeID v, reduce_algorithm *br_alg, std::vector<NodeID> &neighborhood_vec);
    bool try_neighborhood_reduction(NodeID v, reduce_algorithm *br_alg, NodeWeight neighborhood_weight);
    bool solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, std::vector<NodeID> &nodes_vec, const fast_set &nodes_set);
    bool solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, NodeID v);
    bool is_reduced(NodeID v, reduce_algorithm *br_alg);
    virtual bool is_suited(NodeID v, reduce_algorithm *br_alg);
    // for generating training data
    bool solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, std::vector<NodeID> &nodes_vec, const fast_set &nodes_set, int &label);
    bool solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, NodeID v, int &label);
    virtual int generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label) { return 0; }
    virtual void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index) { return; }

    template <typename F>
    inline void for_each_changed_vertex(reduce_algorithm *br_alg, F f)
    {
        for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++)
        {
            NodeID v = marker.current_vertex(v_idx);

            if (v < br_alg->status.n && !is_reduced(v, br_alg))
            {
                f(v);
            }
        }
    }
};
#endif // GENERAL_REDUCTIONS_H
