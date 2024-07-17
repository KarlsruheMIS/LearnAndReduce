#pragma once
#include "reductions.h"

struct high_degree_reduction : public general_reduction
{
    high_degree_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~high_degree_reduction() {}
    virtual high_degree_reduction *clone() const final { return new high_degree_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::high_degree; }
    virtual std::string get_reduction_name() final { return "high_degree"; }
    virtual std::string get_model_path() final { return ""; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    bool reduce_in_component(branch_and_reduce_algorithm *br_alg, NodeID v, graph_access &G, NodeWeight &MWIS_weight);
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;
    bool solve_graph(NodeWeight &solution, graph_access &graph, ReductionConfig &config, bool apply_solution = false);

private:
    struct fold_nodes
    {
        NodeID main;
        std::vector<NodeID> MWIS;
        bool is_component;
    };

    struct restore_data
    {
        fold_nodes nodes;
        NodeWeight main_weight;
        NodeWeight MWIS_weight;
        std::vector<bool> component_nodes;
    };

    void fold(branch_and_reduce_algorithm *br_alg, fold_nodes &nodes, NodeWeight MWIS_weight, std::vector<NodeID> &excluded_nodes, std::vector<bool> &component_nodes);

    std::vector<restore_data> restore_vec;
};