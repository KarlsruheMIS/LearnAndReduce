
#ifndef GENERALIZED_FOLD_REDUCTION_H
#define GENERALIZED_FOLD_REDUCTION_H

// local includes
#include "definitions.h"
#include "dynamic_graph.h"
#include "general_reduction.h"
#include "fast_set.h"

class branch_and_reduce_algorithm;
struct generalized_fold_reduction : public general_reduction
{
    generalized_fold_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~generalized_fold_reduction() {}
    virtual generalized_fold_reduction *clone() const final { return new generalized_fold_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::generalized_fold; }
    virtual std::string get_reduction_name() final { return "generalized_fold"; }
    virtual std::string get_model_path() final { return "models/heavy_vertex.gnn"; }
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
#endif // GENERALIZED_FOLD_REDUCTIONS_H
