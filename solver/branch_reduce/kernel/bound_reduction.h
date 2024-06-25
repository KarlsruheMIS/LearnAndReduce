#pragma once
#include "reductions.h"


struct bound_reduction : public general_reduction
{
    bound_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = true;
    }
    ~bound_reduction() {}
    virtual bound_reduction *clone() const final { return new bound_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::bound; }
    virtual std::string get_reduction_name() final { return "bound"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

    private:
    void build_component_graphs(graph_access &G, ReductionConfig &config);
    NodeWeight solve_partition(size_t i, ReductionConfig &config, bool apply_solution = false);

    // recompute partition weight for node beeing excluded
    NodeWeight recompute_partition_weight_including_node(size_t i, NodeID node, ReductionConfig &config);
    NodeWeight recompute_partition_weight_excluding_node(size_t i, NodeID node, ReductionConfig &config);

    NodeWeight recompute_partition_weight_excluding_nodes(size_t i, std::vector<NodeID> &nodes, ReductionConfig &config);

    void clear_partition_graphs();

    std::vector<graph_access*> partition_graphs;
    std::vector<std::vector<NodeID>> partition_mappings;
    std::vector<NodeWeight> partition_solution_weights;
    int k = 2;
};
