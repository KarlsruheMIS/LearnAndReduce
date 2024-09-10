
#ifndef FOLD2_REDUCTION_H
#define FOLD2_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"
#include "dynamic_graph.h"

#include <array>

class branch_and_reduce_algorithm;

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

    int generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);

private:
    struct fold_nodes
    {
        NodeID deg2_node;
        std::vector<NodeID> neighbors;
    };
    enum fold2_case
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
        fold2_case fold_case;
        dynamic_graph::neighbor_list neighbor_list;
        std::array<std::vector<NodeID>, 2> node_vecs;
    };

    void fold_triangle_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_triangle_min_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_max_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_mid_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);
    void fold_v_shape_min_weight(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);

    size_t v_shape_min_count = 0;

    std::vector<restore_data> restore_vec;
};
#endif // FOLD2_REDUCTION_H
