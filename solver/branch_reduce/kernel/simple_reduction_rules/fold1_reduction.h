
#ifndef FOLD1_REDUCTION_H
#define FOLD1_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;

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
    };

    void fold(branch_and_reduce_algorithm *br_alg, const fold_nodes &nodes);

    std::vector<restore_data> restore_vec;
};
#endif // FOLD1_REDUCTION_H
