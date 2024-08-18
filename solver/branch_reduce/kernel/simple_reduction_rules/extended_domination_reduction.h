
#ifndef EXTENDED_DOMINATION_REDUCTION_H
#define EXTENDED_DOMINATION_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;

struct extended_domination_reduction : public general_reduction
{
    extended_domination_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~extended_domination_reduction() {}
    virtual extended_domination_reduction *clone() const final { return new extended_domination_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::extended_domination; }
    virtual std::string get_reduction_name() final { return "extended_domination"; }
    virtual std::string get_model_path() final { return ""; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

    struct restore_data
    {
        NodeID v;
        NodeID neighbor;
        NodeWeight neighborWeight;
    };
    size_t edge_count = 0;
    void fold(branch_and_reduce_algorithm *br_alg, NodeID v, NodeID neighbor);

    std::vector<restore_data> restore_vec;
};

#endif // EXTENDED_DOMINATION_REDUCTION_H
