
#ifndef TWIN_REDUCTION_H
#define TWIN_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;

struct twin_reduction : public general_reduction
{
    twin_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
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

#endif // TWIN_REDUCTION_H
