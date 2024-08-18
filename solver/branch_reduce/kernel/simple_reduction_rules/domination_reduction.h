
#ifndef DOMINATION_REDUCTION_H
#define DOMINATION_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;

struct domination_reduction : public general_reduction
{
    domination_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~domination_reduction() {}
    virtual domination_reduction *clone() const final { return new domination_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::domination; }
    virtual std::string get_reduction_name() final { return "domination"; }
    virtual std::string get_model_path() final { return "models/domination.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
};

#endif // DOMINATION_REDUCTION_H
