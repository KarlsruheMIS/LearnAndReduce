
#ifndef EXTENDED_TWIN_REDUCTION_H
#define EXTENDED_TWIN_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;

struct extended_twin_reduction : public general_reduction
{
    extended_twin_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~extended_twin_reduction() {}
    virtual extended_twin_reduction *clone() const final { return new extended_twin_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::extended_twin; }
    virtual std::string get_reduction_name() final { return "extended_twin"; }
    virtual std::string get_model_path() final { return ""; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;

    int generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);

    struct restore_data
    {
        NodeID v;
        NodeID cand;
        NodeWeight candWeight;
    };
    size_t edge_count = 0;
    void fold(reduce_algorithm *br_alg, NodeID v, NodeID neighbor);

    std::vector<restore_data> restore_vec;
};

#endif // EXTENDED_TWIN_REDUCTION_H
