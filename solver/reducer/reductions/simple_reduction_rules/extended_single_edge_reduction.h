
#ifndef EXTENDED_SE_REDUCTION_H
#define EXTENDED_SE_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;

struct extended_single_edge_reduction : public general_reduction
{
    extended_single_edge_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~extended_single_edge_reduction() {}
    virtual extended_single_edge_reduction *clone() const final { return new extended_single_edge_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::extended_single_edge; }
    virtual std::string get_reduction_name() final { return "extended_single_edge"; }
    virtual std::string get_model_path() final { return ""; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;
    int generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label);
};

#endif // EXTENDED_SE_REDUCTION_H
