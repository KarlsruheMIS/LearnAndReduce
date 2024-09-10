
#ifndef SINDLE_EDGE_REDUCTION_H
#define SINDLE_EDGE_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;
struct single_edge_reduction : public general_reduction
{
    single_edge_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~single_edge_reduction() {}
    virtual single_edge_reduction *clone() const final { return new single_edge_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::single_edge; }
    virtual std::string get_reduction_name() final { return "single_edge"; }
    virtual std::string get_model_path() final { return "models/single_edge.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    int generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label);
};

#endif // SINDLE_EDGE_REDUCTION_H
