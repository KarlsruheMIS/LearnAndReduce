
#ifndef NEIGHBORHOOD_REDUCTION_H
#define NEIGHBORHOOD_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"

class reduce_algorithm;

struct neighborhood_reduction : public general_reduction
{
    neighborhood_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = true;
    }
    ~neighborhood_reduction() {}
    virtual neighborhood_reduction *clone() const final { return new neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::neighborhood; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;
    virtual std::string get_reduction_name() final { return "neighborhood"; }

    bool is_suited(NodeID v, reduce_algorithm *br_alg);
    int generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);
};

#endif // NEIGHBORHOOD_REDUCTION_H
