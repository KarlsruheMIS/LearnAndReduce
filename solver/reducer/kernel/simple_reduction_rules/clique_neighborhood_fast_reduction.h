
#ifndef CLIQUE_NEIGHBORHOOD_FAST_REDUCTION_H
#define CLIQUE_NEIGHBORHOOD_FAST_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;

struct clique_neighborhood_reduction_fast : public general_reduction
{
    clique_neighborhood_reduction_fast(size_t n) : general_reduction(n) {}
    ~clique_neighborhood_reduction_fast() {}
    virtual clique_neighborhood_reduction_fast *clone() const final { return new clique_neighborhood_reduction_fast(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood_fast; }
    virtual std::string get_reduction_name() final { return "clique_nbh_fast"; }
    virtual bool reduce(reduce_algorithm *r_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *r_alg, NodeID v) final;
    int generate_data(reduce_algorithm *r_alg, NodeID v, std::vector<NodeID> &label);
};

#endif // CLIQUE_NEIGHBORHOOD_FAST_REDUCTION_H
