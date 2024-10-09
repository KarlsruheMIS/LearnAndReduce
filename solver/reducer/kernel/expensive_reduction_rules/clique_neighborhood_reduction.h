
#ifndef CLIQUE_NEIGHBORHOOD_REDUCTION_H
#define CLIQUE_NEIGHBORHOOD_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;

struct clique_neighborhood_reduction : public general_reduction
{
    clique_neighborhood_reduction(size_t n) : general_reduction(n) {}
    ~clique_neighborhood_reduction() {}
    virtual clique_neighborhood_reduction *clone() const final { return new clique_neighborhood_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique_neighborhood; }
    virtual std::string get_reduction_name() final { return "clique_nbh"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;

    bool partition_into_cliques(NodeID v, reduce_algorithm *br_alg);
    bool expand_clique(NodeID max_neighbor, std::vector<NodeID> &neighbors_vec, fast_set &clique_neighbors_set, reduce_algorithm *br_alg);

    NodeWeight target_weight;
    NodeWeight neighbor_weights;
};

#endif // CLIQUE_NEIGHBORHOOD_REDUCTION_H
