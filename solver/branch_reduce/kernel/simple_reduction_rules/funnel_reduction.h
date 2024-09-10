
#ifndef FUNNEL_REDUCTION_H
#define FUNNEL_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"


class branch_and_reduce_algorithm;

struct funnel_reduction : public general_reduction
{
    funnel_reduction(size_t n) : general_reduction(n) {}
    ~funnel_reduction() {}
    virtual funnel_reduction *clone() const final { return new funnel_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::funnel; }
    virtual std::string get_reduction_name() final { return "funnel"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;
    int generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);

private:
    struct fold_data
    {
        NodeID node;
        NodeID funnel_neighbor;
    };

    struct restore_data
    {
        NodeID node;
        NodeID funnel_neighbor;
        std::vector<NodeID> remaining_neighbors;
        std::vector<std::vector<NodeID>> node_vecs;
    };

    bool is_funnel(NodeID v, NodeID &funnel_neighbor, branch_and_reduce_algorithm *br_alg, fast_set &funnel_set, std::vector<NodeID> &funnel_nodes);
    bool is_clique(branch_and_reduce_algorithm *br_alg, fast_set &clique_set, std::vector<NodeID> &clique_nodes);
    void fold(const fold_data &data, fast_set &funnel_set, branch_and_reduce_algorithm *br_alg);

    std::vector<restore_data> restore_vec;
};


#endif // FUNNEL_REDUCTION_H
