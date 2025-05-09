
#ifndef FUNNEL_FOLD_REDUCTION_H
#define FUNNEL_FOLD_REDUCTION_H
// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;
struct funnel_fold_reduction : public general_reduction
{
    funnel_fold_reduction(size_t n) : general_reduction(n) {}
    ~funnel_fold_reduction() {}
    virtual funnel_fold_reduction *clone() const final { return new funnel_fold_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::funnel_fold; }
    virtual std::string get_reduction_name() final { return "funnel fold"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(reduce_algorithm *br_alg) final;
    virtual void apply(reduce_algorithm *br_alg) final;
    int generate_data(reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);

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

    bool is_funnel(NodeID v, NodeID &funnel_neighbor, reduce_algorithm *br_alg, fast_set &funnel_set, std::vector<NodeID> &funnel_nodes);
    bool is_clique(reduce_algorithm *br_alg, fast_set &clique_set, std::vector<NodeID> &clique_nodes);
    void fold(const fold_data &data, fast_set &funnel_set, reduce_algorithm *br_alg);

    std::vector<restore_data> restore_vec;
};

#endif // FUNNEL_FOLD_REDUCTION_H
