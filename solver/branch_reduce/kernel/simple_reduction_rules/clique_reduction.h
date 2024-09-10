
#ifndef CLIQUE_REDUCTION_H
#define CLIQUE_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class branch_and_reduce_algorithm;
struct clique_reduction : public general_reduction
{
    clique_reduction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~clique_reduction() {}
    virtual clique_reduction *clone() const final { return new clique_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::clique; }
    virtual std::string get_reduction_name() final { return "clique"; }
    virtual std::string get_model_path() final { return "models/clique.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final;

    int generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID> &label);

private:
    struct weighted_node
    {
        NodeID node;
        NodeWeight weight;
    };

    struct restore_data
    {
        weighted_node isolated;
        std::vector<NodeID> non_isolated;

        restore_data() = default;
        restore_data(const weighted_node &isolated, std::vector<NodeID> &&non_isolated) : isolated(isolated), non_isolated(std::move(non_isolated)){};
    };

    void fold(branch_and_reduce_algorithm *br_alg, const weighted_node &isolated, std::vector<NodeID> &&non_isolated);

    std::vector<restore_data> restore_vec;
};

#endif // CLIQUE_REDUCTION_H
