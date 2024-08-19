#ifndef HEAVY_SET_3_H 
#define HEAVY_SET_3_H 

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "dynamic_graph.h"
#include "reduction_config.h"

// system includes
#include <vector>
#include <memory>
#include <array>

struct heavy_set3_reduction : public general_reduction
{
    heavy_set3_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heavy_set3_reduction() {}
    virtual heavy_set3_reduction *clone() const final { return new heavy_set3_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set3; }
    virtual std::string get_reduction_name() final { return "heavy_set3"; }
    virtual std::string get_model_path() final { return "~/projects/MWIS_learn_and_reduce/models/heavy_set3.gnn"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;

private:
    enum v_combination
    {
        ooo,
        uoo,
        ovo,
        oow,
        uvo,
        uow,
        ovw,
        uvw
    }; // uoo = u fixed as included and v, w (=o) are excluded (start with no increasing to all vertices)

    bool is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, fast_set &u_neighbors_set, NodeID w, branch_and_reduce_algorithm *br_alg);
    bool check_u_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_v_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_w_combination(std::vector<NodeWeight> &MWIS_weights);
    void unset_weights(graph_access &graph, std::vector<NodeID> &nodes, std::vector<NodeID> &reverse_mapping);
    void set_weights(graph_access &graph, std::vector<NodeID> &nodes, std::vector<NodeID> &reverse_mapping, std::vector<NodeWeight> &weights);
};
#endif // HEAVY_SET_3_H
