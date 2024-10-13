#ifndef HEAVY_SET_3_H
#define HEAVY_SET_3_H

// local includes
#include "definitions.h"
#include "dynamic_graph.h"
#include "fast_set.h"
#include "general_reduction.h"
#include "reduction_config.h"
#include "tiny_solver.h"

// system includes
#include <array>
#include <memory>
#include <vector>

struct heavy_set3_reduction : public general_reduction
{
    heavy_set3_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heavy_set3_reduction() {}
    virtual heavy_set3_reduction *clone() const final { return new heavy_set3_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set3; }
    virtual std::string get_reduction_name() final { return "heavy_set3"; }
    virtual std::string get_model_path() final { return "models/heavy_set3.lr_gcn"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;

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

    bool is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, fast_set &u_neighbors_set, NodeID w, reduce_algorithm *br_alg);
    bool check_u_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_v_combination(std::vector<NodeWeight> &MWIS_weights);
    bool check_w_combination(std::vector<NodeWeight> &MWIS_weights);
    void unset_weights(tiny_solver *solver, std::vector<NodeID> &nodes);
    void set_weights(tiny_solver *solver, std::vector<NodeID> &nodes, std::vector<NodeWeight> &weights);

    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);
};
#endif // HEAVY_SET_3_H
