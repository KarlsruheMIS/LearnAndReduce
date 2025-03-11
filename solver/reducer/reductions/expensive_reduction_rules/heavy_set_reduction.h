
#ifndef HEAVY_SET_H
#define HEAVY_SET_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "reduction_config.h"
#include "fast_set.h"
#include "tiny_solver.h"

class reduce_algorithm;

struct heavy_set_reduction : public general_reduction
{
    heavy_set_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heavy_set_reduction() {}
    virtual heavy_set_reduction *clone() const final { return new heavy_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heavy_set; }
    virtual std::string get_reduction_name() final { return "heavy_set"; }
    virtual std::string get_model_path() final { return "models/heavy_set.lr_gcn"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(reduce_algorithm *br_alg, NodeID v) final;

private:
    enum v_combination
    {
        oo,
        uo,
        ov,
        uv
    }; // uo = u fixed as included and v is excluded (start with no increasing to all vertices)
    bool is_heavy_set(NodeID v, fast_set &v_neighbors_set, NodeID u, reduce_algorithm *br_alg);
    void set_weights(tiny_solver *solver, std::vector<NodeID> &nodes, std::vector<NodeWeight> &weights);
    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);
    void unset_weights(tiny_solver *solver, std::vector<NodeID> &nodes);
    bool compute_common_vertices = true;
};
#endif // HEAVY_SET_H
