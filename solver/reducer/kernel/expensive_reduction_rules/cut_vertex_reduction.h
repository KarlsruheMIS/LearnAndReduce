
#ifndef CUT_VERTEX_REDUCTION_H
#define CUT_VERTEX_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"

class reduce_algorithm;

struct cut_vertex_reduction : public general_reduction
{
    cut_vertex_reduction(size_t n) : general_reduction(n) {}
    ~cut_vertex_reduction() {}
    virtual cut_vertex_reduction *clone() const final { return new cut_vertex_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::cut_vertex; }
    virtual std::string get_reduction_name() final { return "cut_vertex"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual void restore(reduce_algorithm *br_alg) final;
    virtual void apply(reduce_algorithm *br_alg) final;

    void get_mappings_to_remaining_graph(reduce_algorithm *br_alg, std::vector<NodeID> &map, std::vector<NodeID> &reverse_map);
    bool get_fold_data(reduce_algorithm *br_alg, NodeID cut_v, std::vector<NodeID> &cut_v_included_i, std::vector<NodeID> &cut_v_included_e, std::vector<NodeID> &cut_v_excluded_i, std::vector<NodeID> &cut_v_excluded_e, NodeWeight &large_cutMWIS_weight, NodeWeight &small_cutMWIS_weight);
    void get_articulation_points(reduce_algorithm *br_alg, std::vector<NodeID> &articulation_points, std::vector<NodeID> &reverse_map, std::vector<NodeID> &map);
    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);

private:
    struct fold_data
    {
        NodeID cut_vertex;
        NodeWeight cut_vertex_weight;
        NodeWeight large_cutMWIS_weight;
        NodeWeight small_cutMWIS_weight;
        std::vector<NodeID> cut_component;
    };

    struct restore_data
    {
        fold_data data;
        std::vector<NodeID> case_cut_v_included_nodes_to_include;
        std::vector<NodeID> case_cut_v_included_nodes_to_exclude;
        std::vector<NodeID> case_cut_v_excluded_nodes_to_include;
        std::vector<NodeID> case_cut_v_excluded_nodes_to_exclude;
    };

    bool check_components(reduce_algorithm *br_alg, NodeID u, std::vector<NodeID> &smallComponent);
    bool build_small_component(NodeID u, reduce_algorithm *br_alg, std::vector<NodeID> &component, std::vector<bool> &component_visited);
    void dfs_fill_visited(NodeID u, reduce_algorithm *br_alg, std::vector<bool> &component_visited);
    void fold(reduce_algorithm *br_alg, fold_data &data, std::vector<NodeID> &cut_v_included_i, std::vector<NodeID> &cut_v_included_e, std::vector<NodeID> &cut_v_excluded_i, std::vector<NodeID> &cut_v_excluded_e);

    std::vector<restore_data> restore_vec;
};
#endif // CUT_VERTEX_REDUCTION_H
