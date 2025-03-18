/**
 * configuration.h
 * Purpose: Contains preset configurations for the reduction algorithms.
 *
 *****************************************************************************/

#pragma once

#include "definitions.h"
#include "reduction_config.h"

class configuration_reduction
{
public:
    /**
     * Default Constructor.
     */
    configuration_reduction() {};

    /**
     * Default Destructor.
     */
    virtual ~configuration_reduction() {};

    /**
     * Set the standard configuration.
     *
     * @param config Config to be initialized.
     */
    void standard(ReductionConfig &config);
    void disable_new_reductions(ReductionConfig &config);
    void disable_gnn_filter_reductions(ReductionConfig &config);
    void enable_new_reductions(ReductionConfig &config);

    void cyclicFast(ReductionConfig &config);
    void cyclicStrong(ReductionConfig &config);

    void all_reductions(ReductionConfig &config);
    void no_gnn_reductions(ReductionConfig &config);

    void generate_training_data_initial_reductions(ReductionConfig &config);
    void generate_training_data_expensive_reductions(ReductionConfig &config);
};

inline void configuration_reduction::standard(ReductionConfig &config)
{
    // Basic
    config.time_limit = 1000.0;
    // Randomization
    config.seed = 0;
    // Output
    config.console_log = false;
    config.check_sorted = true;
    config.sort_freenodes = true;
    // Reductions
    config.perform_reductions = true;
    // Weights
    config.weight_source = ReductionConfig::Weight_Source::FILE;

    config.set_limit = 1024;
    config.struction_degree = 256;
    config.struction_type = ReductionConfig::Struction_Type::EXTENDED;
    config.key_type = ReductionConfig::Key_Type::APPROXIMATE_INCREASE;
    config.key_reinsert_factor = 2;
    config.global_blow_up_factor = 2;
    config.phase_blow_up_factor = 2;
    config.phase_blow_ups = 1;
    config.max_unimproving_phases = 25;
    config.backtrack_style = ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    config.reduce_and_peel = false;
    config.disable_generalized_fold = false;
    config.disable_clique_neighborhood = true;
    config.subgraph_node_limit = 128;
    config.disable_blow_up = true;
    config.plain_struction = false;
    config.perform_hils = true;
    config.gnn_filter = ReductionConfig::GNN_Filter_Type::INITIAL_TIGHT;

    // struction
    cyclicFast(config);
    config.struction_config_name = "cyclicFast"; ; 
}

inline void configuration_reduction::disable_new_reductions(ReductionConfig &config)
{
    config.disable_fold1 = true;
    config.disable_v_shape_min = true;
    config.disable_v_shape_mid = true;
    config.disable_triangle_mid = true;
    config.disable_triangle_min = true;
    config.disable_single_edge = true;
    config.disable_extended_single_edge = true;
    config.disable_extended_domination = true;
    config.disable_extended_domination_reverse = true;
    config.disable_extended_twin = true;
    config.disable_cut_vertex = true;
    config.disable_funnel = true;
    config.disable_funnel_fold = true;
    config.disable_heavy_set = true;
    config.disable_heavy_set3 = true;
    config.disable_unconfined = true;
}

inline void configuration_reduction::enable_new_reductions(ReductionConfig &config)
{
    config.disable_fold1 = false;
    config.disable_v_shape_min = false;
    config.disable_v_shape_mid = false;
    config.disable_triangle_mid = false;
    config.disable_triangle_min = false;
    config.disable_single_edge = false;
    config.disable_extended_single_edge = false;
    config.disable_extended_domination = false;
    config.disable_extended_domination_reverse = false;
    config.disable_extended_twin = false;
    config.disable_cut_vertex = false;
    config.disable_funnel = false;
    config.disable_funnel_fold = false;
    config.disable_heavy_set = false;
    config.disable_heavy_set3 = false;
    config.disable_unconfined = false;
}

inline void configuration_reduction::disable_gnn_filter_reductions(ReductionConfig &config)
{
    config.disable_heavy_set = true;
    config.disable_heavy_set3 = true;
    config.disable_generalized_fold = true;
    config.disable_critical_set = true;
    config.disable_cut_vertex = true;
    config.disable_unconfined = true;
}

inline void configuration_reduction::cyclicFast(ReductionConfig &config)
{
    config.struction_type = ReductionConfig::Struction_Type::EXTENDED;
    config.key_type = ReductionConfig::Key_Type::APPROXIMATE_INCREASE;
    config.backtrack_style = ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    config.key_reinsert_factor = 2;
    config.phase_blow_up_factor = 2;
    config.phase_blow_ups = 1;
    config.reduce_and_peel = false;
    config.plain_struction = false;
    config.disable_blow_up = false;
    config.global_blow_up_factor = 9999;
    config.struction_degree = 64;
    config.max_unimproving_phases = 25;
    config.set_limit = 512;
}

inline void configuration_reduction::cyclicStrong(ReductionConfig &config)
{
    config.struction_type = ReductionConfig::Struction_Type::EXTENDED;
    config.key_type = ReductionConfig::Key_Type::APPROXIMATE_INCREASE;
    config.backtrack_style = ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    config.key_reinsert_factor = 2;
    config.phase_blow_up_factor = 2;
    config.phase_blow_ups = 1;
    config.reduce_and_peel = false;
    config.plain_struction = false;
    config.global_blow_up_factor = 9999;
    config.struction_degree = 512;
    config.max_unimproving_phases = 64;
    config.set_limit = 2048;
}

inline void configuration_reduction::all_reductions(ReductionConfig &config)
{
    standard(config);
    config.disable_generalized_fold = false;
    config.disable_clique_neighborhood_fast = false;
    enable_new_reductions(config);
}

inline void configuration_reduction::no_gnn_reductions(ReductionConfig &config)
{
    all_reductions(config);
    disable_gnn_filter_reductions(config);
}

inline void configuration_reduction::generate_training_data_initial_reductions(ReductionConfig &config)
{
    standard(config);
    config.disable_blow_up = true;
    enable_new_reductions(config);

    // not used for data generation on original instance
    config.disable_cut_vertex = true;
    config.disable_generalized_fold = true;
    config.disable_heavy_set = true;
    config.disable_heavy_set3 = true;

    // in general not used for data generation
    config.disable_extended_domination_reverse = true;
    config.disable_extended_domination = true;
    config.disable_funnel_fold = true;
    config.disable_struction_decrease = true;
    config.disable_struction_plateau = true;
    config.disable_clique_neighborhood = true;
    config.disable_clique_neighborhood_fast = true;
    config.disable_high_degree = true;
    config.disable_bound_reduction = true;
    config.subgraph_node_limit = 128;
    config.time_limit = 1000000;
}