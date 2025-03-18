/**
 * reduction_config.h
 * Purpose: Configuration used for the reducttions for maximum independent set algorithms.
 *
 *****************************************************************************/

#pragma once

#include <string>
#include <cctype>

#include "definitions.h"
#include "config.h"

// Configuration for the calculation of the MIS
struct ReductionConfig : public Config {
    enum Struction_Type {ORIGINAL, MODIFIED, EXTENDED, EXTENDED_REDUCED, NONE};
    enum Backtrack_Type {IMMEDIATE_TIE_BREAKING, IMMEDIATE_EXCLUDE, END_MIN_KERNEL, NO_BACKTRACK};
    enum GNN_Filter_Type {NEVER, ALWAYS, INITIAL, INITIAL_TIGHT};
    enum Key_Type {RANDOM, DEGREE, INCREASE, APPROXIMATE_INCREASE};

    // Name of the kernel file.
    std::string kernel_filename;
    // disable single reductions
    bool disable_fold1=false;
    bool disable_fold2=false;
    bool disable_v_shape_min=false;
    bool disable_v_shape_mid=false;
    bool disable_v_shape_max=false;
    bool disable_v_shape=false;
    bool disable_triangle_mid = false;
    bool disable_triangle_min = false;
    bool disable_single_edge = false;
    bool disable_extended_single_edge = false;
    bool disable_twin= false;
    bool disable_simplicial_vertex = false;
    bool disable_clique_neighborhood= false;
    bool disable_clique_neighborhood_fast= false;
    bool disable_generalized_fold= false;
    bool disable_critical_set = false;
    bool disable_neighborhood= false;
    bool disable_domination = false;
    bool disable_extended_domination = false;
    bool disable_extended_domination_reverse = false;
    bool disable_extended_twin = false;
    bool disable_struction_decrease= false;
    bool disable_struction_plateau= false;
    bool disable_path= false;
    bool disable_cut_vertex = false;
    bool disable_high_degree = false;
    bool disable_funnel = false;
    bool disable_funnel_fold = false;
    bool disable_heavy_set = false;
    bool disable_heavy_set3 = false;
    bool disable_bound_reduction = false;
    bool disable_unconfined = false;
    bool disable_heuristic_include = false;
    bool disable_heuristic_exclude = false;
    bool use_hils_intersection = false;
    // bound for number of nodes in subgraphs
    NodeID subgraph_node_limit = 10;
    // Write the kernel into a file
    bool write_kernel;
    bool print_reduction_info = false;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
	// Perform reduction
	bool perform_reductions;
    // Choose reduction order and amount for given graph type
    std::string reduction_config_name = "full";
    std::string struction_config_name = "cyclicFast";
    double subgraph_time_limit=1;

    // early terminate solving subgraphs is best weight found is already to large for reduction to be applied
    bool disable_early_termination = false;

    // Type of filter reductions before adding to marker
    GNN_Filter_Type gnn_filter = GNN_Filter_Type::INITIAL_TIGHT;
    std::string gnn_filter_name = "initial tight";

    bool perform_hils;

    // Struction specific parameters
    Struction_Type struction_type;
    unsigned int struction_degree;
    unsigned int set_limit;
    double global_blow_up_factor;
    double phase_blow_up_factor;
    unsigned int phase_blow_ups;
    unsigned int max_unimproving_phases;
    Backtrack_Type backtrack_style;
    // StructionReduction_Style struction_reduction_style;
    Key_Type key_type;
    double key_reinsert_factor;
    bool reduce_and_peel;
    bool plain_struction;
    bool disable_blow_up;
    
    void setKeyType(const std::string & k_type) {
        if (strCompare(k_type, "random")) {
            key_type = Key_Type ::RANDOM;
        } else if (strCompare(k_type, "degree")) {
            key_type = Key_Type ::DEGREE;
        } else if (strCompare(k_type, "increase")) {
            key_type = Key_Type::INCREASE;
        } else if (strCompare(k_type, "approximate_increase")) {
            key_type = Key_Type ::APPROXIMATE_INCREASE;
        }
    }

    void setGNNFilterStyle(const std::string & filter) {
        if (strCompare(filter, "never")) {
            gnn_filter = GNN_Filter_Type::NEVER;
            std::string gnn_filter_name = "never";
        } else if (strCompare(filter, "always")) {
            gnn_filter = GNN_Filter_Type::ALWAYS;
            std::string gnn_filter_name = "always";
        } else if (strCompare(filter, "initial")) {
            gnn_filter = GNN_Filter_Type::INITIAL;
            std::string gnn_filter_name = "initial";
        } else if (strCompare(filter, "initial_tight")) {
            gnn_filter = GNN_Filter_Type::INITIAL_TIGHT;
            std::string gnn_filter_name = "initial tight";
        }
        else
        {
            gnn_filter = GNN_Filter_Type::INITIAL_TIGHT;
            std::string gnn_filter_name = "initial tight";
        }
        gnn_filter_name = filter;
    }

    void setBacktrackType(const std::string & back_type) {
        if (strCompare(back_type, "immediate_tie_breaking")) {
            backtrack_style = Backtrack_Type ::IMMEDIATE_TIE_BREAKING;
        } else if (strCompare(back_type, "immediate_exclude")) {
            backtrack_style = Backtrack_Type ::IMMEDIATE_EXCLUDE;
        } else if (strCompare(back_type, "end_min_kernel")) {
            backtrack_style = Backtrack_Type::END_MIN_KERNEL;
        } else if (strCompare(back_type, "no_backtrack")) {
            backtrack_style = Backtrack_Type ::NO_BACKTRACK;
        }
    }

    void setStructionType(const std::string & s_type) {
        if (strCompare(s_type, "original")) {
            struction_type = Struction_Type ::ORIGINAL;
        } else if (strCompare(s_type, "modified")) {
            struction_type = Struction_Type ::MODIFIED;
        } else if (strCompare(s_type, "extended")) {
            struction_type = Struction_Type::EXTENDED;
        } else if (strCompare(s_type, "extended_reduced")) {
            struction_type = Struction_Type ::EXTENDED_REDUCED;
        } else if (strCompare(s_type, "none")) {
            struction_type = Struction_Type ::NONE;
        }
    }
};

