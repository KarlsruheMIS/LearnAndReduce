/**
 * reduction_config.h
 * Purpose: Configuration used for the evolutionary maximum independent set algorithms.
 *
 *****************************************************************************/

#pragma once

#include <string>
#include <cctype>

#include "definitions.h"
#include "config.h"

// Configuration for the calculation of the MIS
struct ReductionConfig : public Config {
    enum Reduction_Style {initial, time_ordering, weight_ordering, time_and_weight_ordering};
    enum Struction_Type {ORIGINAL, MODIFIED, EXTENDED, EXTENDED_REDUCED, NONE};
    enum Backtrack_Type {IMMEDIATE_TIE_BREAKING, IMMEDIATE_EXCLUDE, END_MIN_KERNEL, NO_BACKTRACK};
    enum StructionReduction_Style {NORMAL, DENSE};
    enum Key_Type {RANDOM, DEGREE, INCREASE, APPROXIMATE_INCREASE};

    // Name of the kernel file.
    std::string kernel_filename;
    std::string kernel_csv_filename;
    // disable single reductions
    bool disable_fold1=false;
    bool disable_fold2=false;
    bool disable_v_shape_min=false;
    bool disable_v_shape_mid=false;
    bool disable_v_shape_max=false;
    bool disable_v_shape=false;
    bool disable_triangle = false;
    bool disable_basic_se = false;
    bool disable_extended_se = false;
    bool disable_twin= false;
    bool disable_clique = false;
    bool disable_clique_neighborhood= false;
    bool disable_clique_neighborhood_fast= false;
    bool disable_generalized_fold= false;
    bool disable_generalized_neighborhood = false;
    bool disable_critical_set = false;
    bool disable_neighborhood= false;
    bool disable_heavy_set = false;
    bool disable_domination = true;
    bool disable_decreasing_struction= false;
    bool disable_path= false;
    bool disable_cut_vertex = false;
    bool disable_funnel = false;
    // bound for number of nodes in heavy_set neighborhood graph =0 is disabled reduction completely
    int  heavy_set= 20;
    int  cut_vertex_max_component_size= 20;
    // Write the kernel into a file
    bool write_kernel;
    bool print_reduction_info = false;
    // Apply all reductions to reduce the graph size
    bool all_reductions;
	// Perform reduction
	bool perform_reductions;
    // Threshold for reductions
    int reduction_threshold;
    // Choose reduction order and amount for given graph type
    Reduction_Style reduction_style;
    std::string reduction_config_name = "decreasing";
    std::string reduction_style_name = "initial";
    double reduction_time_limit;

    // Lower bound for the change rate between best individuals
    double best_limit;
    // Number of candidates for forced insertion.
    unsigned int force_cand;

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
    StructionReduction_Style struction_reduction_style;
    Key_Type key_type;
    double key_reinsert_factor;
    bool reduce_and_peel;
    bool plain_struction;
    bool disable_blow_up;

    void setReductionStyle(const std::string & redu_style) {
        if (strCompare(redu_style, "time"))                    { 
            reduction_style = Reduction_Style::time_ordering;
            reduction_style_name = "time";
        } else if (strCompare(redu_style, "weight"))           { 
            reduction_style = Reduction_Style::weight_ordering;
            reduction_style_name = "weight";
        } else if (strCompare(redu_style, "time_and_weight"))  { 
            reduction_style = Reduction_Style::time_and_weight_ordering;
            reduction_style_name = "time_and_weight";
        } else { 
            reduction_style = Reduction_Style::initial;
            reduction_style_name = "initial";
            }
    }
 
  void setStructionReductionType(const std::string & r_type) {
        if (strCompare(r_type, "dense")) {
            struction_reduction_style = StructionReduction_Style::DENSE;
        } else  {
            struction_reduction_style = StructionReduction_Style::NORMAL;
        }
    }
    
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

