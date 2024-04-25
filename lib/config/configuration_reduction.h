/**
 * configuration.h
 * Purpose: Contains preset configurations for the evolutionary algorithms.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * H
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "definitions.h"
#include "reduction_config.h"

class configuration_reduction {
    public:
        /**
         * Default Constructor.
         */
        configuration_reduction() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_reduction() {};

        /**
         * Set the standard configuration.
         *
         * @param config Config to be initialized.
         */
        void standard( ReductionConfig & config );
        void disable_new_reductions( ReductionConfig & config );
        void enable_new_reductions( ReductionConfig & config );
        void original_cyclicFast( ReductionConfig & config );
        void original_cyclicStrong( ReductionConfig & config );
        void original_kamis( ReductionConfig & config );
        void mmwis( ReductionConfig & config );
        void all_reductions_cyclicFast( ReductionConfig & config );
        void all_reductions_cyclicStrong( ReductionConfig & config );
        void all_decreasing( ReductionConfig & config );
};




inline void configuration_reduction::standard( ReductionConfig & config ) {
    // Basic
    config.time_limit                             = 1000.0;
    config.reduction_time_limit                   = 1000.0;
    // Randomization
    config.seed                                   = 0;
    // Output
    config.console_log                            = false;
    config.check_sorted                           = true;
    // ILS
    config.ils_iterations                         = 15000;
    config.force_cand                             = 4;
	config.sort_freenodes                         = true;
    // Reductions
	config.perform_reductions                     = true;
    config.reduction_style                        = ReductionConfig::Reduction_Style::initial;
    // Weights
    config.weight_source                          = ReductionConfig::Weight_Source::FILE;

    config.set_limit                              = 1024;
    config.struction_degree                       = 256;
    config.struction_type                         = ReductionConfig::Struction_Type::NONE;
    config.key_type                               = ReductionConfig::Key_Type::APPROXIMATE_INCREASE;
    config.key_reinsert_factor                    = 2;
    config.global_blow_up_factor                  = 2;
    config.phase_blow_up_factor                   = 2;
    config.phase_blow_ups                         = 1;
    config.max_unimproving_phases                 = 100;
    config.backtrack_style                        = ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE;
    config.reduce_and_peel                        = false;
    config.disable_generalized_fold               = false;
    config.disable_clique_neighborhood            = true;
    config.disable_clique_neighborhood_fast       = false;
    config.subgraph_node_limit                    = 20;
    config.disable_blow_up                        = true;
    config.plain_struction                        = false;
    config.perform_hils                           = true;
    config.disable_domination                     = true;
}

inline void configuration_reduction::disable_new_reductions( ReductionConfig & config ) {
    config.disable_fold1                          = true;
    config.disable_v_shape_min                    = true;
    config.disable_v_shape_mid                    = true;
    config.disable_triangle                       = true;
    config.disable_basic_se                       = true;
    config.disable_extended_se                    = true;
    config.disable_domination                     = false;
    config.disable_cut_vertex                     = true;
    config.disable_funnel                         = true;
    config.disable_funnel_fold                    = true;
    config.disable_heavy_vertex                   = true;
    config.disable_heavy_set                      = true;
    config.disable_heavy_set3                     = true;
}

inline void configuration_reduction::enable_new_reductions( ReductionConfig & config ) {
    config.disable_domination                     = true;
    config.disable_fold1                          = false;
    config.disable_v_shape_min                    = false;
    config.disable_v_shape_mid                    = false;
    config.disable_triangle                       = false;
    config.disable_basic_se                       = false;
    config.disable_extended_se                    = false;
    config.disable_cut_vertex                     = false;
    config.disable_funnel                         = false;
    config.disable_funnel_fold                    = false;
    config.disable_heavy_vertex                   = false;
    config.disable_heavy_set                      = false;
    config.disable_heavy_set3                     = false;
}

inline void configuration_reduction::original_cyclicFast( ReductionConfig & config ) {
    standard(config);
    disable_new_reductions(config);
    config.disable_blow_up                        = false;
    config.struction_type                         = ReductionConfig::Struction_Type::EXTENDED;
    config.struction_reduction_style              = ReductionConfig::StructionReduction_Style::NORMAL;
    config.disable_generalized_fold               = true;
    config.disable_clique_neighborhood            = true;
    config.global_blow_up_factor                  = 9999;
    config.struction_degree                       = 64;
    config.max_unimproving_phases                 = 512;
    config.set_limit                              = 512;
}


inline void configuration_reduction::original_cyclicStrong( ReductionConfig & config ) {
    standard(config);
    disable_new_reductions(config);
    config.struction_type                         = ReductionConfig::Struction_Type::EXTENDED;
    config.struction_reduction_style              = ReductionConfig::StructionReduction_Style::NORMAL;
    config.disable_blow_up                        = false;
    config.disable_generalized_fold               = true;
    config.global_blow_up_factor                  = 9999;
    config.struction_degree                       = 512;
    config.max_unimproving_phases                 = 64;
    config.set_limit                              = 2048;
}

inline void configuration_reduction::original_kamis( ReductionConfig & config ) {
    standard(config);
    disable_new_reductions(config);
}

inline void configuration_reduction::mmwis( ReductionConfig & config ) {
    standard(config);
    config.disable_clique_neighborhood_fast       = false;
    config.disable_cut_vertex                     = true;
    config.disable_funnel                         = true;
}

inline void configuration_reduction::all_reductions_cyclicFast( ReductionConfig & config ) {
    original_cyclicFast(config);
    config.struction_reduction_style              = ReductionConfig::StructionReduction_Style::FULL;
    enable_new_reductions(config);
}

inline void configuration_reduction::all_reductions_cyclicStrong( ReductionConfig & config ) {
    original_cyclicStrong(config);
    config.struction_reduction_style              = ReductionConfig::StructionReduction_Style::FULL;
    enable_new_reductions(config);
}

inline void configuration_reduction::all_decreasing( ReductionConfig & config ) {
    standard(config);
}
