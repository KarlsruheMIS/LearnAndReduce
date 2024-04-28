/**
 * configuration.h
 * Purpose: Contains preset configurations for the maximum weight indpeneden set algorithms.
 *
 *****************************************************************************/

#pragma once

#include "definitions.h"
#include "config.h"
// #include "kaHIP_interface.h"

#ifdef OMIS
const int FAST           = 0;
const int ECO            = 1;
const int STRONG         = 2;
const int FASTSOCIAL     = 3;
const int ECOSOCIAL      = 4;
const int STRONGSOCIAL   = 5;
#endif


class configuration_mis {
    public:
        /**
         * Default Constructor.
         */
        configuration_mis() {} ;

        /**
         * Default Destructor.
         */
        virtual ~configuration_mis() {};

        void standard( Config & config );

};

inline void configuration_mis::standard( Config & mis_config ) {
    // Basic
    mis_config.time_limit                             = 1000.0;
    mis_config.seed                                   = 0;

    // Weights
    mis_config.weight_source                          = Config::Weight_Source::FILE;

    // Output
    mis_config.console_log                            = false;
    mis_config.check_sorted                           = true;

    // ILS
    mis_config.ils_iterations                         = 15000;
    /* mis_config.ils_time_limit                         = 15000; */
    mis_config.local_search_threshold                 = mis_config.ils_iterations;
    mis_config.force_k                                = 1;
    mis_config.force_cand                             = 4;
    mis_config.sort_freenodes                         = true;
    // Initial solution
    mis_config.start_greedy_adaptive                  = false;

}

