/******************************************************************************
 * mmwis_log.h
 *
 *****************************************************************************/
#pragma once

#include <sstream>

#include "definitions.h"
#include "timer.h"
#include "reduction_config.h"
#include "mis_log.h"
#include "data_structure/graph_access.h"

class reduction_log : public mis_log {
    public:
        reduction_log();
        ~reduction_log();

        /**
         * Get the singleton logger instance.
         * 
         * @return Instance of the logger.
         */
        static reduction_log *instance() {
            static reduction_log inst;
            return &inst;
        };

        /**
         * Print the current config.
         */
        void print_config();

        /**
         * Print information about the reduction step.
         * Includes number of extracted nodes and resulting kernel size
         *
         * @param mis_config Config for the logger.
         * @param extracted_nodes Number of removed nodes.
         * @param kernel_size Number of remaining nodes.
         */
        void print_reduction(ReductionConfig & mis_config, NodeWeight offset, NodeID kernel_size);

        /**
         * Print a title.
         */
        void print_title();

    protected:
        // General information
        ReductionConfig log_config;

        // Reduction information
        NodeID kernel_size;
        NodeWeight offset;

};
