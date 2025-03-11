/******************************************************************************
 * log.h
 *
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
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <sstream>

#include "timer.h"
#include "reduction_config.h"
#include "data_structure/graph_access.h"

class log {
    public:
        /**
         * Get the singleton logger instance.
         * 
         * @return Instance of the logger.
         */
        static log *instance() {
            static log inst;
            return &inst;
        };

        /**
         * Set the config.
         *
         * @param config Config for the evolutionary algorithm.
         */
        void set_config(ReductionConfig & config);

        /**
         * Set the graph.
         *
         * @param G Graph representation.
         */
        void set_graph(graph_access & G);

        /**
         * Write the log to a file.
         */
        void write_log();

        /**
         * Print information about the graph.
         */
        void print_graph();


        /**
         * Print information in one line for experiments.
         */
        void print_one_line_kernel_data(ReductionConfig &mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m);

        /**
         * Print the current config.
         */
        void print_config();

        /**
         * Print information about the reduction step.
         * Includes number of extracted nodes and resulting kernel size
         *
         * @param mis_config Config for the logger.
         * @param time Time needed for reduction.
         * @param offset Reduction offset  
         * @param kernel_size_n Number of remaining nodes.
         * @param kernel_size_m Number of remaining edges.
         * @param component_count Number of connected components in the kernel.
         * @param max_component Size of the largest connected component in the kernel.
         */
        void print_reduction(ReductionConfig &mis_config, double time, NodeWeight offset, NodeID kernel_size_n, size_t max_component);
        void print_full_reduction(ReductionConfig &mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m, size_t component_count, size_t max_component);

        /**
         * Print the final results.
         * @param optimal Is the computed solution optimal.
         */
        void print_full_reduction(ReductionConfig &mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m);
        void print_results(bool optimal);

        /**
         * Print a title.
         */
        // void print_title();
        void print_data_generation_title();
        void print_reduction_title();
        
        /**
         * Restart the timer for the total time including IO, etc.
         */
        void restart_total_timer();

        /**
         * Update the size of the best solution.
         * @param size Candidate to replace the best solution size.
         */
        void set_best_size(unsigned int size);

        /**
         * Update the size of the best solution.
         * @param time Candidate to replace the best solution time.
         */
        void set_best_time(double time);

        void set_best(unsigned int size, double time);

        /**
         * Reset the size of the best solution.
         */
        void reset_best_size();

    private:
        // General information
        timer total_timer;
        std::stringstream filebuffer_string;
        ReductionConfig log_config;

        // Graph informations
        std::string graph_name;
        unsigned int number_of_nodes;
        unsigned int number_of_edges;
        unsigned int arc_scans;
        double avg_degree;
        double density;

        // Reduction information
        unsigned int current_is_size;
        unsigned int optimum_size;

        // Results information
        double total_time_taken;
        double time_taken_best;

        /**
         * Default Constructor.
         */
        log();

        /**
         * Default Destructor.
         */
        virtual ~log();
};

