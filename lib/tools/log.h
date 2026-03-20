/******************************************************************************
 * log.h
 *
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
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
         * Print information in one line for experiments.
         */
        void print_one_line_solution_data(ReductionConfig &mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m, NodeWeight solution, double full_time);

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

