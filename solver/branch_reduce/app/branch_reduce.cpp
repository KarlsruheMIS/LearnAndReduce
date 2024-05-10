/******************************************************************************
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

#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <argtable3.h>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>

#include "timer.h"
#include "cout_handler.h"
#include "struction_log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "reduction_config.h"
#include "parse_parameters.h"
#include "branch_and_reduce_algorithm.h"
#include "graph_operations.h"
#include "solution_check.h"


int main(int argn, char **argv) {
    bool output_convergence = false;
    bool output_best_solution = true;
    struction_log::instance()->print_title();
    ReductionConfig config;

    // Parse the command line parameters;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config); 
    if (ret_code) return 0;

    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string name = config.graph_filename.substr(0,config.graph_filename.find_last_of('.'));
    config.graph_filename = name;
    struction_log::instance()->set_config(config);
    struction_log::instance()->print_config();

    // Read the graph
    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);
    struction_log::instance()->set_graph(G);
    struction_log::instance()->print_graph();
    struction_log::instance()->restart_total_timer();


    branch_and_reduce_algorithm reducer(G, config);
    auto start = std::chrono::system_clock::now();

    reducer.run_branch_reduce();

    auto end = std::chrono::system_clock::now();

    std::chrono::duration<float> branch_reduce_time = end - start;

    reducer.apply_branch_reduce_solution(G);

    NodeWeight is_weight = 0;
    solution_check<graph_access> sc(G);
    if (!sc.check_IS()) {
        std::cerr << "ERROR: graph after inverse reduction is not independent" << std::endl;
        exit(1);
    } else {
        forall_nodes(G, node) {
            if (G.getPartitionIndex(node) == 1) {
                is_weight += G.getNodeWeight(node);
            }
        } endfor
    }

    struction_log::instance()->print_reduction(config, reducer.kernelization_time, reducer.kernelization_offset, reducer.min_kernel, reducer.max_min_kernel_comp);
    struction_log::instance()->print_results(!reducer.timeout && reducer.get_heuristically_reduced_vertices() == 0);

#ifndef OUTPUT_WEIGHT_CONVERGENCE
    std::cout << name << "," << G.number_of_nodes() << "," << G.number_of_edges() << ","
              << is_weight << "," << reducer.best_is_time << "," << branch_reduce_time.count() << ","
              << reducer.timeout << ","
              << reducer.min_kernel << "," << reducer.max_min_kernel_comp << "," << reducer.kernelization_time << std::endl;
#endif

    return 0;
}
