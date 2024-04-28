/**
 * reduction_evomis.cpp
 * Purpose: Main program for the evolutionary algorithm.
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
#include "struction_log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "cout_handler.h"
#include "reduction_config.h"
#include "parse_parameters.h"
#include "branch_and_reduce_algorithm.h"
#include "strongly_connected_components.h"
#include "solution_check.h"
#include "graph_operations.h"


int main(int argn, char **argv) {
    struction_log::instance()->restart_total_timer();
    struction_log::instance()->print_reduction_title();

    ReductionConfig config;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code) return 0;

    // read and set the graph
    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string path = graph_filepath.substr(0, graph_filepath.find_last_of('/'));
    std::string name = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    config.graph_filename = name;
    std::string path_and_file = graph_filepath.substr(0,graph_filepath.find_last_of('-'));
    struction_log::instance()->set_config(config);
    struction_log::instance()->print_config();

    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);
    struction_log::instance()->set_graph(G);
    struction_log::instance()->print_graph();

    branch_and_reduce_algorithm reducer(G, config);
    timer t;

    cout_handler ch;
    ch.enable_cout();
    graph_access &g = reducer.kernelize();
    double time = t.elapsed();

    ch.disable_cout();
    std::vector<int> comp_map(g.number_of_nodes(), 0);
    size_t comp_count = strongly_connected_components().strong_components(g, comp_map);
    std::vector<size_t> comp_size(comp_count, 0);
    for (auto comp_id : comp_map)
        comp_size[comp_id]++;
    size_t max_component = comp_count != 0 ? *std::max_element(comp_size.begin(), comp_size.end()) : 0;

    // cout_handler::enable_cout();
    ch.enable_cout();
    struction_log::instance()->print_full_reduction(config, time, reducer.get_current_is_weight(), g.number_of_nodes(), g.number_of_edges() / 2, comp_count, max_component);

    // std::cout << g.number_of_nodes() << "," << g.number_of_edges() / 2 << ",";
    // std::cout << comp_count << "," << max_component << ",";
    // std::cout << time << "," << reducer.get_current_is_weight() << std::endl;

    if (config.write_kernel && g.number_of_nodes() > 0)
    {   
        go.writeGraphWeighted(g, config.kernel_filename + ".kernel_graph");
        go.writeGraphWeighted_to_csv(g, config.kernel_filename + ".csv");
        go.writeWeights_to_csv(g, config.kernel_filename + "_weight.csv");
    }
    return 0;
}
