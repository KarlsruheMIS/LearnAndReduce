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

bool write_reduction_data(std::vector<std::vector<bool>> &reduction_data, std::string filename) {
    int count_data_per_graph = 0;
    for (size_t i = 0; i < reduction_data.size(); i++) {
        // check at least 10% of applications
        if (std::count(reduction_data[i].begin(), reduction_data[i].end(), true) < 0.10 * reduction_data[i].size()) {
            continue;
        }
        count_data_per_graph++;
    }
    if (count_data_per_graph < 2) return false; // specify how much data per graph needed

    for (size_t i = 0; i < reduction_data.size(); i++) {
        // check at least 10% of applications
        if (std::count(reduction_data[i].begin(), reduction_data[i].end(), true) < 0.010 * reduction_data[i].size()) {
            continue;
        }

        std::cout << "Data is suitable for reduction " << i << std::endl;
        std::string reduction_filename = filename + "_reduction" + std::to_string(i) + ".txt";
        std::ofstream file;
        file.open(reduction_filename);
        for (int j = 0; j < reduction_data[i].size(); j++) {
            file << reduction_data[i][j] << std::endl;
        }
        file.close();
    }
    return true;
}

int main(int argn, char **argv) {
    struction_log::instance()->restart_total_timer();
    bool output_convergence = false;
    bool output_best_solution = true;
    struction_log::instance()->print_title();
    ReductionConfig config;

    // Parse the command line parameters;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config); 
    if (ret_code) return 0;

    config.generate_training_data = true;
    config.disable_critical_set = true;
    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string name = config.graph_filename.substr(0,config.graph_filename.find_last_of('.'));
    struction_log::instance()->set_config(config);
    struction_log::instance()->print_config();

    // Read the graph
    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);
    struction_log::instance()->set_graph(G);
    struction_log::instance()->print_graph();

    branch_and_reduce_algorithm reducer(G, config);

    int num_of_reductions = 15;

    if (config.size_of_subgraph > G.number_of_nodes()) {
        config.size_of_subgraph = G.number_of_nodes();
        config.num_of_subgraphs = 1;
    }
    std::vector<std::vector<bool>> reduction_data(num_of_reductions, std::vector<bool>(config.size_of_subgraph, false));
    graph_access subgraph; 
    for (int i = 0; i < config.num_of_subgraphs; i++) { // number of different subgraphs out of one graph
        reducer.get_training_data_for_graph_size(subgraph, config.size_of_subgraph, reduction_data); // size of the subgraph
        if(write_reduction_data(reduction_data, "training_data/reduction_data/" + name + "seed"+ std::to_string(config.seed) +"_training_data_" + std::to_string(i)))
            graph_io::writeGraphWeighted(subgraph, "training_data/graph/" + name + "seed"+ std::to_string(config.seed) +"_training_data_" + std::to_string(i) + ".graph");
        reduction_data = std::vector<std::vector<bool>>(num_of_reductions, std::vector<bool>(config.size_of_subgraph, false));
    }

    return 0;
}
