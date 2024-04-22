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
#include "graph_access.h"
#include "graph_io.h"
#include "cout_handler.h"
#include "reduction_config.h"
#include "parse_parameters.h"
#include "solution_check.h"
#include "graph_operations.h"

bool reduction_file_exists(const std::string & filename) {
    std::ifstream file(filename);
    return file.good();
}
void parse_reduction_data(const std::string & reduction_file, std::vector<bool> & reduction_data) {
    // each row has a 1 or 0 to be stored in reduction_data
    std::ifstream file(reduction_file);
    std::string line;
    while (std::getline(file, line)) {
        reduction_data.push_back(std::stoi(line));
    }
}
int write_weights_and_reductions_to_csv(graph_access & G, const std::string & resulting_filename, const std::string & filename) {
    // if reduction file with filename exists, add the content, otherwise all 0
    int num_of_reductions = 15;
    //reduction names:
    std::vector<std::string> reduction_names = {"fold1", "neighborhood", "fold2", "clique", "funnel", "single_edge", "extended_single_edge", "twin", "clique_nbh_fast", "heavy_vertex", "heavy_set", "generalized_fold", "cut_vertex"};
    std::vector<std::vector<bool>> reduction_data(num_of_reductions);
    for (size_t i = 0; i < num_of_reductions; i++) {
        std:: string reduction_file = "training_data/reduction_data/" +filename + "_reduction" + std::to_string(i) + ".txt";
        if (reduction_file_exists(reduction_file)) {
            parse_reduction_data(reduction_file, reduction_data[i]);
        } else {
            reduction_data[i] = std::vector<bool>(G.number_of_nodes(), false);
        }
    }

    std::ofstream f(resulting_filename.c_str());
    f << "id;node_weight";
    for (size_t i = 0; i < num_of_reductions; i++) {
        f << ";" << reduction_names[i];
    }
    f << std::endl;

    forall_nodes(G, node) {
        f << node << ";" << G.getNodeWeight(node);
        for (size_t i = 0; i < num_of_reductions; i++) {
            f << ";" << reduction_data[i][node];
        }
        f << std::endl;
    } endfor

    f.close();
    return 0;
}

int main(int argn, char **argv) {

    ReductionConfig config;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code) return 0;

    if (config.output_filename.empty()) {
        config.output_filename = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    }
    // read and set the graph
    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string path = graph_filepath.substr(0, graph_filepath.find_last_of('/'));
    std::string name = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    std::string path_and_file = graph_filepath.substr(0,graph_filepath.find_last_of('-'));

    config.check_sorted = false;
    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);


    go.writeGraphWeighted_to_csv(G, config.output_filename + ".csv");
    // go.writeWeights_to_csv(G, config.output_filename + "_weight.csv");
    write_weights_and_reductions_to_csv(G, config.output_filename + "_weights_and_reduction.csv" , name);
    return 0;
}