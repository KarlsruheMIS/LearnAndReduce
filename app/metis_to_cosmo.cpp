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
int write_weights_and_reductions_to_csv(graph_access & G, const std::string & resulting_filename, const std::string & filename, ReductionConfig & config) {
    // if reduction file with filename exists, add the content, otherwise all 0
    //reduction names:
     // {"fold1", "neighborhood", "fold2", "clique", "funnel", "funnel_fold", "single_edge", "extended_single_edge", "twin", "clique_nbh_fast", "clique_neighborhood", "decreasing_struction", "heavy_vertex", "generalized_fold", "heavy_set", "heavy_set3", "cut_vertex"};
    std::vector<std::string> reduction_names;
    if (!config.disable_fold1) reduction_names.push_back("fold1");
    if (!config.disable_neighborhood) reduction_names.push_back("neighborhood");
    if (!config.disable_fold2) reduction_names.push_back("fold2");
    if (!config.disable_clique) reduction_names.push_back("clique");
    if (!config.disable_funnel) reduction_names.push_back("funnel");
    if (!config.disable_funnel_fold) reduction_names.push_back("funnel_fold");
    if (!config.disable_domination) reduction_names.push_back("domination");
    if (!config.disable_basic_se) reduction_names.push_back("single_edge");
    if (!config.disable_extended_se) reduction_names.push_back("extended_single_edge");
    if (!config.disable_twin) reduction_names.push_back("twin");
    if (!config.disable_clique_neighborhood_fast) reduction_names.push_back("clique_nbh_fast");
    if (!config.disable_clique_neighborhood) reduction_names.push_back("clique_neighborhood");
    if (!config.disable_decreasing_struction) reduction_names.push_back("decreasing_struction");
    if (!config.disable_generalized_fold) reduction_names.push_back("generalized_fold");
    if (!config.disable_heavy_set) reduction_names.push_back("heavy_set");
    if (!config.disable_heavy_set3) reduction_names.push_back("heavy_set3");
    if (!config.disable_cut_vertex) reduction_names.push_back("cut_vertex");
    if (!config.disable_cut_vertex) reduction_names.push_back("component");

    int num_of_reductions = reduction_names.size();
    std::vector<std::vector<bool>> reduction_data(num_of_reductions);
    for (int i = 0; i < num_of_reductions; i++) {
        std:: string reduction_file = "training_data/reduction_data/" +filename + "_" + reduction_names[i] + ".txt";
        if (reduction_file_exists(reduction_file)) {
            parse_reduction_data(reduction_file, reduction_data[i]);
        } else {
            reduction_data[i] = std::vector<bool>(G.number_of_nodes(), false);
        }
    }

    std::ofstream f(resulting_filename.c_str());
    f << "id;node_weight";
    for (int i = 0; i < num_of_reductions; i++) {
        f << ";" << reduction_names[i];
    }
    f << std::endl;

    forall_nodes(G, node) {
        f << node << ";" << G.getNodeWeight(node);
        for (int i = 0; i < num_of_reductions; i++) {
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
    write_weights_and_reductions_to_csv(G, config.output_filename + "_weights_and_reduction.csv" , name, config);
    return 0;
}
