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


int main(int argn, char **argv) {

    ReductionConfig config;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code) return 0;

    if (config.output_filename.empty()) {
        std::cerr << "Output filename is empty." << std::endl;
        return 0;
    }
    // read and set the graph
    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string path = graph_filepath.substr(0, graph_filepath.find_last_of('/'));
    std::string name = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    std::string path_and_file = graph_filepath.substr(0,graph_filepath.find_last_of('-'));

    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);
    // mis_log::instance()->set_graph(G);
    // mis_log::instance()->print_graphj);

    go.writeGraphWeighted_to_csv(G, config.output_filename + ".csv");
    go.writeWeights_to_csv(G, config.output_filename + "_weight.csv");
    return 0;
}
