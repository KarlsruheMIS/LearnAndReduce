//
// Created by alex on 25.05.20.
//

/**
 * reduction_evomis.cpp
 * Purpose: Main program for the evolutionary algorithm.
 *
 ******************************************************************************
 * Copyright (C) 2015-2017 Sebastian Lamm <lamm@ira.uka.de>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either
            std::cout << t.elapsed() << "," << min_kernel << std::endl;version 2 of the License, or (at your option)
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
#include "parse_parameters.h"
#include "branch_and_reduce_algorithm.h"
#include "strongly_connected_components.h"

int main(int argn, char **argv) {
    struction_log::instance()->restart_total_timer();

    ReductionConfig config;

    // Parse the command line parameters;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config); 
    if (ret_code) return 0;
    // int ret_code = setConfig(config); return 0;
    // Parse the command line parameters;
    // int ret_code = parse_struction_parameters(argn, argv, config, graph_filepath);
    // if (ret_code) {
    // }

    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    struction_log::instance()->set_config(config);

    // Read the graph
    graph_access G;
    graph_io::readGraphWeighted(G, graph_filepath);

    branch_and_reduce_algorithm reducer(G, config);
    std::cout << "time,graph_size" << std::endl;
    reducer.kernelize();


    return 0;
}
