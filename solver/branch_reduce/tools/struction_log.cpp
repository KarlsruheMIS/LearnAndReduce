/******************************************************************************
 * struction_log.cpp
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

#include "struction_log.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

struction_log::struction_log() {
    number_of_nodes = 0;
    number_of_edges = 0;
    avg_degree = 0.0;
    density = 0.0;
    arc_scans = 0;
    total_time_taken = 0.0;
    time_taken_best = 0.0;
    optimum_size = 0;
}

struction_log::~struction_log() { }

void struction_log::set_config(ReductionConfig & config) {
    log_config = config; 
}

void struction_log::set_graph(graph_access & G) {
    number_of_nodes = G.number_of_nodes();
    number_of_edges = G.number_of_edges();
    avg_degree = (double) number_of_edges / number_of_nodes;
    density = (double) (2 * number_of_edges) / (number_of_nodes * (number_of_nodes - 1));
}

void struction_log::write_log() {
    std::stringstream filename_stream;
    filename_stream << "./logs/log_"<<  log_config.graph_filename <<   
                       "_seed_" <<  log_config.seed;
    std::ofstream f(filename_stream.str());
    f << filebuffer_string.str();
    f.close();
}

// void struction_log::print_newline() {
//     filebuffer_string << std::endl; 

//     if (log_config.console_log) {
//         std::cout << std::endl; 
//     }
// }
void struction_log::print_title() {
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "\t WeightedMIS Branch&Reduce"                                         << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;

    std::cout << "=========================================="                                   << std::endl;
    std::cout << "\t WeightedMIS Branch&Reduce"                                                 << std::endl;
    std::cout << "=========================================="                                   << std::endl;
}

void struction_log::print_reduction_title() {
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "\t WeightedMIS Reductions"                                            << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;

    std::cout << "=========================================="                                   << std::endl;
    std::cout << "\t WeightedMIS Reductions"                                                    << std::endl;
    std::cout << "=========================================="                                   << std::endl;
}

void struction_log::print_graph() {
    filebuffer_string << "\t\tGraph"                                                            << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "IO time:\t\t\t\t"         << total_timer.elapsed()                    << std::endl;
    filebuffer_string << "Filename:\t\t\t\t"        << log_config.graph_filename                << std::endl;
    filebuffer_string << "|-Nodes:\t\t\t\t"         << number_of_nodes                          << std::endl;
    filebuffer_string << "|-Edges:\t\t\t\t"         << number_of_edges                          << std::endl;
    filebuffer_string << std::endl;

    std::cout << "\t\tGraph"                                                                    << std::endl;
    std::cout << "=========================================="                                   << std::endl;
    std::cout << "IO time:\t\t\t"           << total_timer.elapsed()                            << std::endl;
    std::cout << "Filename:\t\t\t"          << log_config.graph_filename                        << std::endl;
    std::cout << "|-Nodes:\t\t\t"           << number_of_nodes                                  << std::endl;
    std::cout << "|-Edges:\t\t\t"           << number_of_edges                                  << std::endl;
    std::cout << std::endl;
}

 void struction_log::print_config() {
    filebuffer_string << "\t\tConfiguration"                                                    << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Time limit:\t\t\t"         << log_config.time_limit                   << std::endl; 
    filebuffer_string << "Seed:\t\t\t\t"             << log_config.seed                         << std::endl; 
    filebuffer_string << "Reduction Style:\t\t"      << log_config.reduction_style_name         << std::endl; 
    filebuffer_string << "Reduction Config: \t\t"    << log_config.reduction_config_name        << std::endl;
    filebuffer_string                                                                           << std::endl;
    
    std::cout << "\t\tConfiguration"                                                            << std::endl;
    std::cout << "==========================================="                                  << std::endl;
    std::cout << "Time limit:\t\t\t"         << log_config.time_limit                           << std::endl; 
    std::cout << "Seed:\t\t\t\t"             << log_config.seed                                 << std::endl; 
    std::cout << "Reduction Style:\t\t"      << log_config.reduction_style_name                 << std::endl;
    std::cout << "Reduction Config: \t\t"    << log_config.reduction_config_name                << std::endl;
    std::cout                                                                                   << std::endl;
} 

void struction_log::print_reduction(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, size_t max_component) {
    filebuffer_string << "\t\tReduction"                                                        << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"       << offset                                      << std::endl;
    filebuffer_string << "Time:\t\t\t\t"         << time                                        << std::endl;
    filebuffer_string << "Kernel nodes:\t\t\t"   << kernel_size_n                               << std::endl;
    filebuffer_string << "MaxComponent:\t\t\t"   << max_component <<  "\n"                      << std::endl;

    std::cout << "=========================================="                               << std::endl;
    std::cout << "\t\tReduction"                                                            << std::endl;
    std::cout << "=========================================="                               << std::endl;
    std::cout << "Offset:\t\t\t\t"           << offset                                      << std::endl;
    std::cout << "Time:\t\t\t\t"             << time                                        << std::endl;
    std::cout << "Kernel nodes:\t\t\t"       << kernel_size_n                               << std::endl;
    std::cout << "MaxComponent:\t\t\t"       << max_component <<  "\n"                      << std::endl;
}


void struction_log::print_full_reduction(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m, size_t component_count, size_t max_component) {
    filebuffer_string << "\t\tReduction"                                                        << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"       << offset                                      << std::endl;
    filebuffer_string << "Time:\t\t\t\t"         << time                                        << std::endl;
    filebuffer_string << "Kernel nodes:\t\t\t"   << kernel_size_n                               << std::endl;
    filebuffer_string << "Kernel edges:\t\t\t"   << kernel_size_m                               << std::endl;
    filebuffer_string << "Components:\t\t\t"     << component_count                             << std::endl;
    filebuffer_string << "MaxComponent:\t\t\t"   << max_component <<  "\n"                      << std::endl;

    std::cout << "=========================================="                               << std::endl;
    std::cout << "\t\tReduction"                                                            << std::endl;
    std::cout << "=========================================="                               << std::endl;
    std::cout << "Offset:\t\t\t\t"           << offset                                      << std::endl;
    std::cout << "Time:\t\t\t\t"             << time                                        << std::endl;
    std::cout << "Kernel nodes:\t\t\t"       << kernel_size_n                               << std::endl;
    std::cout << "Kernel edges:\t\t\t"       << kernel_size_m                               << std::endl;
    std::cout << "Components:\t\t\t"         << component_count                             << std::endl;
    std::cout << "MaxComponent:\t\t\t"       << max_component <<  "\n"                      << std::endl;
}

void struction_log::print_results(bool optimal) {
    filebuffer_string                                                                           << std::endl;
    filebuffer_string << "\t\tStatistics"                                                       << std::endl;
    filebuffer_string << "=========================================="                           << std::endl;
    filebuffer_string << "Final Weight:\t\t\t"  << optimum_size                                 << std::endl;
    filebuffer_string << "Total time:\t\t\t\t"  << time_taken_best                              << std::endl;
    filebuffer_string << "Is optimal:\t\t\t"    << optimal                                      << std::endl;
    filebuffer_string << "Total time taken:\t\t"<< total_timer.elapsed()                        << std::endl;
    std::cout                                                                                   << std::endl;
    std::cout << "\t\tBest"                                                                     << std::endl;
    std::cout << "=========================================="                                   << std::endl;
    std::cout << "Final Weight:\t\t\t"          << optimum_size                                 << std::endl;
    std::cout << "Time found:\t\t\t"            << time_taken_best                              << std::endl;
    std::cout << "Is optimal:\t\t\t"            << optimal                                      << std::endl;
    std::cout << "Total time taken:\t\t"        << total_timer.elapsed()                        << std::endl;
    std::cout                                                                                   << std::endl;
}

void struction_log::restart_total_timer() {
    total_timer.restart();
}

void struction_log::reset_best_size() {
    optimum_size = 0;
}

void struction_log::set_best_size(unsigned int size) {
    if (size > optimum_size) {
        optimum_size = size;
    }
}
void struction_log::set_best_time(double time) {
        time_taken_best = time;
}

void struction_log::set_best(unsigned int size, double time) {
    if (size > optimum_size) {
        optimum_size = size;
        time_taken_best = time;
    }
}
