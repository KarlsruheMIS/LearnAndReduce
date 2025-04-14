
#include "log.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

log::log() {
    number_of_nodes = 0;
    number_of_edges = 0;
    avg_degree = 0.0;
    density = 0.0;
    arc_scans = 0;
    total_time_taken = 0.0;
    time_taken_best = 0.0;
    optimum_size = 0;
}

log::~log() { }

void log::set_config(ReductionConfig & config) {
    log_config = config; 
}

void log::set_graph(graph_access & G) {
    number_of_nodes = G.number_of_nodes();
    number_of_edges = G.number_of_edges();
    avg_degree = (double) number_of_edges / number_of_nodes;
    density = (double) (2 * number_of_edges) / (number_of_nodes * (number_of_nodes - 1));
}

void log::write_log() {
    std::stringstream filename_stream;
    filename_stream << "./logs/log_"<<  log_config.graph_filename <<   
                       "_seed_" <<  log_config.seed;
    std::ofstream f(filename_stream.str());
    f << filebuffer_string.str();
    f.close();
}


void log::print_data_generation_title() {
    filebuffer_string << "=============================================="                       << std::endl;
    filebuffer_string << "\t Generate Reduction Data "                                          << std::endl;
    filebuffer_string << "=============================================="                       << std::endl;

    std::cout << "=============================================="                               << std::endl;
    std::cout << "\t Generate Reduction Data"                                                   << std::endl;
    std::cout << "=============================================="                               << std::endl;
}

void log::print_reduction_title() {
    filebuffer_string << "=============================================="                      << std::endl;
    filebuffer_string << "\t\t LearnAndReduce MWIS"                                            << std::endl;
    filebuffer_string << "=============================================="                      << std::endl;

    std::cout << "=============================================="                              << std::endl;
    std::cout << "\t\t LearnAndReduce MWIS"                                                    << std::endl;
    std::cout << "=============================================="                              << std::endl;
}

void log::print_graph() {
    filebuffer_string << "\t\tGraph"                                                            << std::endl;
    filebuffer_string << "=============================================="                       << std::endl;
    filebuffer_string << "IO time:\t\t\t\t"         << total_timer.elapsed()                    << std::endl;
    filebuffer_string << "Filename:\t\t\t\t"        << log_config.graph_filename                << std::endl;
    filebuffer_string << "|-Nodes:\t\t\t\t"         << number_of_nodes                          << std::endl;
    filebuffer_string << "|-Edges:\t\t\t\t"         << number_of_edges/2                        << std::endl;
    filebuffer_string << std::endl;

    std::cout << "\t\tGraph"                                                                    << std::endl;
    std::cout << "=============================================="                               << std::endl;
    std::cout << "IO time:\t\t\t"           << total_timer.elapsed()                            << std::endl;
    std::cout << "Filename:\t\t\t"          << log_config.graph_filename                        << std::endl;
    std::cout << "|-Nodes:\t\t\t"           << number_of_nodes                                  << std::endl;
    std::cout << "|-Edges:\t\t\t"           << number_of_edges/2                                << std::endl;
    std::cout << std::endl;
}

void log::print_one_line_kernel_data(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m) {
    filebuffer_string << log_config.graph_filename << ","
                      << log_config.struction_config_name << ","
                      << log_config.reduction_config_name << ","
                      << log_config.time_limit << ","
                      << log_config.seed << ","
                      << log_config.gnn_filter_name << ","
                      << number_of_nodes << ","
                      << number_of_edges/2 << ","
                      << kernel_size_n << ","
                      << kernel_size_m << ","
                      << offset << ","
                      << time 
                      << std::endl;

    std::cout         << log_config.graph_filename << ","
                      << log_config.struction_config_name << ","
                      << log_config.reduction_config_name << ","
                      << log_config.time_limit << ","
                      << log_config.seed << ","
                      << log_config.gnn_filter_name << ","
                      << number_of_nodes << ","
                      << number_of_edges/2 << ","
                      << kernel_size_n << ","
                      << kernel_size_m << ","
                      << offset << ","
                      << time 
                      << std::endl;
}

 void log::print_config() {
    filebuffer_string << "\t\tConfiguration"                                                    << std::endl;
    filebuffer_string << "============================================="                        << std::endl;
    filebuffer_string << "Time limit:\t\t\t"         << log_config.time_limit                   << std::endl; 
    filebuffer_string << "Seed:\t\t\t\t"             << log_config.seed                         << std::endl; 
    filebuffer_string << "Struction Config:\t\t"      << log_config.struction_config_name       << std::endl; 
    filebuffer_string << "Reduction Config: \t\t"    << log_config.reduction_config_name        << std::endl;
    filebuffer_string << "GNN Filter:\t\t\t"         << log_config.gnn_filter_name              << std::endl;
    filebuffer_string                                                                           << std::endl;
    
    std::cout << "\t\tConfiguration"                                                            << std::endl;
    std::cout << "=============================================="                               << std::endl;
    std::cout << "Time limit:\t\t\t"         << log_config.time_limit                           << std::endl; 
    std::cout << "Seed:\t\t\t\t"             << log_config.seed                                 << std::endl; 
    std::cout << "Struction Config:\t\t"     << log_config.struction_config_name                << std::endl;
    std::cout << "Reduction Config: \t\t"    << log_config.reduction_config_name                << std::endl;
    std::cout << "GNN Filter:\t\t\t"         << log_config.gnn_filter_name                      << std::endl;
    std::cout                                                                                   << std::endl;
} 

void log::print_reduction(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, size_t max_component) {
    filebuffer_string << "\t\tReduction"                                                        << std::endl;
    filebuffer_string << "=============================================="                       << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"       << offset                                      << std::endl;
    filebuffer_string << "Time:\t\t\t\t"         << time                                        << std::endl;
    filebuffer_string << "Kernel nodes:\t\t\t"   << kernel_size_n                               << std::endl;
    filebuffer_string << "MaxComponent:\t\t\t"   << max_component <<  "\n"                      << std::endl;

    std::cout << "=============================================="                           << std::endl;
    std::cout << "\t\tReduction"                                                            << std::endl;
    std::cout << "=============================================="                           << std::endl;
    std::cout << "Offset:\t\t\t\t"           << offset                                      << std::endl;
    std::cout << "Time:\t\t\t\t"             << time                                        << std::endl;
    std::cout << "Kernel nodes:\t\t\t"       << kernel_size_n                               << std::endl;
    std::cout << "MaxComponent:\t\t\t"       << max_component <<  "\n"                      << std::endl;
}


// void log::print_full_reduction(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m, size_t component_count, size_t max_component) {
void log::print_full_reduction(ReductionConfig & mis_config, double time, NodeWeight offset, NodeID kernel_size_n, EdgeID kernel_size_m) {
    filebuffer_string << "\t\tReduction"                                                        << std::endl;
    filebuffer_string << "=============================================="                       << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"       << offset                                      << std::endl;
    filebuffer_string << "Time:\t\t\t\t"         << time                                        << std::endl;
    filebuffer_string << "Kernel nodes:\t\t\t"   << kernel_size_n                               << std::endl;
    filebuffer_string << "Kernel edges:\t\t\t"   << kernel_size_m                               << std::endl;
    // filebuffer_string << "Components:\t\t\t"     << component_count                             << std::endl;
    // filebuffer_string << "MaxComponent:\t\t\t"   << max_component <<  "\n"                      << std::endl;

    std::cout << "=============================================="                           << std::endl;
    std::cout << "\t\tReduction"                                                            << std::endl;
    std::cout << "=============================================="                           << std::endl;
    std::cout << "Offset:\t\t\t\t"           << offset                                      << std::endl;
    std::cout << "Time:\t\t\t\t"             << time                                        << std::endl;
    std::cout << "Kernel nodes:\t\t\t"       << kernel_size_n                               << std::endl;
    std::cout << "Kernel edges:\t\t\t"       << kernel_size_m                               << std::endl;
    // std::cout << "Components:\t\t\t"         << component_count                             << std::endl;
    // std::cout << "MaxComponent:\t\t\t"       << max_component <<  "\n"                      << std::endl;
}

void log::print_results(bool optimal) {
    filebuffer_string                                                                           << std::endl;
    filebuffer_string << "\t\tStatistics"                                                       << std::endl;
    filebuffer_string << "=============================================="                       << std::endl;
    filebuffer_string << "Final Weight:\t\t\t"  << optimum_size                                 << std::endl;
    filebuffer_string << "Total time:\t\t\t\t"  << time_taken_best                              << std::endl;
    filebuffer_string << "Is optimal:\t\t\t"    << optimal                                      << std::endl;
    filebuffer_string << "Total time taken:\t\t"<< total_timer.elapsed()                        << std::endl;
    std::cout                                                                                   << std::endl;
    std::cout << "\t\tBest"                                                                     << std::endl;
    std::cout << "=============================================="                               << std::endl;
    std::cout << "Final Weight:\t\t\t"          << optimum_size                                 << std::endl;
    std::cout << "Time found:\t\t\t"            << time_taken_best                              << std::endl;
    std::cout << "Is optimal:\t\t\t"            << optimal                                      << std::endl;
    std::cout << "Total time taken:\t\t"        << total_timer.elapsed()                        << std::endl;
    std::cout                                                                                   << std::endl;
}

void log::restart_total_timer() {
    total_timer.restart();
}

void log::reset_best_size() {
    optimum_size = 0;
}

void log::set_best_size(unsigned int size) {
    if (size > optimum_size) {
        optimum_size = size;
    }
}
void log::set_best_time(double time) {
        time_taken_best = time;
}

void log::set_best(unsigned int size, double time) {
    if (size > optimum_size) {
        optimum_size = size;
        time_taken_best = time;
    }
}
