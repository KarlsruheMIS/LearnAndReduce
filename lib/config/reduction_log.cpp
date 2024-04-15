/******************************************************************************
 * reduction_log.cpp
 *****************************************************************************/

#include "reduction_log.h"

#include <stdio.h>
#include <fstream>
#include <iostream>

reduction_log::reduction_log() : mis_log() { 
    kernel_size = number_of_nodes;
    offset = 0;
}

reduction_log::~reduction_log() {}

void reduction_log::print_config() {
    filebuffer_string << "\t\tConfiguration"        << std::endl;
    filebuffer_string << "=========================================="                            << std::endl;
    filebuffer_string << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    filebuffer_string << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    filebuffer_string << "Reduction Style:\t\t"           << log_config.reduction_style_name     << std::endl; 
    filebuffer_string << std::endl;
    
    std::cout << "\t\tConfiguration"        << std::endl;
    std::cout << "==========================================="                           << std::endl;
    std::cout << "Time limit:\t\t\t"         << log_config.time_limit                    << std::endl; 
    std::cout << "Seed:\t\t\t\t"             << log_config.seed                          << std::endl; 
    std::cout << "Reduction Style:\t\t"      << log_config.reduction_style_name          << std::endl;
    std::cout << std::endl;
}


void reduction_log::print_title() {
    filebuffer_string << "==========================================="                   << std::endl;
    filebuffer_string << "\tMaximum Weight Independent Set Reduction"                  << std::endl;
    filebuffer_string << "==========================================="                   << std::endl;

    std::cout << "==========================================="                           << std::endl;
    std::cout << "\tMaximum Weight Independent Set Reduction"                      << std::endl;
    std::cout << "==========================================="                           << std::endl;
}

void reduction_log::print_reduction(ReductionConfig & mis_config, NodeWeight offset, NodeID kernel_size) {
    filebuffer_string << "\t\tReduction"                                                         << std::endl;
    filebuffer_string << "==========================================="                           << std::endl;
    filebuffer_string << "Offset:\t\t\t\t"            << offset                                  << std::endl;
    filebuffer_string << "Kernel size:\t\t\t"         << kernel_size                             << std::endl;
    filebuffer_string << "Time:\t\t\t\t"              << total_timer.elapsed()           << "\n" << std::endl;

    std::cout << "==========================================="                               << std::endl;
    std::cout << "\t\tReduction"                                                             << std::endl;
    std::cout << "==========================================="                               << std::endl;
    std::cout << "Offset:\t\t\t\t"                << offset                                  << std::endl;
    std::cout << "Kernel size:\t\t\t"             << kernel_size                             << std::endl;
    std::cout << "Time:\t\t\t\t"                  << total_timer.elapsed()           << "\n" << std::endl;
}
