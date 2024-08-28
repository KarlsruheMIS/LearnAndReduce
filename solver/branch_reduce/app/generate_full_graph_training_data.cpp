
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <argtable3.h>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>
#include <iomanip>

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
#include "LRConv.h"
#include "configuration_reduction.h"

bool write_reduction_data_csv(graph_access &G, std::vector<std::vector<bool>> &reduction_data,
                              std::string filename, std::vector<std::string> &reduction_names)
{
    std::ofstream file;

    // write header
    file.open(filename + "_reduction_data.csv");
    file << "id";
    for (NodeID i = 0; i < reduction_data.size(); i++)
        file << ";" << reduction_names[i];
    file << std::endl;
    file << std::fixed << std::setprecision(6);

    // write data
    for (NodeID u = 0; u < G.number_of_nodes(); u++)
    {
        file << u ;
        for (NodeID i = 0; i < reduction_data.size(); i++)
        {
            file << ";" << reduction_data[i][u];
        }
        file << std::endl;
    }

    file.close();
    return true;
}

int main(int argn, char **argv)
{
    struction_log::instance()->restart_total_timer();
    struction_log::instance()->print_title();
    ReductionConfig config;

    // Parse the command line parameters;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code)
        return 0;

    configuration_reduction cfg;
    cfg.generate_training_data_initial_reductions(config);

    std::string graph_filepath = config.graph_filename;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string name = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    struction_log::instance()->set_config(config);
    struction_log::instance()->print_config();

    // Read the graph
    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);

    go.assign_weights(G, config);
    go.writeGraphWeighted_to_csv(G, "training_data/csv/" + name + "_original_graph.csv");
    go.writeWeights_to_csv(G, "training_data/csv/" + name + "_original_weights.csv");

    struction_log::instance()->set_graph(G);
    struction_log::instance()->print_graph();

    //initially reduce graph with simple reductions
    std::vector<std::string> transformation_names;
    branch_and_reduce_algorithm reducer(G, config);
    reducer.get_transformation_names(transformation_names);
    size_t num_of_reductions = transformation_names.size();
    std::vector<std::vector<bool>> reduction_data(num_of_reductions, std::vector<bool>(G.number_of_nodes(), false));

    reducer.generate_initial_reduce_data(G, reduction_data); 
    write_reduction_data_csv(G, reduction_data, "training_data/csv/" + name + "_initial" , transformation_names);

    // build and write the kernel
    graph_access &kernel = reducer.kernelize();
    go.writeGraphWeighted_to_csv(kernel, "training_data/csv/" + name + "_kernel_graph.csv");
    go.writeWeights_to_csv(kernel, "training_data/csv/" + name + "_kernel_weights.csv");

    // config for full reductions
    if (kernel.number_of_nodes() == 0)
    {
        std::cout << "Kernel is empty, no further reductions possible." << std::endl;
        return 0;
    }

    // enable expensive reductions
    config.disable_generalized_fold           = false;
    config.disable_clique_neighborhood        = false;
    config.disable_clique_neighborhood_fast   = false;
    config.disable_cut_vertex                 = false;
    config.disable_heavy_set                  = false;
    config.disable_heavy_set3                 = false;
    config.disable_critical_set               = false;
  
    std::vector<std::string> transformation_names_full;
    branch_and_reduce_algorithm kernel_reducer(kernel, config);
    kernel_reducer.get_transformation_names(transformation_names_full);
    size_t num_of_reductions_full = transformation_names_full.size();
    std::vector<std::vector<bool>> reduction_data_full(num_of_reductions_full, std::vector<bool>(kernel.number_of_nodes(), false));

    kernel_reducer.generate_initial_reduce_data(kernel, reduction_data_full); 
    write_reduction_data_csv(kernel, reduction_data_full, "training_data/csv/" + name + "_kernel" , transformation_names_full);

    return 0;
}
