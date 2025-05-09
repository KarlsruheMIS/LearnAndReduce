
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
#include "log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "reduction_config.h"
#include "parse_parameters.h"
#include "reduce_algorithm.h"
#include "graph_operations.h"
#include "solution_check.h"
#include "LRConv.h"
#include "configuration_reduction.h"

bool write_reduction_data_csv(graph_access &G, std::vector<std::vector<int>> &reduction_data,
                              std::string filename, std::vector<std::string> &reduction_names)
{
    std::ofstream file;

    // write header
    file.open(filename + "_reduction_data.csv");
    file << "id;w";
    for (NodeID i = 0; i < reduction_data.size(); i++)
        file << ";" << reduction_names[i];
    file << std::endl;
    file << std::fixed << std::setprecision(6);

    // write data
    for (NodeID u = 0; u < G.number_of_nodes(); u++)
    {
        file << u << ";" << G.getNodeWeight(u);
        for (NodeID i = 0; i < reduction_data.size(); i++)
        {
            file << ";" << reduction_data[i][u];
        }
        file << std::endl;
    }

    file.close();
    return true;
}

void generate_data(reduce_algorithm *reducer, graph_access &G, ReductionConfig &config, std::string name)
{
    std::vector<std::string> transformation_names;
    reducer->get_transformation_names(transformation_names);
    size_t num_of_reductions = transformation_names.size();
    std::vector<std::vector<int>> reduction_data(num_of_reductions, std::vector<int>(G.number_of_nodes(), 0));

    reducer->generate_initial_reduce_data(G, reduction_data); 
    write_reduction_data_csv(G, reduction_data, "training_data/csv/" + name , transformation_names);
}

int main(int argn, char **argv)
{
    log::instance()->restart_total_timer();
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
    log::instance()->set_config(config);
    if (config.verbose)
    {
        log::instance()->print_data_generation_title();
        log::instance()->print_config();
    }

    // Read the graph
    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    if (G.number_of_nodes() < 1000)
    {
        if (config.verbose)
            std::cout << "Graph " << name << " is too small for training data generation" << std::endl;
        return 0;
    }

    go.assign_weights(G, config);
    go.writeGraphWeighted_to_csv(G, "training_data/csv/" + name + "_original_graph.csv");

    log::instance()->set_graph(G);
    if (config.verbose)
        log::instance()->print_graph();

    reduce_algorithm reducer(G, config);
    generate_data(&reducer, G, config, name+"_original");


    // disable reductions for kernel computation
    config.disable_unconfined                 = true;
    config.disable_extended_single_edge                = true;
    config.disable_single_edge                   = true;
    config.disable_funnel                     = true;
    config.disable_critical_set               = true;


    reduce_algorithm r(G, config);
    graph_access& kernel = r.kernelize();
    if (kernel.number_of_nodes() == 0 ) return 0;
    go.writeGraphWeighted_to_csv(kernel, "training_data/csv/" + name + "_kernel_graph.csv");

    // enable all reductions
    config.disable_unconfined                 = false;
    config.disable_extended_single_edge                = false;
    config.disable_single_edge                   = false;
    config.disable_funnel                     = false;
    config.disable_generalized_fold           = false;
    config.disable_cut_vertex                 = false;
    config.disable_heavy_set                  = false;
    config.disable_heavy_set3                 = false;
    config.disable_critical_set               = false;

    reduce_algorithm reducer2(kernel, config);
    generate_data(&reducer2, kernel, config, name+"_kernel");
    return 0;
}
