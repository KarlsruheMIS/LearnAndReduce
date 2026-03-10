
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
#include "log.h"
#include "graph_access.h"
#include "graph_io.h"
#include "cout_handler.h"
#include "reduction_config.h"
#include "parse_parameters.h"
#include "reduce_algorithm.h"
#include "solution_check.h"
#include "graph_operations.h"

int main(int argn, char **argv)
{
    log::instance()->restart_total_timer();

    ReductionConfig config;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code)
        return 0;

    // read and set the graph
    std::string graph_filepath = config.graph_filename;
    std::string solution_filepath;
    config.graph_filename = graph_filepath.substr(graph_filepath.find_last_of('/') + 1);
    std::string path = graph_filepath.substr(0, graph_filepath.find_last_of('/'));
    std::string name = config.graph_filename.substr(0, config.graph_filename.find_last_of('.'));
    config.graph_filename = name;
    std::string path_and_file = graph_filepath.substr(0, graph_filepath.find_last_of('-'));
    log::instance()->set_config(config);
    config.solution_from_file = true;
    if (config.kernel_filename.empty())
    {
        std::cout <<  "Please provide a reduced instance file name using --kernel=<file_name> and run again." << std::endl;
        exit(-1);
    }
        

    graph_access G;
    graph_operations go;
    graph_io::readGraphWeighted(G, graph_filepath);
    go.assign_weights(G, config);
    log::instance()->set_graph(G);

    if (config.verbose)
    {
        log::instance()->print_reduction_title();
        log::instance()->print_config();
        log::instance()->print_graph();
    }
    else
    {
        config.print_reduction_info = false;
    }

    reduce_algorithm reducer(G, config);
    timer t;

    graph_access &g = reducer.kernelize();
    double time = t.elapsed();
    if (config.verbose)
    {
        log::instance()->print_full_reduction(config, time, reducer.get_current_is_weight(), g.number_of_nodes(), g.number_of_edges() / 2);
    }

    if (g.number_of_nodes() > 0) // always write reduced instance
    {
        go.writeGraphWeighted(g, config.kernel_filename + ".kernel_graph");
        // go.writeGraphWeighted_to_csv(g, config.kernel_filename + ".csv");
        // go.writeWeights_to_csv(g, config.kernel_filename + "_weight.csv");
    }

    std::vector<bool> reduced_solution(g.number_of_nodes(), false);
    bool computed_reduced_solution = false;
    if (g.number_of_nodes() > 0)
    {
        std::cout << "Please give the path to a solution for the reduced graph "<< config.kernel_filename <<".kernel_graph... " << std::endl;
        std::cin >> solution_filepath;

        computed_reduced_solution = graph_io::readVector<bool>(reduced_solution, solution_filepath);
    }

    if (g.number_of_nodes() == 0 || computed_reduced_solution)
    {
        std::vector<bool> full_solution(G.number_of_nodes(), false);
        NodeWeight solution_weight = reducer.lift_solution(reduced_solution, full_solution);

        if (config.write_solution)
            graph_io::writeVector(full_solution, config.output_filename);

        if (config.verbose)
            std::cout << "final solution weight is \t\t" << solution_weight << std::endl;

        log::instance()->print_one_line_solution_data(config, time, reducer.get_current_is_weight(), g.number_of_nodes(), g.number_of_edges() / 2,  solution_weight, time);
    }
    return 0;
}
