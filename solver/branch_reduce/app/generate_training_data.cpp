/******************************************************************************
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


bool write_reduction_data(std::vector<std::vector<bool>> &reduction_data, std::string filename, std::vector<std::string> &reduction_names)
{
    int count_data_per_graph = 0;
    for (size_t i = 0; i < reduction_data.size(); i++)
    {
        // check at least 10% of applications
        if (std::count(reduction_data[i].begin(), reduction_data[i].end(), true) < 0.10 * reduction_data[i].size())
        {
            continue;
        }
        count_data_per_graph++;
    }
    if (count_data_per_graph < 2)
        return false; // specify how much data per graph needed

    for (size_t i = 0; i < reduction_data.size(); i++)
    {
        // check at least 10% of applications
        if (std::count(reduction_data[i].begin(), reduction_data[i].end(), true) < 0.10 * reduction_data[i].size())
        {
            continue;
        }

        std::cout << "Data is suitable for reduction " << reduction_names[i] << std::endl;
        std::string reduction_filename = filename + reduction_names[i] + ".txt";
        std::ofstream file;
        file.open(reduction_filename);
        for (size_t j = 0; j < reduction_data[i].size(); j++)
        {
            file << reduction_data[i][j] << std::endl;
        }
        file.close();
    }
    return true;
}

bool write_reduction_data_csv(graph_access &G, std::vector<std::vector<bool>> &reduction_data,
                              std::string filename, std::vector<std::string> &reduction_names,
                              std::vector<bool> &exclude_data, std::vector<bool> &include_data)
{
    std::vector<bool> used_reduction(reduction_data.size(), 0);
    int count_data_per_graph = 0;
    for (size_t i = 0; i < reduction_data.size(); i++)
    {
        // check at least 10% of applications
        if (std::count(reduction_data[i].begin(), reduction_data[i].end(), true) < 0.10 * reduction_data[i].size())
        {
            continue;
        }
        count_data_per_graph++;
        used_reduction[i] = true;
    }
    if (count_data_per_graph < 2 || G.number_of_nodes() < 100)
        return false; // specify how much data per graph needed

    for (size_t i = 0; i < reduction_data.size(); i++)
    {
        if (used_reduction[i])
            std::cout << reduction_names[i] << std::endl;
    }

    float *node_attr = NULL, *edge_attr = NULL;
    // LRConv::compute_node_attr(&node_attr, G);
    LRConv::compute_attr(&node_attr, &edge_attr, G);

    std::ofstream file;
    file.open(filename + ".csv");
    file << "source;target;uc;vc;ic;uw;vw;iw;twin;dom" << std::endl;
    for (NodeID u = 0; u < G.number_of_nodes(); u++)
    {
        for (NodeID i = G.get_first_edge(u); i != G.get_first_invalid_edge(u); i++)
        {
            float *ed = edge_attr + (edge_features * i);
            file << u << ";" << G.getEdgeTarget(i);
            for (NodeID j = 0; j < edge_features; j++)
                file << ";" << ed[j];
            file << std::endl;
        }
    }
    file.close();

    file.open(filename + "_meta.csv");
    file << "id;d;w;nw;l;i;e";
    for (NodeID i = 0; i < reduction_data.size(); i++)
        if (used_reduction[i])
            file << ";" << reduction_names[i];
    file << std::endl;

    for (NodeID u = 0; u < G.number_of_nodes(); u++)
    {
        file << u;
        float *nd = node_attr + (u * node_features);
        for (NodeID i = 0; i < node_features; i++)
            file << ";" << nd[i];
        file << ";" << include_data[u] << ";" << exclude_data[u];
        for (NodeID i = 0; i < reduction_data.size(); i++)
            if (used_reduction[i])
                file << ";" << reduction_data[i][u];
        file << std::endl;
    }
    file.close();
    free(node_attr);
    free(edge_attr);
    return true;
}

int main(int argn, char **argv)
{
    struction_log::instance()->restart_total_timer();
    /* bool output_convergence = false; */
    /* bool output_best_solution = true; */
    struction_log::instance()->print_title();
    ReductionConfig config;

    // Parse the command line parameters;
    ReductionArguments arguments(argn, argv);
    int ret_code = arguments.setConfig(config);
    if (ret_code)
        return 0;

    // config.generate_training_data = true;
    config.disable_critical_set = true;
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
    struction_log::instance()->set_graph(G);
    struction_log::instance()->print_graph();

    branch_and_reduce_algorithm reducer(G, config);

    if (config.size_of_subgraph > G.number_of_nodes())
    {
        config.size_of_subgraph = G.number_of_nodes();
        config.num_of_subgraphs = 1;
    }

    std::vector<std::string> transformation_names;
    reducer.get_transformation_names(transformation_names);
    size_t num_of_reductions = transformation_names.size();

    std::vector<std::vector<bool>> reduction_data(num_of_reductions, std::vector<bool>(config.size_of_subgraph, false));
    std::vector<bool> exclude_data(config.size_of_subgraph, false);
    std::vector<bool> include_data(config.size_of_subgraph, false);
    graph_access subgraph;
    config.seed = time(NULL);
    for (int i = 0; i < config.num_of_subgraphs; i++)
    {
        // number of different subgraphs out of one graph
        reducer.get_training_data_for_graph_size(subgraph, config.size_of_subgraph, reduction_data, include_data, exclude_data, i); // size of the subgraph

        write_reduction_data_csv(subgraph, reduction_data, "training_data/csv/" + name + "_seed" + std::to_string(config.seed) + "_training_data_graph_" + std::to_string(i), transformation_names, exclude_data, include_data);

        // if (write_reduction_data(reduction_data, "training_data/reduction_data/" + name + "_seed" + std::to_string(config.seed) + "_training_data_graph_" + std::to_string(i) + "_", reduction_names))
        //     graph_io::writeGraphWeighted(subgraph, "training_data/graph/" + name + "_seed" + std::to_string(config.seed) + "_training_data_graph_" + std::to_string(i) + ".graph");

        reduction_data = std::vector<std::vector<bool>>(num_of_reductions, std::vector<bool>(config.size_of_subgraph, false));
    }

    return 0;
}
