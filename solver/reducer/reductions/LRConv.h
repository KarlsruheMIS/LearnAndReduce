#pragma once

#include <string>
#include <unordered_map>

#include "graph_access.h"

// const int node_features = 4;
// const int edge_features = 8;
const int hidden_dim = 16;
// const float scale = 100.0f;
const double id_scale = 1000000.0;

const int total_node_features = 16;
const int total_edge_features = 8;

class reduce_algorithm;

class LRConv
{
public:
    LRConv(int N);
    // LRConv(const LRConv &gnn);
    ~LRConv();

    void change_parameters(const std::string path);

    const float *predict(reduce_algorithm *r_alg);

    // static void compute_attr(float **node_attr, float **edge_attr, graph_access &g);

    static void compute_attr_norm(float **node_attr, float **edge_attr, graph_access &g);

    // void compute_node_attr_dynamic(branch_and_reduce_algorithm *g);

    // void compute_node_attr_dynamic_norm(branch_and_reduce_algorithm *g);

private:
    std::unordered_map<std::string, std::vector<float>> params;

    float **W;
    float **B;
    float *param;

    float *x, *y, *e;
};