#pragma once

#include <string>

#include "graph_access.h"
#include "branch_and_reduce_algorithm.h"

const int node_features = 4;
const int edge_features = 8;
const int hidden_dim = 16;
const float scale = 100.0f;
const double id_scale = 1000000.0;

class LRConv
{
public:
    LRConv(const std::string path);
    ~LRConv();

    void change_parameters(const std::string path);

    const float *predict(const float *node_attr, const float *edge_attr, graph_access &g);

    const float *predict_light(graph_access &g);

    const float *predict_light_blas(graph_access &g);

    const float *predict_light_dynamic_blas(branch_and_reduce_algorithm *g);

    static void compute_attr(float **node_attr, float **edge_attr, graph_access &g);

    static void compute_node_attr(float **node_attr, graph_access &g);

    static void compute_node_attr_dynamic(float **node_attr, branch_and_reduce_algorithm *g);

private:
    int params;
    float *param;

    int allocated;
    float *x, *y;
};