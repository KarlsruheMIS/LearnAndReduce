#pragma once

#include <string>

#include "graph_access.h"

const int node_features = 4;
const int edge_features = 8;
const int hidden_dim = 16;
const float scale = 100.0f;

class LRConv
{
public:
    LRConv(const std::string path);
    ~LRConv();

    void change_parameters(const std::string path);

    const float *predict(const float *node_attr, const float *edge_attr, graph_access &g);

    const float *predict_light(graph_access &g);

    static void compute_attr(float **node_attr, float **edge_attr, graph_access &g);

private:
    int params;
    float *param;

    int allocated;
    float *x, *y;
};