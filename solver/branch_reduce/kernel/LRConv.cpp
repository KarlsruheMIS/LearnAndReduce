#include "LRConv.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

LRConv::LRConv(const std::string path)
    : params(0), param(NULL), allocated(0), x(NULL), y(NULL)
{
    change_parameters(path);
}

LRConv::~LRConv()
{
    free(param);
    free(x);
    free(y);
}

void LRConv::change_parameters(const std::string path)
{
    free(param);

    std::ifstream file;
    file.open(path);

    file >> params;
    param = (float *)aligned_alloc(32, sizeof(float) * params);

    for (int i = 0; i < params; i++)
        file >> param[i];

    file.close();
}

void LRConv_message(const float *x, int dim_in,
                    float *y, int dim_out,
                    const float *param1, const float *bias1,
                    const float *param2, const float *bias2,
                    graph_access &g)
{
    float *buff = (float *)aligned_alloc(32, sizeof(float) * dim_in * 2);
    float *tmp1 = (float *)aligned_alloc(32, sizeof(float) * dim_out);
    float *tmp2 = (float *)aligned_alloc(32, sizeof(float) * dim_out);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        float *agg = y + dim_out * u;
        for (int i = 0; i < dim_out; i++)
            agg[i] = 0.0f;

        for (int i = 0; i < dim_in; i++)
            buff[i] = x[u * dim_in + i];

        for (int e = g.get_first_edge(u); e != g.get_first_invalid_edge(u); e++)
        {
            int v = g.getEdgeTarget(e);
            for (int i = 0; i < dim_in; i++)
                buff[dim_in + i] = x[v * dim_in + i];

            for (int i = 0; i < dim_out; i++)
            {
                tmp1[i] = bias1[i];
                tmp2[i] = bias2[i];
            }

            // Layer 1
            for (int i = 0; i < dim_out; i++)
                for (int j = 0; j < dim_in * 2; j++)
                    tmp1[i] += buff[j] * param1[i * (dim_in * 2) + j];

            for (int i = 0; i < dim_out; i++)
                tmp1[i] = tmp1[i] >= 0.0f ? tmp1[i] : 0.0f;

            // Layer 2
            for (int i = 0; i < dim_out; i++)
                for (int j = 0; j < dim_out; j++)
                    tmp2[i] += tmp1[j] * param2[i * dim_out + j];

            for (int i = 0; i < dim_out; i++)
                tmp2[i] = tmp2[i] >= 0.0f ? tmp2[i] : 0.0f;

            // Aggregate
            for (int i = 0; i < dim_out; i++)
                agg[i] = agg[i] < tmp2[i] ? tmp2[i] : agg[i];
        }
    }

    free(buff);
    free(tmp1);
    free(tmp2);
}

const float *LRConv::predict_light(graph_access &g)
{
    if (allocated < g.number_of_nodes() * hidden_dim)
    {
        allocated = g.number_of_nodes() * hidden_dim;
        free(x);
        free(y);
        x = (float *)aligned_alloc(32, sizeof(float) * allocated);
        y = (float *)aligned_alloc(32, sizeof(float) * allocated);
    }

    // Compute vertex features
    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        int Wn = 0;
        for (int e = g.get_first_edge(u); e != g.get_first_invalid_edge(u); e++)
        {
            int v = g.getEdgeTarget(e);
            Wn += g.getNodeWeight(v);
        }

        x[u * node_features + 0] = (float)g.getNodeDegree(u) / scale;
        x[u * node_features + 1] = (float)g.getNodeWeight(u) / scale;
        x[u * node_features + 2] = (float)Wn / scale;
        x[u * node_features + 3] = (float)u / scale;
    }

    int p_input = (node_features * 2) * hidden_dim;
    int b = hidden_dim;
    int p_hidden_large = (hidden_dim * 2) * hidden_dim;
    int p_hidden = hidden_dim * hidden_dim;

    float *p1 = param;
    float *p2 = param + p_input + b;
    LRConv_message(x, node_features, y, hidden_dim,
                   p1, p1 + p_input, p2, p2 + p_hidden, g);

    float *p3 = p2 + p_hidden + b;
    float *p4 = p3 + p_hidden_large + b;
    LRConv_message(y, hidden_dim, x, hidden_dim,
                   p3, p3 + p_hidden_large, p4, p4 + p_hidden, g);

    float *p5 = p4 + p_hidden + b;
    float *p6 = p5 + p_hidden + b;
    float *tmp = (float *)aligned_alloc(32, sizeof(float) * hidden_dim);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        y[u] = *(p6 + hidden_dim);

        for (int i = 0; i < hidden_dim; i++)
            tmp[i] = p5[hidden_dim + i];

        // Layer 1
        for (int i = 0; i < hidden_dim; i++)
            for (int j = 0; j < hidden_dim; j++)
                tmp[i] += x[u * hidden_dim + j] * p5[i * hidden_dim + j];

        for (int i = 0; i < hidden_dim; i++)
            tmp[i] = tmp[i] >= 0.0f ? tmp[i] : 0.0f;

        // Layer 2
        for (int i = 0; i < hidden_dim; i++)
            y[u] += tmp[i] * p6[i];
    }

    free(tmp);

    return y;
}

const float *LRConv::predict_light_blas(graph_access &g)
{
    
}

const float *LRConv::predict(const float *node_attr, const float *edge_attr, graph_access &g)
{
    // iterate over edges
    // kernel matmul (bias)
    // max aggr
    // transform twice
    // output

    int dim = node_features * 2 + edge_features;
    float *buff = (float *)aligned_alloc(32, sizeof(float) * dim);
    float *tmp1 = (float *)aligned_alloc(32, sizeof(float) * dim);
    float *tmp2 = (float *)aligned_alloc(32, sizeof(float) * dim);
    float *agg = (float *)aligned_alloc(32, sizeof(float) * dim);
    if (allocated < g.number_of_nodes())
    {
        allocated = g.number_of_nodes();
        free(y);
        y = (float *)aligned_alloc(32, sizeof(float) * g.number_of_nodes());
    }

    float *w1 = param;
    float *b1 = w1 + dim * dim;
    float *w2 = b1 + dim;
    float *b2 = w2 + dim * dim;
    float *w3 = b2 + dim;
    float *b3 = w3 + dim * dim;
    float *w4 = b3 + dim;
    float *b4 = w4 + dim;

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        for (int i = 0; i < dim; i++)
            agg[i] = 0.0f;

        for (int i = 0; i < node_features; i++)
            buff[i] = node_attr[u * node_features + i];

        for (int e = g.get_first_edge(u); e != g.get_first_invalid_edge(u); e++)
        {
            int v = g.getEdgeTarget(e);
            for (int i = 0; i < node_features; i++)
                buff[node_features + i] = node_attr[v * node_features + i];

            for (int i = 0; i < edge_features; i++)
                buff[2 * node_features + i] = edge_attr[e * edge_features + i];

            for (int i = 0; i < dim; i++)
            {
                tmp1[i] = b1[i];
                tmp2[i] = b2[i];
            }

            // Layer 1
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    tmp1[i] += buff[j] * w1[i * dim + j];

            for (int i = 0; i < dim; i++)
                tmp1[i] = tmp1[i] >= 0.0f ? tmp1[i] : 0.0f;

            // Layer 2
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    tmp2[i] += tmp1[j] * w2[i * dim + j];

            for (int i = 0; i < dim; i++)
                tmp2[i] = tmp2[i] >= 0.0f ? tmp2[i] : 0.0f;

            // Aggregate
            for (int i = 0; i < dim; i++)
                agg[i] = agg[i] < tmp2[i] ? tmp2[i] : agg[i];
        }

        for (int i = 0; i < dim; i++)
            tmp1[i] = b3[i];

        y[u] = *b4;

        // Layer 3
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                tmp1[i] += agg[j] * w3[i * dim + j];

        for (int i = 0; i < dim; i++)
            tmp1[i] = tmp1[i] >= 0.0f ? tmp1[i] : 0.0f;

        // Layer 4
        for (int i = 0; i < dim; i++)
            y[u] += tmp1[i] * w4[i];
    }

    free(buff);
    free(tmp1);
    free(tmp2);
    free(agg);

    return y;
}

void LRConv::compute_attr(float **node_attr, float **edge_attr, graph_access &g)
{
    free(*node_attr);
    *node_attr = (float *)aligned_alloc(32, sizeof(float) * g.number_of_nodes() * node_features);
    free(*edge_attr);
    *edge_attr = (float *)aligned_alloc(32, sizeof(float) * g.number_of_edges() * edge_features);

    std::vector<int> marks(g.number_of_nodes(), -1);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        int d = g.getNodeDegree(u), Wn = 0, C = 0;

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);
            marks[v] = u;
            Wn += g.getNodeWeight(v);
        }

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);

            int uc = d, vc = 0, ic = 0;
            int uw = Wn, vw = 0, iw = 0;

            for (int j = g.get_first_edge(v); j != g.get_first_invalid_edge(v); j++)
            {
                int w = g.getEdgeTarget(j);
                if (marks[w] == u)
                {
                    uc--;
                    uw -= g.getNodeWeight(w);
                    ic++;
                    iw += g.getNodeWeight(w);

                    C++;
                }
                else
                {
                    vc++;
                    vw += g.getNodeWeight(w);
                }
            }
            uc--;
            uw -= g.getNodeWeight(v);
            vc--;
            vw -= g.getNodeWeight(u);

            float *ea = *edge_attr + (i * edge_features);
            ea[0] = (float)uc / scale;
            ea[1] = (float)vc / scale;
            ea[2] = (float)ic / scale;
            ea[3] = (float)uw / scale;
            ea[4] = (float)vw / scale;
            ea[5] = (float)iw / scale;
            ea[6] = (uc == 0 && vc == 0) ? 1.0f : 0.0f;
            ea[7] = (uc == 0 && vc > 0) ? 1.0f : 0.0f;
        }

        float *ua = *node_attr + (u * node_features);
        ua[0] = (float)d / scale;
        ua[1] = (float)g.getNodeWeight(u) / scale;
        ua[2] = (float)Wn / scale;
        ua[3] = (float)u / scale;
        // ua[3] = (d < 2) ? 1.0f : (float)C / (float)((d * d) - d);
    }
}

unsigned int hash(unsigned int x)
{
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

void LRConv::compute_node_attr(float **node_attr, graph_access &g)
{
    free(*node_attr);
    *node_attr = (float *)aligned_alloc(32, sizeof(float) * g.number_of_nodes() * node_features);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        int d = g.getNodeDegree(u), Wn = 0;

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);
            Wn += g.getNodeWeight(v);
        }

        float *ua = *node_attr + (u * node_features);
        ua[0] = (float)d / scale;
        ua[1] = (float)g.getNodeWeight(u) / scale;
        ua[2] = (float)Wn / scale;
        ua[3] = (float)(hash(u) % 1000000) / id_scale;
    }
}