#include "lr_gcn.h"

#include <stdlib.h>
#include <immintrin.h>
#include <math.h>
#include <string.h>

#include <vector>

const int param_offsets[4] = {0, 272, 800, 1456};
const int layer_inputs[4] = {8, 16, 40, 16};
const int layer_outputs[4] = {16, 16, 16, 1};

const int max_input = 16;
const int cat_input = 40;

const int layer_output = 16;

const int num_features = 8;

const int node_weight = 0;
const int diff_neighborhood_weight = 1;
const int min_neighborhood_weight = 2;
const int max_neighborhood_weight = 3;
const int node_degree = 4;
const int avg_neighborhood_degree = 5;
const int min_neighborhood_degree = 6;
const int max_neighborhood_degree = 7;

int lr_gcn_parse(FILE *f, float **paramv)
{
    int paramc;
    int np = fscanf(f, "%d\n", &paramc);

    *paramv = (float *)malloc(sizeof(float) * paramc);
    for (int i = 0; i < paramc; i++)
        np = fscanf(f, "%f\n", (*paramv) + i);

    return paramc;
}

void compute_features(graph_access &g, float *x)
{
    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        float *ux = x + u * num_features;

        for (int i = 0; i < num_features; i++)
            ux[i] = 0.0;

        ux[node_weight] = g.getNodeWeight(u);
        ux[diff_neighborhood_weight] = g.getNodeWeight(u);
        ux[node_degree] = g.getNodeDegree(u);

        float scale = 1.0f / ux[node_degree];

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);

            float d = g.getNodeDegree(v);
            float w = g.getNodeWeight(v);

            ux[diff_neighborhood_weight] -= w;
            if (w < ux[min_neighborhood_weight] || ux[min_neighborhood_weight] == 0.0f)
                ux[min_neighborhood_weight] = w;
            if (w > ux[max_neighborhood_weight])
                ux[max_neighborhood_weight] = w;

            ux[avg_neighborhood_degree] += d * scale;
            if (d < ux[min_neighborhood_degree] || ux[min_neighborhood_degree] == 0.0f)
                ux[min_neighborhood_degree] = d;
            if (d > ux[max_neighborhood_degree])
                ux[max_neighborhood_degree] = d;
        }
    }
}

void matrix_multiplication_avx(const float *A, const float *B, float *C, const float *b, int n, int k)
{
    const int m = layer_output;
    for (int i1 = 0; i1 < n; i1++)
    {
        __m256 c1 = _mm256_load_ps(b);
        __m256 c2 = _mm256_load_ps(b + 8);

        for (int i3 = 0; i3 < k; i3++)
        {
            __m256 a = _mm256_broadcast_ss(A + i1 * k + i3);

            __m256 b1 = _mm256_load_ps(B + i3 * m);
            __m256 b2 = _mm256_load_ps(B + i3 * m + 8);

            c1 = _mm256_fmadd_ps(a, b1, c1);
            c2 = _mm256_fmadd_ps(a, b2, c2);
        }

        _mm256_store_ps(C + i1 * m, c1);
        _mm256_store_ps(C + i1 * m + 8, c2);
    }
}

void lr_conv_combined_avx(graph_access &g, int dim, float *in, float *out, const float *p, const float *b)
{
    float *in_buff = (float *)malloc(sizeof(float) * dim * 2);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        float *x = in + u * dim;
        float *y = out + u * layer_output;

        __m256 y1 = _mm256_setzero_ps();
        __m256 y2 = _mm256_setzero_ps();

        for (int j = 0; j < layer_output; j++)
            y[j] = 0.0f;

        for (int j = 0; j < dim; j++)
            in_buff[j] = x[j];

        float degree = g.getNodeDegree(u);
        __m256 scale = _mm256_set1_ps(1.0f / degree);

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);
            float *xv = in + v * dim;

            for (int j = 0; j < dim; j++)
                in_buff[dim + j] = xv[j];

            __m256 out1 = _mm256_load_ps(b);
            __m256 out2 = _mm256_load_ps(b + 8);

            for (int j = 0; j < dim * 2; j++)
            {
                __m256 a = _mm256_broadcast_ss(in_buff + j);

                __m256 b1 = _mm256_load_ps(p + j * layer_output);
                __m256 b2 = _mm256_load_ps(p + j * layer_output + 8);

                out1 = _mm256_fmadd_ps(a, b1, out1);
                out2 = _mm256_fmadd_ps(a, b2, out2);
            }

            y1 = _mm256_fmadd_ps(out1, scale, y1);
            y2 = _mm256_fmadd_ps(out2, scale, y2);
        }

        _mm256_store_ps(y, y1);
        _mm256_store_ps(y + 8, y2);
    }

    free(in_buff);
}

void ReLU(float *x, int n)
{
    for (int u = 0; u < n; u++)
        for (int i = 0; i < layer_output; i++)
            x[u * layer_output + i] = x[u * layer_output + i] > 0.0f ? x[u * layer_output + i] : 0.0f;
}

void Sigmoid(float *x, int n)
{
    for (int u = 0; u < n; u++)
        for (int i = 0; i < layer_output; i++)
            x[u * layer_output + i] = 1.0f / (1.0f + expf(-x[u * layer_output + i]));
}

void prep_params(const float *p, int n, int m, float *out_p, float *out_b)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < layer_output; j++)
            out_p[i * layer_output + j] = 0.0f;

    for (int j = 0; j < layer_output; j++)
        out_b[j] = 0.0f;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            out_p[i * layer_output + j] = p[j * n + i];

    for (int j = 0; j < m; j++)
        out_b[j] = p[n * m + j];
}

float *lr_gcn_predict(graph_access &g, std::string path)
{
    float *params = NULL;
    FILE *f = fopen(path.c_str(), "r");
    lr_gcn_parse(f, &params);
    fclose(f);

    int N = g.number_of_nodes();
    float *x = (float *)aligned_alloc(32, sizeof(float) * max_input * N);
    float *y = (float *)aligned_alloc(32, sizeof(float) * max_input * N);
    float *cat = (float *)aligned_alloc(32, sizeof(float) * cat_input * N);
    float *p = (float *)aligned_alloc(32, sizeof(float) * cat_input * layer_output);
    float *b = (float *)aligned_alloc(32, sizeof(float) * layer_output);

    memset(x, 0, sizeof(float) * max_input * N);
    memset(y, 0, sizeof(float) * max_input * N);
    memset(cat, 0, sizeof(float) * cat_input * N);

    // Compute input features
    compute_features(g, x);

    // Store input in cat for output layers
    for (int i = 0; i < N; i++)
        for (int j = 0; j < layer_inputs[0]; j++)
            cat[i * cat_input + j] = x[i * layer_inputs[0] + j];

    // Copy and transpose first layer parameters and run lr_conv
    prep_params(params + param_offsets[0], layer_inputs[0] * 2, layer_outputs[0], p, b);
    lr_conv_combined_avx(g, layer_inputs[0], x, y, p, b);
    ReLU(y, N);

    // Store first layer output in cat for output layers
    for (int i = 0; i < N; i++)
        for (int j = 0; j < layer_output; j++)
            cat[i * cat_input + layer_inputs[0] + j] = y[i * layer_output + j];

    // Copy and transpose second layer parameters and run lr_conv
    prep_params(params + param_offsets[1], layer_inputs[1] * 2, layer_outputs[1], p, b);
    lr_conv_combined_avx(g, layer_inputs[1], y, x, p, b);
    ReLU(x, N);

    // Store second layer output in cat for output layers
    for (int i = 0; i < N; i++)
        for (int j = 0; j < layer_output; j++)
            cat[i * cat_input + layer_inputs[0] + layer_inputs[1] + j] = x[i * layer_output + j];

    // Copy and transpose third layer parameters and run linear
    prep_params(params + param_offsets[2], layer_inputs[2], layer_outputs[2], p, b);
    matrix_multiplication_avx(cat, p, y, b, N, layer_inputs[2]);
    ReLU(y, N);

    // Copy and transpose forth layer parameters and run linear
    prep_params(params + param_offsets[3], layer_inputs[3], layer_outputs[3], p, b);
    matrix_multiplication_avx(y, p, x, b, N, layer_inputs[3]);
    Sigmoid(x, N);

    float *res = (float *)malloc(sizeof(float) * N);
    for (int i = 0; i < N; i++)
        res[i] = x[i * layer_output];

    free(params);
    free(x);
    free(y);
    free(cat);
    free(p);
    free(b);

    return res;
}