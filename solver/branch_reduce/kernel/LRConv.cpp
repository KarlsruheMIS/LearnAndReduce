#include "LRConv.h"
#include "branch_and_reduce_algorithm.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <immintrin.h>

const int W_OFF[6] = {0, 144, 416, 944, 1216, 1488};
const int B_OFF[6] = {128, 400, 928, 1200, 1472, 1504};

const int total_params = 1520;

LRConv::LRConv(int N)
    : params(), W(NULL), B(NULL), param(NULL), allocated(0), x(NULL), y(NULL)
{
    W = (float **)malloc(sizeof(float *) * 6);
    B = (float **)malloc(sizeof(float *) * 6);
    param = (float *)aligned_alloc(32, sizeof(float) * total_params);

    for (int i = 0; i < 6; i++)
    {
        W[i] = param + W_OFF[i];
        B[i] = param + B_OFF[i];
    }

    allocated = N * hidden_dim;
    x = (float *)aligned_alloc(32, sizeof(float) * allocated);
    y = (float *)aligned_alloc(32, sizeof(float) * allocated);
}

LRConv::LRConv(const LRConv &gnn)
    : params(), W(NULL), B(NULL), param(NULL), allocated(gnn.allocated), x(NULL), y(NULL)
{
    W = (float **)malloc(sizeof(float *) * 6);
    B = (float **)malloc(sizeof(float *) * 6);
    param = (float *)aligned_alloc(32, sizeof(float) * total_params);

    for (int i = 0; i < total_params; i++)
        param[i] = gnn.param[i];

    for (int i = 0; i < 6; i++)
    {
        W[i] = param + W_OFF[i];
        B[i] = param + B_OFF[i];
    }

    x = (float *)aligned_alloc(32, sizeof(float) * allocated);
    y = (float *)aligned_alloc(32, sizeof(float) * allocated);
}

LRConv::~LRConv()
{
    free(W);
    free(B);
    free(param);
    free(x);
    free(y);
}

void transpose(std::vector<float> &in, int offset, float *out, int N, int M)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            out[i * M + j] = in[offset + j * N + i];
}

void LRConv::change_parameters(const std::string path)
{
    if (params.find(path) == params.end())
    {
        int n_params;
        std::ifstream file;
        file.open(path);

        file >> n_params;
        std::vector<float> res(total_params, 0.0f);

        for (int i = 0; i < n_params; i++)
            file >> res[i];

        file.close();
        params[path] = res;
    }

    auto &p = params[path];

    transpose(p, W_OFF[0], W[0], 8, 16);
    transpose(p, B_OFF[0], B[0], 1, 16);
    transpose(p, W_OFF[1], W[1], 16, 16);
    transpose(p, B_OFF[1], B[1], 1, 16);
    transpose(p, W_OFF[2], W[2], 32, 16);
    transpose(p, B_OFF[2], B[2], 1, 16);
    transpose(p, W_OFF[3], W[3], 16, 16);
    transpose(p, B_OFF[3], B[3], 1, 16);
    transpose(p, W_OFF[4], W[4], 16, 16);
    transpose(p, B_OFF[4], B[4], 1, 16);
    transpose(p, W_OFF[5], W[5], 1, 16);
    transpose(p, B_OFF[5], B[5], 1, 16);
}

static inline void blas_kernel_6_16(const float *A, int lda, const float *B, const float *bias, float *C)
{
    __m256 c00 = _mm256_load_ps(bias);
    __m256 c01 = _mm256_load_ps(bias + 8);
    __m256 c10 = _mm256_load_ps(bias);
    __m256 c11 = _mm256_load_ps(bias + 8);
    __m256 c20 = _mm256_load_ps(bias);
    __m256 c21 = _mm256_load_ps(bias + 8);
    __m256 c30 = _mm256_load_ps(bias);
    __m256 c31 = _mm256_load_ps(bias + 8);
    __m256 c40 = _mm256_load_ps(bias);
    __m256 c41 = _mm256_load_ps(bias + 8);
    __m256 c50 = _mm256_load_ps(bias);
    __m256 c51 = _mm256_load_ps(bias + 8);

    for (int i = 0; i < lda; i++)
    {
        __m256 b0 = _mm256_load_ps(B + i * 16);
        __m256 b1 = _mm256_load_ps(B + 8 + i * 16);

        __m256 a0 = _mm256_broadcast_ss(A + i + (lda * 0));
        c00 = _mm256_fmadd_ps(a0, b0, c00);
        c01 = _mm256_fmadd_ps(a0, b1, c01);

        __m256 a1 = _mm256_broadcast_ss(A + i + (lda * 1));
        c10 = _mm256_fmadd_ps(a1, b0, c10);
        c11 = _mm256_fmadd_ps(a1, b1, c11);

        __m256 a2 = _mm256_broadcast_ss(A + i + (lda * 2));
        c20 = _mm256_fmadd_ps(a2, b0, c20);
        c21 = _mm256_fmadd_ps(a2, b1, c21);

        __m256 a3 = _mm256_broadcast_ss(A + i + (lda * 3));
        c30 = _mm256_fmadd_ps(a3, b0, c30);
        c31 = _mm256_fmadd_ps(a3, b1, c31);

        __m256 a4 = _mm256_broadcast_ss(A + i + (lda * 4));
        c40 = _mm256_fmadd_ps(a4, b0, c40);
        c41 = _mm256_fmadd_ps(a4, b1, c41);

        __m256 a5 = _mm256_broadcast_ss(A + i + (lda * 5));
        c50 = _mm256_fmadd_ps(a5, b0, c50);
        c51 = _mm256_fmadd_ps(a5, b1, c51);
    }

    __m256 zero = _mm256_setzero_ps();

    c00 = _mm256_max_ps(c00, zero);
    c01 = _mm256_max_ps(c01, zero);
    c10 = _mm256_max_ps(c10, zero);
    c11 = _mm256_max_ps(c11, zero);
    c20 = _mm256_max_ps(c20, zero);
    c21 = _mm256_max_ps(c21, zero);
    c30 = _mm256_max_ps(c30, zero);
    c31 = _mm256_max_ps(c31, zero);
    c40 = _mm256_max_ps(c40, zero);
    c41 = _mm256_max_ps(c41, zero);
    c50 = _mm256_max_ps(c50, zero);
    c51 = _mm256_max_ps(c51, zero);

    _mm256_store_ps(C, c00);
    _mm256_store_ps(C + 8, c01);
    _mm256_store_ps(C + 16, c10);
    _mm256_store_ps(C + 24, c11);
    _mm256_store_ps(C + 32, c20);
    _mm256_store_ps(C + 40, c21);
    _mm256_store_ps(C + 48, c30);
    _mm256_store_ps(C + 56, c31);
    _mm256_store_ps(C + 64, c40);
    _mm256_store_ps(C + 72, c41);
    _mm256_store_ps(C + 80, c50);
    _mm256_store_ps(C + 88, c51);
}

void LRConv_message_blas_dyn(const float *x, int dim, float *y,
                             const float *W1, const float *B1,
                             const float *W2, const float *B2,
                             dynamic_graph &g, std::vector<branch_and_reduce_algorithm::IS_status> &status)
{
    float *A = (float *)aligned_alloc(32, sizeof(float) * (dim * 2) * 6);
    float *C1 = (float *)aligned_alloc(32, sizeof(float) * 16 * 6);
    float *C2 = (float *)aligned_alloc(32, sizeof(float) * 16 * 6);

    for (int u = 0; u < g.size(); u++)
    {
        if (status[u] != branch_and_reduce_algorithm::IS_status::not_set)
            continue;

        for (int i = 0; i < 6; i++)
            for (int j = 0; j < dim; j++)
                A[i * (dim * 2) + j] = x[u * dim + j];

        for (int i = 0; i < 16; i++)
            y[u * 16 + i] = 0.0f;

        int t = 0;
        for (auto v : g[u])
        {
            if (status[v] != branch_and_reduce_algorithm::IS_status::not_set)
                continue;

            for (int i = 0; i < dim; i++)
                A[t * (dim * 2) + dim + i] = x[v * dim + i];
            t++;

            if (t == 6)
            {
                t = 0;
                blas_kernel_6_16(A, dim * 2, W1, B1, C1);
                blas_kernel_6_16(C1, 16, W2, B2, C2);

                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 16; j++)
                        if (y[u * 16 + j] < C2[i * 16 + j])
                            y[u * 16 + j] = C2[i * 16 + j];
            }
        }

        if (t > 0)
        {
            blas_kernel_6_16(A, dim * 2, W1, B1, C1);
            blas_kernel_6_16(C1, 16, W2, B2, C2);

            for (int i = 0; i < t; i++)
                for (int j = 0; j < 16; j++)
                    if (y[u * 16 + j] < C2[i * 16 + j])
                        y[u * 16 + j] = C2[i * 16 + j];
        }
    }

    free(A);
    free(C1);
    free(C2);
}

const float *LRConv::predict(branch_and_reduce_algorithm *br_alg)
{
    auto &g = br_alg->status.graph;
    auto &status = br_alg->status.node_status;

    if (br_alg->status.graph.size() * hidden_dim > allocated)
    {
        allocated = br_alg->status.graph.size() * hidden_dim;
        free(x);
        free(y);
        x = (float *)aligned_alloc(32, sizeof(float) * allocated);
        y = (float *)aligned_alloc(32, sizeof(float) * allocated);
    }

    compute_node_attr_dynamic(br_alg);

    LRConv_message_blas_dyn(x, node_features, y, W[0], B[0], W[1], B[1], g, status);
    LRConv_message_blas_dyn(y, hidden_dim, x, W[2], B[2], W[3], B[3], g, status);

    float *tmp = (float *)aligned_alloc(32, sizeof(float) * hidden_dim);

    for (int u = 0; u < g.size(); u++)
    {
        if (status[u] != branch_and_reduce_algorithm::IS_status::not_set)
            continue;

        y[u] = *B[5];

        for (int i = 0; i < hidden_dim; i++)
            tmp[i] = B[4][i];

        // Layer 1
        for (int i = 0; i < hidden_dim; i++)
            for (int j = 0; j < hidden_dim; j++)
                tmp[j] += x[u * hidden_dim + i] * W[4][i * hidden_dim + j];

        for (int i = 0; i < hidden_dim; i++)
            tmp[i] = tmp[i] >= 0.0f ? tmp[i] : 0.0f;

        // Layer 2
        for (int i = 0; i < hidden_dim; i++)
            y[u] += tmp[i] * W[5][i];
    }

    free(tmp);

    return y;
}

static inline unsigned int hash(unsigned int x)
{
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
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
        int d = g.getNodeDegree(u), Dn = 0, Wn = 0, C = 0;

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);
            marks[v] = u;
            Wn += g.getNodeWeight(v);
            Dn += g.getNodeDegree(v);
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
        if (d == 0)
        {
            ua[0] = 0.0f;
            ua[1] = 0.0f;
            ua[2] = 2.0f;
        }
        else
        {
            ua[0] = (float)d / ((float)Dn / (float)d);
            ua[1] = (float)g.getNodeWeight(u) / ((float)Wn / (float)d);
            ua[2] = (float)g.getNodeWeight(u) / (float)Wn;
        }
        ua[3] = (float)(hash(u + 1) % 1000000) / id_scale;
        // ua[0] = (float)d / scale;
        // ua[1] = (float)g.getNodeWeight(u) / scale;
        // ua[2] = (float)Wn / scale;
        // ua[3] = (float)(hash(u + 1) % 1000000) / id_scale;
    }
}

void LRConv::compute_node_attr_dynamic(branch_and_reduce_algorithm *br_alg)
{
    auto &g = br_alg->status.graph;
    auto &status = br_alg->status.node_status;
    auto &weights = br_alg->status.weights;

    for (int u = 0; u < g.size(); u++)
    {
        if (status[u] != branch_and_reduce_algorithm::IS_status::not_set)
            continue;

        int d = g[u].size(), Wn = 0;
        for (auto v : g[u])
            Wn += weights[v];

        float *ua = x + (u * node_features);
        ua[0] = (float)d / scale;
        ua[1] = (float)weights[u] / scale;
        ua[2] = (float)Wn / scale;
        ua[3] = (float)(hash(u + 1) % 1000000) / id_scale;
    }
}

void LRConv::compute_node_attr_dynamic_norm(branch_and_reduce_algorithm *br_alg)
{
    auto &g = br_alg->status.graph;
    auto &status = br_alg->status.node_status;
    auto &weights = br_alg->status.weights;

    for (int u = 0; u < g.size(); u++)
    {
        if (status[u] != branch_and_reduce_algorithm::IS_status::not_set)
            continue;

        int d = g[u].size(), Dn = 0, Wn = 0;
        for (auto v : g[u])
        {
            Wn += weights[v];
            Dn += g[v].size();
        }

        float *ua = x + (u * node_features);
        ua[0] = (float)d / ((float)Dn / (float)d);
        ua[1] = (float)weights[u] / ((float)Wn / (float)d);
        ua[2] = (float)weights[u] / (float)Wn;
        ua[3] = (float)(hash(u + 1) % 1000000) / id_scale;
    }
}