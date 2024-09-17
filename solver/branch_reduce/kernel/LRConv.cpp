#include "LRConv.h"
#include "branch_and_reduce_algorithm.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <immintrin.h>

const int W_OFF[6] = {0, 656, 928, 1584, 1856, 2128};
const int B_OFF[6] = {640, 912, 1568, 1840, 2112, 2144};

const int total_params = 2145;

LRConv::LRConv(int N)
    : params(), W(NULL), B(NULL), param(NULL), x(NULL), y(NULL)
{
    W = (float **)malloc(sizeof(float *) * 6);
    B = (float **)malloc(sizeof(float *) * 6);
    param = (float *)aligned_alloc(32, sizeof(float) * total_params);

    for (int i = 0; i < 6; i++)
    {
        W[i] = param + W_OFF[i];
        B[i] = param + B_OFF[i];
    }

    x = NULL;
    y = NULL;
    e = NULL;
}

// LRConv::LRConv(const LRConv &gnn)
//     : params(), W(NULL), B(NULL), param(NULL), allocated(gnn.allocated), x(NULL), y(NULL)
// {
//     W = (float **)malloc(sizeof(float *) * 6);
//     B = (float **)malloc(sizeof(float *) * 6);
//     param = (float *)aligned_alloc(32, sizeof(float) * total_params);

//     for (int i = 0; i < total_params; i++)
//         param[i] = gnn.param[i];

//     for (int i = 0; i < 6; i++)
//     {
//         W[i] = param + W_OFF[i];
//         B[i] = param + B_OFF[i];
//     }

//     x = (float *)aligned_alloc(32, sizeof(float) * allocated);
//     y = (float *)aligned_alloc(32, sizeof(float) * allocated);
// }

LRConv::~LRConv()
{
    free(W);
    free(B);
    free(param);
    free(x);
    free(y);
    free(e);
}

void transpose(std::vector<float> &in, int offset, float *out, int N, int M)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            out[i * M + j] = in[offset + j * N + i];
}

void LRConv::change_parameters(const std::string path)
{
    // printf("%s\n", path.c_str());
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

    transpose(p, W_OFF[0], W[0], 40, 16);
    transpose(p, B_OFF[0], B[0], 1, 16);
    transpose(p, W_OFF[1], W[1], 16, 16);
    transpose(p, B_OFF[1], B[1], 1, 16);
    transpose(p, W_OFF[2], W[2], 40, 16);
    transpose(p, B_OFF[2], B[2], 1, 16);
    transpose(p, W_OFF[3], W[3], 16, 16);
    transpose(p, B_OFF[3], B[3], 1, 16);
    transpose(p, W_OFF[4], W[4], 16, 16);
    transpose(p, B_OFF[4], B[4], 1, 16);
    transpose(p, W_OFF[5], W[5], 1, 16);
    transpose(p, B_OFF[5], B[5], 1, 1);
}

// static inline void blas_kernel_6_16(const float *A, int lda, const float *B, const float *bias, float *C)
// {
//     __m256 c00 = _mm256_load_ps(bias);
//     __m256 c01 = _mm256_load_ps(bias + 8);
//     __m256 c10 = _mm256_load_ps(bias);
//     __m256 c11 = _mm256_load_ps(bias + 8);
//     __m256 c20 = _mm256_load_ps(bias);
//     __m256 c21 = _mm256_load_ps(bias + 8);
//     __m256 c30 = _mm256_load_ps(bias);
//     __m256 c31 = _mm256_load_ps(bias + 8);
//     __m256 c40 = _mm256_load_ps(bias);
//     __m256 c41 = _mm256_load_ps(bias + 8);
//     __m256 c50 = _mm256_load_ps(bias);
//     __m256 c51 = _mm256_load_ps(bias + 8);

//     for (int i = 0; i < lda; i++)
//     {
//         __m256 b0 = _mm256_load_ps(B + i * 16);
//         __m256 b1 = _mm256_load_ps(B + 8 + i * 16);

//         __m256 a0 = _mm256_broadcast_ss(A + i + (lda * 0));
//         c00 = _mm256_fmadd_ps(a0, b0, c00);
//         c01 = _mm256_fmadd_ps(a0, b1, c01);

//         __m256 a1 = _mm256_broadcast_ss(A + i + (lda * 1));
//         c10 = _mm256_fmadd_ps(a1, b0, c10);
//         c11 = _mm256_fmadd_ps(a1, b1, c11);

//         __m256 a2 = _mm256_broadcast_ss(A + i + (lda * 2));
//         c20 = _mm256_fmadd_ps(a2, b0, c20);
//         c21 = _mm256_fmadd_ps(a2, b1, c21);

//         __m256 a3 = _mm256_broadcast_ss(A + i + (lda * 3));
//         c30 = _mm256_fmadd_ps(a3, b0, c30);
//         c31 = _mm256_fmadd_ps(a3, b1, c31);

//         __m256 a4 = _mm256_broadcast_ss(A + i + (lda * 4));
//         c40 = _mm256_fmadd_ps(a4, b0, c40);
//         c41 = _mm256_fmadd_ps(a4, b1, c41);

//         __m256 a5 = _mm256_broadcast_ss(A + i + (lda * 5));
//         c50 = _mm256_fmadd_ps(a5, b0, c50);
//         c51 = _mm256_fmadd_ps(a5, b1, c51);
//     }

//     __m256 zero = _mm256_setzero_ps();

//     c00 = _mm256_max_ps(c00, zero);
//     c01 = _mm256_max_ps(c01, zero);
//     c10 = _mm256_max_ps(c10, zero);
//     c11 = _mm256_max_ps(c11, zero);
//     c20 = _mm256_max_ps(c20, zero);
//     c21 = _mm256_max_ps(c21, zero);
//     c30 = _mm256_max_ps(c30, zero);
//     c31 = _mm256_max_ps(c31, zero);
//     c40 = _mm256_max_ps(c40, zero);
//     c41 = _mm256_max_ps(c41, zero);
//     c50 = _mm256_max_ps(c50, zero);
//     c51 = _mm256_max_ps(c51, zero);

//     _mm256_store_ps(C, c00);
//     _mm256_store_ps(C + 8, c01);
//     _mm256_store_ps(C + 16, c10);
//     _mm256_store_ps(C + 24, c11);
//     _mm256_store_ps(C + 32, c20);
//     _mm256_store_ps(C + 40, c21);
//     _mm256_store_ps(C + 48, c30);
//     _mm256_store_ps(C + 56, c31);
//     _mm256_store_ps(C + 64, c40);
//     _mm256_store_ps(C + 72, c41);
//     _mm256_store_ps(C + 80, c50);
//     _mm256_store_ps(C + 88, c51);
// }

// void LRConv_message_blas(const float *x, int dim, float *y,
//                          const float *e, int dim_e,
//                          const float *W1, const float *B1,
//                          const float *W2, const float *B2,
//                          graph_access &g)
// {
//     int lda = (dim * 2) + dim_e;
//     float *A = (float *)aligned_alloc(32, sizeof(float) * lda * 6);
//     float *C1 = (float *)aligned_alloc(32, sizeof(float) * 16 * 6);
//     float *C2 = (float *)aligned_alloc(32, sizeof(float) * 16 * 6);

//     for (int u = 0; u < g.number_of_nodes(); u++)
//     {
//         for (int i = 0; i < 6; i++)
//             for (int j = 0; j < dim; j++)
//                 A[i * lda + j] = x[u * dim + j];

//         for (int i = 0; i < 16; i++)
//             y[u * 16 + i] = 0.0f;

//         int t = 0;
//         for (int ei = g.get_first_edge(u); ei != g.get_first_invalid_edge(u); ei++)
//         {
//             int v = g.getEdgeTarget(ei);

//             for (int i = 0; i < dim; i++)
//                 A[t * lda + dim + i] = x[v * dim + i];
//             for (int i = 0; i < dim_e; i++)
//                 A[t * lda + 2 * dim + i] = e[ei * dim_e + i];
//             t++;

//             if (t == 6)
//             {
//                 t = 0;
//                 blas_kernel_6_16(A, lda, W1, B1, C1);
//                 blas_kernel_6_16(C1, 16, W2, B2, C2);

//                 for (int i = 0; i < 6; i++)
//                     for (int j = 0; j < 16; j++)
//                         if (y[u * 16 + j] < C2[i * 16 + j])
//                             y[u * 16 + j] = C2[i * 16 + j];
//             }
//         }

//         if (t > 0)
//         {
//             blas_kernel_6_16(A, lda, W1, B1, C1);
//             blas_kernel_6_16(C1, 16, W2, B2, C2);

//             for (int i = 0; i < t; i++)
//                 for (int j = 0; j < 16; j++)
//                     if (y[u * 16 + j] < C2[i * 16 + j])
//                         y[u * 16 + j] = C2[i * 16 + j];
//         }
//     }

//     free(A);
//     free(C1);
//     free(C2);
// }

void LRConv_message(const float *x, int dim_in,
                    const float *e, int dim_edge,
                    float *y, int dim_out,
                    const float *param1, const float *bias1,
                    const float *param2, const float *bias2,
                    graph_access &g)
{
    float *buff = (float *)aligned_alloc(32, sizeof(float) * (dim_in * 2 + dim_edge));
    float *tmp1 = (float *)aligned_alloc(32, sizeof(float) * dim_out);
    float *tmp2 = (float *)aligned_alloc(32, sizeof(float) * dim_out);

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        float *agg = y + dim_out * u;
        for (int i = 0; i < dim_out; i++)
            agg[i] = 0.0f;

        for (int i = 0; i < dim_in; i++)
            buff[i] = x[u * dim_in + i];

        for (int ei = g.get_first_edge(u); ei != g.get_first_invalid_edge(u); ei++)
        {
            int v = g.getEdgeTarget(ei);
            for (int i = 0; i < dim_in; i++)
                buff[dim_in + i] = x[v * dim_in + i];

            for (int i = 0; i < dim_edge; i++)
                buff[dim_in * 2 + i] = e[ei * dim_edge + i];

            for (int i = 0; i < dim_out; i++)
            {
                tmp1[i] = bias1[i];
                tmp2[i] = bias2[i];
            }

            // Layer 1
            for (int i = 0; i < dim_out; i++)
                for (int j = 0; j < dim_in * 2 + dim_edge; j++)
                    tmp1[i] += buff[j] * param1[i * (dim_in * 2 + dim_edge) + j];

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

void LRConv_dense(const float *x, int dim_in,
                  float *y, int dim_out, int M,
                  const float *param, const float *bias, bool act)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < dim_out; j++)
        {
            y[i * dim_out + j] = bias[j];
            for (int k = 0; k < dim_in; k++)
            {
                y[i * dim_out + j] += x[i * dim_in + k] * param[k * dim_out + j];
            }
            if (act && y[i * dim_out + j] < 0.0f)
                y[i * dim_out + j] = 0.0f;
        }
    }
}

const float *LRConv::predict(branch_and_reduce_algorithm *br_alg)
{
    graph_access g;
    std::vector<NodeID> rm(br_alg->status.n);
    br_alg->build_graph_access(g, rm);

    // printf("%d %d\n", g.number_of_nodes(), g.number_of_edges());

    LRConv::compute_attr_norm(&x, &e, g);
    free(y);
    y = (float *)aligned_alloc(32, sizeof(float) * g.number_of_nodes() * total_node_features);

    LRConv_message(x, total_node_features, e, total_edge_features,
                   y, hidden_dim, W[0], B[0], W[1], B[1], g);

    LRConv_message(y, hidden_dim, e, total_edge_features,
                   x, hidden_dim, W[2], B[2], W[3], B[3], g);

    LRConv_dense(x, hidden_dim, y, hidden_dim, g.number_of_nodes(), W[4], B[4], true);
    LRConv_dense(y, hidden_dim, x, 1, g.number_of_nodes(), W[5], B[5], false);

    free(y);
    y = (float *)aligned_alloc(32, sizeof(float) * br_alg->status.n);

    for (int i = 0; i < br_alg->status.n; i++)
        y[i] = -99999.0f;

    for (int i = 0; i < g.number_of_nodes(); i++)
        y[rm[i]] = x[i];

    return y;
}

static inline unsigned int hash(unsigned int x)
{
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

void LRConv::compute_attr_norm(float **node_attr, float **edge_attr, graph_access &g)
{
    free(*node_attr);
    *node_attr = (float *)aligned_alloc(32, sizeof(float) * g.number_of_nodes() * total_node_features);
    free(*edge_attr);
    *edge_attr = (float *)aligned_alloc(32, sizeof(float) * g.number_of_edges() * total_edge_features);

    int tW = 0, tD = 0, maW = 0, maD = 0;
    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        tW += g.getNodeWeight(u);
        tD += g.getNodeDegree(u);

        maW = std::max(maW, (int)g.getNodeWeight(u));
        maD = std::max(maD, (int)g.getNodeDegree(u));
    }

    float aW = (float)tW / (float)g.number_of_nodes();
    float aD = (float)tD / (float)g.number_of_nodes();

    for (int u = 0; u < g.number_of_nodes(); u++)
    {
        int d = g.getNodeDegree(u), Dn = 0, Wn = 0,
            maDn = 0, miDn = std::numeric_limits<int>::max(),
            maWn = 0, miWn = std::numeric_limits<int>::max();

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);
            Wn += g.getNodeWeight(v);
            Dn += g.getNodeDegree(v);

            maWn = std::max(maWn, (int)g.getNodeWeight(v));
            maDn = std::max(maDn, (int)g.getNodeDegree(v));

            miWn = std::min(miWn, (int)g.getNodeWeight(v));
            miDn = std::min(miDn, (int)g.getNodeDegree(v));
        }

        for (int i = g.get_first_edge(u); i != g.get_first_invalid_edge(u); i++)
        {
            int v = g.getEdgeTarget(i);

            float *ea = *edge_attr + (i * total_edge_features);
            ea[0] = (float)(g.getNodeWeight(u) > g.getNodeWeight(v));
            ea[1] = (float)(g.getNodeDegree(u) > g.getNodeDegree(v));
            ea[2] = (float)(g.getNodeWeight(u) < g.getNodeWeight(v));
            ea[3] = (float)(g.getNodeDegree(u) < g.getNodeDegree(v));
            ea[4] = (float)(g.getNodeWeight(u) == g.getNodeWeight(v));
            ea[5] = (float)(g.getNodeDegree(u) == g.getNodeDegree(v));
            ea[6] = (float)g.getNodeWeight(u) / (float)g.getNodeWeight(v);
            ea[7] = (float)g.getNodeDegree(u) / (float)g.getNodeDegree(v);
        }

        float *ua = *node_attr + (u * total_node_features);
        ua[0] = (float)(hash(u + 1) % (int)id_scale) / id_scale;
        ua[1] = g.getNodeWeight(u);
        ua[2] = g.getNodeDegree(u);
        ua[3] = Wn;
        ua[4] = d > 0 ? miWn : 0.0f;
        ua[5] = maWn;
        ua[6] = d > 0 ? miDn : 0.0f;
        ua[7] = maDn;
        ua[8] = d > 0 ? (float)g.getNodeWeight(u) / (float)Wn : 0.0f;
        ua[9] = d > 0 ? (float)g.getNodeWeight(u) / ((float)Wn / (float)d) : 0.0f;
        ua[10] = d > 0 ? (float)g.getNodeDegree(u) / (float)Dn : 0.0f;
        ua[11] = d > 0 ? (float)g.getNodeDegree(u) / ((float)Dn / (float)d) : 0.0f;

        ua[12] = (float)g.getNodeWeight(u) / aW;
        ua[13] = (float)g.getNodeWeight(u) / (float)maW;
        ua[14] = (float)g.getNodeDegree(u) / aD;
        ua[15] = (float)g.getNodeDegree(u) / (float)maD;
    }
}
