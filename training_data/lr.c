#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <immintrin.h>
#include <math.h>

typedef struct
{
    int N;
    int *V, *E;
    long long *W;
} graph;

static inline void parse_id(char *data, size_t *p, long long *v)
{
    while (data[*p] < '0' || data[*p] > '9')
        (*p)++;

    *v = 0;
    while (data[*p] >= '0' && data[*p] <= '9')
        *v = (*v) * 10 + data[(*p)++] - '0';
}

graph *graph_parse(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);

    char *data = mmap(0, size, PROT_READ, MAP_PRIVATE, fileno_unlocked(f), 0);
    size_t p = 0;

    long long N, M, t;
    parse_id(data, &p, &N);
    parse_id(data, &p, &M);
    parse_id(data, &p, &t);

    int *V = malloc(sizeof(int) * (N + 1));
    int *E = malloc(sizeof(int) * (M * 2));

    long long *W = malloc(sizeof(long long) * N);

    int ei = 0;
    for (int u = 0; u < N; u++)
    {
        parse_id(data, &p, W + u);
        V[u] = ei;
        while (ei < M * 2)
        {
            while (data[p] == ' ')
                p++;
            if (data[p] == '\n')
                break;

            long long e;
            parse_id(data, &p, &e);
            E[ei++] = e - 1;
            ;
        }
        p++;
    }
    V[N] = ei;

    munmap(data, size);

    graph *g = malloc(sizeof(graph));
    *g = (graph){.N = N, .V = V, .E = E, .W = W};

    return g;
}

graph *graph_parse_csv(FILE *fg, FILE *fr)
{
    int rc = fscanf(fg, "source;target\n");

    int allocated = 64;
    int *X = malloc(sizeof(int) * allocated);
    int *Y = malloc(sizeof(int) * allocated);

    int u, v;

    int M = 0, N = 0;
    while (1)
    {
        rc = fscanf(fg, "%d;%d\n", &u, &v);
        if (rc != 2)
            break;

        X[M] = u;
        Y[M] = v;

        if (u + 1 > N)
            N = u + 1;
        if (v + 1 > N)
            N = v + 1;

        M++;

        if (M >= allocated)
        {
            allocated *= 2;
            X = realloc(X, sizeof(int) * allocated);
            Y = realloc(Y, sizeof(int) * allocated);
        }
    }

    int *V = malloc(sizeof(int) * (N + 1));
    int *E = malloc(sizeof(int) * M);

    for (int i = 0; i <= N; i++)
        V[i] = 0;

    for (int i = 0; i < M; i++)
    {
        E[i] = Y[i];
        V[X[i] + 1]++;
    }

    for (int i = 1; i <= N; i++)
        V[i] += V[i - 1];

    rc = fscanf(fr, "id;w;neighborhood;degree_1;degree_2;simplicial_vertex;domination;single_edge;extended_single_edge;twin;extended_twin;funnel;unconfined_csr;critical_set;generalized_fold;heavy_set;heavy_set3;cut_vertex\n");

    long long *W = malloc(sizeof(long long) * N);

    int nc = 0;

    for (int i = 0; i < N; i++)
    {
        int id, w, nb, d1, d2, sv, dom, si, esi, tw, etw, fu, uc, cr, gf, hs, hs3, cv;
        rc = fscanf(fr, "%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d;%d\n",
                    &id, &w, &nb, &d1, &d2, &sv, &dom, &si, &esi, &tw, &etw, &fu, &uc, &cr, &gf, &hs, &hs3, &cv);

        nc += cr;

        W[id] = w;
    }

    printf("%d;", nc);

    free(X);
    free(Y);

    graph *g = malloc(sizeof(graph));
    g->N = N;
    g->V = V;
    g->E = E;
    g->W = W;

    return g;
}

void graph_free(graph *g)
{
    if (g == NULL)
        return;

    free(g->V);
    free(g->E);
    free(g->W);

    free(g);
}

int lr_gcn_parse(FILE *f, float **paramv)
{
    int paramc;
    int np = fscanf(f, "%d\n", &paramc);

    *paramv = malloc(sizeof(float) * paramc);
    for (int i = 0; i < paramc; i++)
        np = fscanf(f, "%f\n", (*paramv) + i);

    return paramc;
}

/*
    m is allways equal to 16
    Assumes A, B, and C are already in aligned memory
 */
void matrix_multiplication_ref(const float *A, const float *B, float *C, const float *b, int n, int k)
{
    const int m = 16;
    for (int i1 = 0; i1 < n; i1++)
    {
        for (int i2 = 0; i2 < m; i2++)
            C[i1 * m + i2] = b[i2];

        for (int i3 = 0; i3 < k; i3++)
        {
            for (int i2 = 0; i2 < m; i2++)
            {
                C[i1 * m + i2] += A[i1 * k + i3] * B[i3 * m + i2];
            }
        }
    }
}

void compute_features(const graph *g, float *x)
{
    for (int u = 0; u < g->N; u++)
    {
        for (int i = 0; i < 8; i++)
            x[u * 8 + i] = 0.0;

        x[u * 8] = g->W[u];
        x[u * 8 + 1] = g->W[u];
        x[u * 8 + 4] = g->V[u + 1] - g->V[u];

        float scale = 1.0f / x[u * 8 + 4];
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            float d = g->V[v + 1] - g->V[v];
            float w = g->W[v];

            x[u * 8 + 1] -= w;
            if (w < x[u * 8 + 2] || x[u * 8 + 2] == 0.0f)
                x[u * 8 + 2] = w;
            if (w > x[u * 8 + 3])
                x[u * 8 + 3] = w;

            x[u * 8 + 5] += d * scale;
            if (d < x[u * 8 + 6] || x[u * 8 + 6] == 0.0f)
                x[u * 8 + 6] = d;
            if (d > x[u * 8 + 7])
                x[u * 8 + 7] = d;
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 8; j++)
            x[(g->N + i) * 8 + j] = 0.0f;
}

void lr_conv(const graph *g, int dim, float *t1, float *t2, const float *p, const float *b)
{
    for (int u = 0; u < g->N; u++)
    {
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            int v = g->E[i];
            for (int j = 0; j < dim; j++)
                t2[i * dim * 2 + j] = t1[u * dim + j];
            for (int j = 0; j < dim; j++)
                t2[i * dim * 2 + dim + j] = t1[v * dim + j];
        }
    }

    matrix_multiplication_ref(t2, p, t1, b, g->V[g->N], dim * 2);

    for (int u = 0; u < g->N; u++)
    {
        float scale = 1.0f / (float)(g->V[u + 1] - g->V[u]);
        for (int j = 0; j < 16; j++)
            t2[u * 16 + j] = 0.0f;
        for (int i = g->V[u]; i < g->V[u + 1]; i++)
        {
            for (int j = 0; j < 16; j++)
                t2[u * 16 + j] += t1[i * 16 + j] * scale;
        }
        for (int j = 0; j < 16; j++)
            t2[u * 16 + j] = t2[u * 16 + j] > 0.0f ? t2[u * 16 + j] : 0.0f;
    }
}

void ReLU(float *x, int n)
{
    for (int u = 0; u < n; u++)
        for (int i = 0; i < 16; i++)
            x[u * 16 + i] = x[u * 16 + i] > 0.0f ? x[u * 16 + i] : 0.0f;
}

void Sigmoid(float *x, int n)
{
    for (int u = 0; u < n; u++)
        for (int i = 0; i < 16; i++)
            x[u * 16 + i] = 1.0f / (1.0f + expf(-x[u * 16 + i]));
}

void prep_params(const float *p, int n, int m, float *out_p, float *out_b)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 16; j++)
            out_p[i * 16 + j] = 0.0f;

    for (int j = 0; j < 16; j++)
        out_b[j] = 0.0f;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            out_p[i * 16 + j] = p[j * n + i];

    for (int j = 0; j < m; j++)
        out_b[j] = p[n * m + j];
}

// Hard-coded implementation of lr_gcn
float *evaluate_model(const graph *g, const float *params)
{
    long long maxc = g->N > g->V[g->N] ? g->N : g->V[g->N];
    maxc += 4;

    float *t1 = aligned_alloc(32, sizeof(float) * 32ll * maxc);
    float *t2 = aligned_alloc(32, sizeof(float) * 32ll * maxc);
    float *t3 = aligned_alloc(32, sizeof(float) * 40ll * (g->N + 4));
    float *p = aligned_alloc(32, sizeof(float) * 40ll * 16);
    float *b = aligned_alloc(32, sizeof(float) * 16);

    for (int i = 0; i < (g->N + 4) * 40; i++)
        t3[i] = 0.0f;

    for (int i = 0; i < (g->N + 4) * 16; i++)
        t1[i] = 0.0f;

    // Compute input features
    compute_features(g, t1);

    // Store input in t3 for output layers
    for (int i = 0; i < g->N; i++)
        for (int j = 0; j < 8; j++)
            t3[i * 40 + j] = t1[i * 8 + j];

    // Copy and transpose first layer parameters and run lr_conv
    prep_params(params, 16, 16, p, b);
    lr_conv(g, 8, t1, t2, p, b);

    // Store first layer output in t3 for output layers
    for (int i = 0; i < g->N; i++)
        for (int j = 0; j < 16; j++)
            t3[i * 40 + 8 + j] = t2[i * 16 + j];

    // Copy and transpose second layer parameters and run lr_conv
    prep_params(params + (16 * 16 + 16), 32, 16, p, b);
    lr_conv(g, 16, t2, t1, p, b);

    // Store second layer output in t3 for output layers
    for (int i = 0; i < g->N; i++)
        for (int j = 0; j < 16; j++)
            t3[i * 40 + 24 + j] = t1[i * 16 + j];

    // Copy and transpose third layer parameters and run linear
    prep_params(params + (16 * 16 + 16) + (32 * 16 + 16), 40, 16, p, b);
    matrix_multiplication_ref(t3, p, t2, b, g->N, 40);
    ReLU(t2, g->N);

    // Copy and transpose forth layer parameters and run linear
    prep_params(params + (16 * 16 + 16) + (32 * 16 + 16) + (40 * 16 + 16), 16, 1, p, b);
    matrix_multiplication_ref(t2, p, t1, b, g->N, 16);
    Sigmoid(t1, g->N);

    float *res = malloc(sizeof(float) * g->N);
    for (int i = 0; i < g->N; i++)
        res[i] = t1[i * 16];

    free(t1);
    free(t2);
    free(t3);
    free(p);
    free(b);

    return res;
}

int main(int argc, char **argv)
{
    FILE *fg = fopen(argv[1], "r");
    FILE *fr = fopen(argv[2], "r");
    graph *g = graph_parse_csv(fg, fr);
    fclose(fg);
    fclose(fr);

    // printf("%d %d\n", g->N, g->V[g->N] / 2);

    float *paramv = NULL;

    FILE *f = fopen(argv[3], "r");
    int paramc = lr_gcn_parse(f, &paramv);
    fclose(f);

    // printf("%d %f\n", paramc, paramv[0]);

    float *y = evaluate_model(g, paramv);

    int t = 0;
    for (int i = 0; i < g->N; i++)
        if (y[i] > 0.5f)
            t++;

    printf("%d\n", t);

    graph_free(g);
    free(paramv);
    free(y);

    return 0;
}