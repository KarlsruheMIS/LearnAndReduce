#include "bounds.h"

#include <vector>
#include <algorithm>

void greedy_coloring_weighted(int N, std::vector<int> &V, std::vector<int> &E,
                              std::vector<int> &W, std::vector<int> &WC,
                              std::vector<int> &p, std::vector<int> &color)
{
    int n_next = 0, n_prev = 0, max_degree = 0, K = 0;
    std::vector<int> E_split(E.size()), count(N), pred(N), prev(N), next(N), marks(N + 1);
    for (int u = 0; u < N; u++)
    {
        int j = V[u], k = V[u + 1] - 1;

        for (int i = V[u]; i < V[u + 1]; i++)
        {
            int v = E[i];
            if (p[v] > p[u] || (p[u] == p[v] && u > v))
                E_split[j++] = v;
            else
                E_split[k--] = v;
        }

        color[u] = 0;
        count[u] = j - V[u];
        pred[u] = j;

        if (count[u] == 0)
            prev[n_prev++] = u;

        if (V[u + 1] - V[u] > max_degree)
            max_degree = V[u + 1] - V[u];
    }

    for (int i = 0; i < max_degree + 1; i++)
        marks[i] = N;

    while (n_prev > 0)
    {
        for (int i = 0; i < n_prev; i++)
        {
            int u = prev[i];

            // Find lowest available color
            for (int j = V[u]; j < pred[u]; j++)
                marks[color[E_split[j]]] = u;

            int c = 1;
            while (c < N && (marks[c] == u || (marks[c + 1] != u && WC[c + 1] >= W[u])))
                c++;
            color[u] = c;

            if (c > K)
            {
                K = c;
                WC[K] = W[u];
            }

            // Updete uncolored neighbors
            for (int j = pred[u]; j < V[u + 1]; j++)
            {
                int v = E_split[j];

                count[v]--;
                if (count[v] == 0)
                    next[n_next++] = v;
            }
        }

        n_prev = n_next;
        n_next = 0;

        std::swap(next, prev);
    }
}

int greedy_ub(graph_access &g)
{
    int N = g.number_of_nodes();
    std::vector<int> V(N + 1), E, p(N), M(N + 1, -1), W(N), CW(N + 1);
    V[0] = 0;
    for (int u = 0; u < N; u++)
    {
        M[u] = u;
        W[u] = g.getNodeWeight(u);
        CW[u] = 0;
        for (int e = g.get_first_edge(u); e != g.get_first_invalid_edge(u); e++)
            M[g.getEdgeTarget(e)] = u;
        for (int v = 0; v < N; v++)
            if (M[v] != u)
                E.push_back(v);
        V[u + 1] = E.size();
    }
    CW[N] = 0;

    for (int u = 0; u < N; u++)
        p[u] = W[u];

    std::vector<int> color(N);
    greedy_coloring_weighted(N, V, E, W, CW, p, color);

    int ub = 0;
    for (int i = 0; i < N; i++)
        M[i] = 0;

    int max_c = 1;
    for (int u = 0; u < N; u++)
    {
        if (color[u] > max_c)
            max_c = color[u];
        if (g.getNodeWeight(u) > M[color[u]])
            M[color[u]] = g.getNodeWeight(u);
    }
    // printf("%d\n", max_c);
    for (int i = 1; i < N; i++)
        ub += M[i];

    return ub;
}

int greedy_lb(graph_access &g)
{
    std::vector<int> used(g.number_of_nodes(), 0);
    std::vector<std::pair<int, int>> node_weight(g.number_of_nodes());
    for (int u = 0; u < g.number_of_nodes(); u++)
        node_weight[u] = {g.getNodeWeight(u), u};
    std::sort(node_weight.begin(), node_weight.end(), [](auto a, auto b)
              { return a.first > b.first; });

    int lb = 0;
    for (int i = 0; i < g.number_of_nodes(); i++)
    {
        auto [w, u] = node_weight[i];
        if (used[u])
            continue;
        lb += w;
        used[u] = 1;
        for (int e = g.get_first_edge(u); e != g.get_first_invalid_edge(u); e++)
        {
            int v = g.getEdgeTarget(e);
            used[v] = 1;
        }
    }
    return lb;
}