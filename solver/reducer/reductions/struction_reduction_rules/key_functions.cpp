//
// Created by alex on 31.01.20.
//
#include "key_functions.h"
#include "reduce_algorithm.h"


Gain ApproximateIncreaseKey::key(reduce_algorithm* r_alg, NodeID n) {
    auto &status = r_alg->status;
    auto &neighbors = status.graph[n];
    auto &set = r_alg->set_1;

    size_t new_nodes_estimate = 0;
    for (NodeID u : status.graph[n]) {
        new_nodes_estimate += status.weights[n] < status.weights[u];
    }
    set.resize(status.graph.size());
    for (NodeID u : neighbors) {
        set.clear();
        for (NodeID v : status.graph[u]) {
            set.add(v);
        }
        for (NodeID v : neighbors) {
            new_nodes_estimate += u > v && !set.get(v) && status.weights[n] < status.weights[u] + status.weights[v];
        }
    }

    return key_by_set_estimate(r_alg, n, new_nodes_estimate);
}

Gain ApproximateIncreaseKey::key_by_set_estimate(reduce_algorithm* r_alg, NodeID n, size_t set_estimate) {
    return (int) r_alg->deg(n) - (int) set_estimate;
}

size_t ApproximateIncreaseKey::set_limit(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return key_reinsert_factor * set_estimate(r_alg, n, key) + 1;
}

size_t ApproximateIncreaseKey::set_estimate(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return r_alg->deg(n) - key;
}



Gain DegreeKey::key(reduce_algorithm* r_alg, NodeID n) {
    return -r_alg->deg(n);
}

Gain DegreeKey::key_by_set_estimate(reduce_algorithm* r_alg, NodeID n, size_t set_estimate) {
    return std::numeric_limits<Gain>::min();
}

size_t DegreeKey::set_limit(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return std::numeric_limits<size_t>::max() - 1;
}

size_t DegreeKey::set_estimate(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return 0;
}


Gain IncreaseKey::key(reduce_algorithm* r_alg, NodeID n) {
    auto &status = r_alg->status;
    auto &neighbors = status.graph[n];
    finder.findAllMWIS<false>(r_alg, neighbors, status.weights[n], 1000);
    size_t new_nodes_estimate = finder.get_sets().size();

    return key_by_set_estimate(r_alg, n, new_nodes_estimate);
}

Gain IncreaseKey::key_by_set_estimate(reduce_algorithm* r_alg, NodeID n, size_t set_estimate) {
    return (int) r_alg->deg(n) - (int) set_estimate;
}

size_t IncreaseKey::set_limit(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return set_estimate(r_alg, n, key) + 1;
}

size_t IncreaseKey::set_estimate(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return r_alg->deg(n) - key;
}


Gain RandomKey::key(reduce_algorithm* r_alg, NodeID n) {
    return distribution(generator);
}

Gain RandomKey::key_by_set_estimate(reduce_algorithm* r_alg, NodeID n, size_t set_estimate) {
    return std::numeric_limits<Gain>::min();
}

size_t RandomKey::set_limit(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return std::numeric_limits<size_t>::max() - 1;
}

size_t RandomKey::set_estimate(reduce_algorithm* r_alg, NodeID n, Gain key) {
    return 0;
}
