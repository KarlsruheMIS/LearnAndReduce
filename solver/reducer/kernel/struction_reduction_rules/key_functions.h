//
// Created by alex on 31.01.20.
//

#ifndef COMPONENTS_KEY_FUNCTIONS_H
#define COMPONENTS_KEY_FUNCTIONS_H

#include <definitions.h>
#include <random>
#include "reduction_config.h"
#include "mwis_finder.h"

class reduce_algorithm;
struct ApproximateIncreaseKey {
    ApproximateIncreaseKey(const ReductionConfig &config) : key_reinsert_factor(config.key_reinsert_factor) {}
    Gain key(reduce_algorithm* br_alg, NodeID n);
    Gain key_by_set_estimate(reduce_algorithm* br_alg, NodeID n, size_t set_estimate);
    size_t set_limit(reduce_algorithm* br_alg, NodeID n, Gain key);
    size_t set_estimate(reduce_algorithm* br_alg, NodeID n, Gain key);

private:
    const double key_reinsert_factor;
};
struct IncreaseKey {
    IncreaseKey(const ReductionConfig &config) {}
    Gain key(reduce_algorithm* br_alg, NodeID n);
    Gain key_by_set_estimate(reduce_algorithm* br_alg, NodeID n, size_t set_estimate);
    size_t set_limit(reduce_algorithm* br_alg, NodeID n, Gain key);
    size_t set_estimate(reduce_algorithm* br_alg, NodeID n, Gain key);
private:
    mwis_finder finder;
};
struct DegreeKey {
    DegreeKey(const ReductionConfig &config) {}
    Gain key(reduce_algorithm* br_alg, NodeID n);
    Gain key_by_set_estimate(reduce_algorithm* br_alg, NodeID n, size_t set_estimate);
    size_t set_limit(reduce_algorithm* br_alg, NodeID n, Gain key);
    size_t set_estimate(reduce_algorithm* br_alg, NodeID n, Gain key);
};
struct RandomKey {
    RandomKey(const ReductionConfig &config) : generator(config.seed), distribution(std::numeric_limits<Gain>::min() + 1, std::numeric_limits<Gain>::max() - 1) {}
    Gain key(reduce_algorithm* br_alg, NodeID n);
    Gain key_by_set_estimate(reduce_algorithm* br_alg, NodeID n, size_t set_estimate);
    size_t set_limit(reduce_algorithm* br_alg, NodeID n, Gain key);
    size_t set_estimate(reduce_algorithm* br_alg, NodeID n, Gain key);
private:
    std::default_random_engine generator;
    std::uniform_int_distribution<Gain> distribution;
};
#endif //COMPONENTS_KEY_FUNCTIONS_H
