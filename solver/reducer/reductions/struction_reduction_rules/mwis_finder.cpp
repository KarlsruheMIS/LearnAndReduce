//
// Created by alex on 14.12.19.
//

#include "mwis_finder.h"
#include "reduce_algorithm.h"

std::vector<mwis>& mwis_finder::get_sets() {
    return independent_sets;
}
template<bool just_greater>
bool mwis_finder::findAllMWIS(reduce_algorithm* r_alg, neighbor_list &nodes, NodeWeight min_weight, size_t max_sets) {
    const auto &status = r_alg->status;
    auto &node_state = r_alg->zero_vec;
    independent_sets.clear();

    auto &set = r_alg->set_1;
    set.clear();
    NodeWeight neighbor_weight = 0;
    for (const NodeID n : nodes) {
        node_state[n] = INCLUDABLE;
        neighbor_weight += status.weights[n];
        set.add(n);
    }
    this->r_alg = r_alg;
    this->min_weight = min_weight;
    this->max_sets = max_sets + 1;

    findAllMWISRec<just_greater>(nodes, 0, neighbor_weight);
    for (const NodeID n : nodes) {
        node_state[n] = NOT_INCLUDABLE;
    }

    return independent_sets.size() <= max_sets;
}

template<bool just_greater>
void mwis_finder::findAllMWISRec(neighbor_list &nodes, size_t start_idx, NodeWeight remaining_neighbor_weight) {
    if (independent_sets.size() == max_sets || (just_greater && tmp_mwis.weight > min_weight) || tmp_mwis.weight + remaining_neighbor_weight < min_weight) {
        return;
    }

    const auto &status = r_alg->status;
    auto &node_state = r_alg->zero_vec;
    /* auto &order = r_alg->buffers[1]; */

    NodeID cur_node;
    while (true) {
        if (start_idx == nodes.size())
            return;
        cur_node = nodes[start_idx];
        if (node_state[cur_node] == INCLUDABLE)
            break;
        ++start_idx;
    }

    tmp_mwis.push_back(cur_node, status.weights[cur_node]);
    if (tmp_mwis.weight > min_weight) {
        independent_sets.push_back(tmp_mwis);
    }

    ++node_state[cur_node];
    remaining_neighbor_weight -= status.weights[cur_node];
    for (NodeID u : status.graph[cur_node]) {
        if (node_state[u] == NOT_INCLUDABLE || node_state[u]++ != INCLUDABLE) continue;

        remaining_neighbor_weight -= status.weights[u];
    }

    //Find all maximum independent sets containing cur_node.node recusively
    findAllMWISRec<just_greater>(nodes, start_idx + 1, remaining_neighbor_weight);
    tmp_mwis.pop_back(status.weights[cur_node]);

    for (NodeID u : status.graph[cur_node]) {
        if (node_state[u] == NOT_INCLUDABLE || --node_state[u] != INCLUDABLE) continue;

        remaining_neighbor_weight += status.weights[u];
    }

    //Find all maximum independent sets not containing cur_node.node recusively
    findAllMWISRec<just_greater>(nodes, start_idx + 1, remaining_neighbor_weight);

    --node_state[cur_node];
}

template bool mwis_finder::findAllMWIS<true>(reduce_algorithm* r_alg, neighbor_list &nodes, NodeWeight min_weight, size_t max_sets);
template bool mwis_finder::findAllMWIS<false>(reduce_algorithm* r_alg, neighbor_list &nodes, NodeWeight min_weight, size_t max_sets);
