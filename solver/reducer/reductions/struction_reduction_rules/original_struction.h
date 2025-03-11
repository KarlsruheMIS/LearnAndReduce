//
// Created by alex on 16.01.20.
//

#ifndef COMPONENTS_ORIGINAL_H
#define COMPONENTS_ORIGINAL_H
class reduce_algorithm;
template<bool modified = false>
class original_struction {

public:
    bool reduce(reduce_algorithm* r_alg, NodeID n, size_t max_nodes);
    void apply(reduce_algorithm* r_alg) {}
    void restore(reduce_algorithm* r_alg) {}
    size_t removed_vertices(reduce_algorithm* r_alg, NodeID n) {return 1;}

private:
    reduce_algorithm* r_alg;
    struct restore_data {
        explicit restore_data(NodeID main, size_t set_nodes) : main(main), set_nodes(set_nodes) {}
        NodeID main;
        size_t set_nodes;
    };
    std::vector<restore_data> restore_vec;

    bool find_additional_nodes(std::vector<NodeID> &additional_nodes, neighbor_list &neighbors, std::vector<NodeID> &layer_prefix_count, size_t max_nodes) const;
};

#endif //COMPONENTS_ORIGINAL_H
