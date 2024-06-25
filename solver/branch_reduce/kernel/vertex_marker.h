#pragma once

#include <definitions.h>
#include "fast_set.h"

class vertex_marker {
private:

public:

    std::vector<NodeID> current;
    std::vector<NodeID> next;
    fast_set added_vertices;

    vertex_marker() : vertex_marker(10) {}
	vertex_marker(size_t size) : added_vertices(size) {
		current.reserve(size); 
		next.reserve(size); 
	}

	void add(NodeID vertex) {
		if (!added_vertices.get(vertex)) {
			next.push_back(vertex);
			added_vertices.add(vertex);
		}
	}

	void add(NodeID vertex, std::vector<NodeID>& reverse_mapping) {
		if (next.size() == next.capacity()) {
			std::cout << "current.size() == current.capacity()" << std::endl;
			assert(false);
		}
		if (!added_vertices.get(reverse_mapping[vertex])) {
			next.push_back(vertex);
			assert(reverse_mapping[vertex] < current.capacity());
			added_vertices.add(reverse_mapping[vertex]);
		}
	}

	void get_next() {
			current.swap(next);
			clear_next();
	}

	void clear_next() {
		next.clear();
		added_vertices.clear();
	}

	void fill_current_partition_ascending(size_t n, std::vector<NodeID>& reverse_mapping) {
		current.clear();
		for (size_t i = 0; i < n; i++) {
			current.push_back(reverse_mapping[i]);
		} 
	}

	void fill_current_ascending(size_t n) {
		current.clear();
		for (size_t i = 0; i < n; i++) {
			current.push_back(i);
		}
	}

	NodeID current_vertex(size_t index) {
		return current[index];
	}

	size_t current_size() {
		return current.size();
	}

    void resize(size_t size) {
        added_vertices.resize(size);
    }

    std::vector<NodeID> &current_vec() {
	    return current;
	}
};

