#pragma once

#include "graph_access.h"
#include "fast_set.h"
#include "sized_vector.h"
#include "config.h"
#include <random>

class graph_operations {

    public:
    graph_operations() {};
    ~graph_operations() {};

    template <typename Graph, typename vec>
    void induce_subgraph(Graph &G, Graph &subG, const vec &nodes, vec &reverse_mapping);

    template <typename Graph, typename C>
    void assign_weights(Graph &G, const C &config);

    template <typename Graph>
    bool is_free(Graph &G, NodeID node);

    template <typename Graph>
    bool check_graph(Graph &G);

    template <typename Graph>
    void print_G(Graph &G);

    template <typename Graph>
    int writeGraphWeighted(Graph &G, const std::string& filename);

    template <typename Graph>
    int writeGraphWeighted_to_csv(Graph &G, const std::string& filename);

};
template <typename Graph, typename vec>
void graph_operations::induce_subgraph(Graph& G, Graph& sub_G, const vec& nodes, vec& reverse_mapping) {

	assert(nodes.size() == reverse_mapping.size() && "induce_subgraph: nodes and reverse_mapping must have the same size");
    std::vector<NodeID> mapping(G.number_of_nodes(), -1);

    //fill neighbor set
    fast_set set(G.number_of_nodes());
	for (size_t i = 0; i < nodes.size(); i++) {
        set.add(nodes[i]);
    }

    // count edges and set mappings
	size_t edge_count = 0;
	for (size_t i = 0; i < nodes.size(); i++) {
		NodeID node = nodes[i];
        for (EdgeID e = G.get_first_edge(node); e < G.get_first_invalid_edge(node); ++e) {
            NodeID neighbor = G.getEdgeTarget(e);
            if (set.get(neighbor))
				edge_count++;
		} 
		reverse_mapping[i] = node;
		mapping[node] = i;
	}

	sub_G.start_construction(nodes.size(), edge_count);

	// add nodes and edges
	for (size_t i = 0; i < nodes.size(); i++) {
		NodeID node = nodes[i];
		NodeID new_node = sub_G.new_node();

		sub_G.setNodeWeight(new_node, G.getNodeWeight(node));

        for (EdgeID e = G.get_first_edge(node); e < G.get_first_invalid_edge(node); ++e) {
            NodeID neighbor = G.getEdgeTarget(e);
			if (set.get(neighbor)) {
				sub_G.new_edge(new_node, mapping[neighbor]);
            }
		}
	}

	sub_G.finish_construction();
}

template <typename Graph>
void graph_operations::print_G(Graph& G) {
  forall_nodes(G, node) { std::cout << node << " weight("<<G.getNodeWeight(node)<<")" <<" partition("<<G.getPartitionIndex(node)<<"): ";
    forall_out_edges(G, e, node) {
        NodeID target = G.getEdgeTarget(e);
        std::cout << target <<", ";
    } endfor std::cout<<std::endl; } endfor
}

template <typename Graph, typename C>
void graph_operations::assign_weights(Graph& G, const C& mis_config) {
	constexpr NodeWeight MAX_WEIGHT = 200;

	if (mis_config.weight_source == Config::Weight_Source::HYBRID) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, (node + 1) % MAX_WEIGHT + 1);
		} endfor
	} else if (mis_config.weight_source == Config::Weight_Source::UNIFORM) {
		std::default_random_engine generator(mis_config.seed);
  		std::uniform_int_distribution<NodeWeight> distribution(1,MAX_WEIGHT);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	} else if (mis_config.weight_source == Config::Weight_Source::GEOMETRIC) {
		std::default_random_engine generator(mis_config.seed);
  		std::binomial_distribution<int> distribution(MAX_WEIGHT / 2);

		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	} else if (mis_config.weight_source == Config::Weight_Source::UNIT) {
		forall_nodes(G, node) {
			G.setNodeWeight(node, 1);
		} endfor
	} else if (mis_config.weight_source == Config::Weight_Source::SMALL_UNIFORM) {
		std::default_random_engine generator(mis_config.seed);
  		std::uniform_int_distribution<NodeWeight> distribution(2,20);
		forall_nodes(G, node) {
			G.setNodeWeight(node, distribution(generator));
		} endfor
	}
}

template<typename Graph>
bool graph_operations::is_free(Graph &G, NodeID node) {
    forall_out_edges(G, e, node) {
        NodeID target = G.getEdgeTarget(e);

        // if target in independent set, node is not
        if (G.getPartitionIndex(target) == 1) {
            return false;
        }
    } endfor
    return true;
}

template<typename Graph>
int graph_operations::writeGraphWeighted(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 <<  " 10" <<  std::endl;

        forall_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1);
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

template<typename Graph>
int graph_operations::writeGraphWeighted_to_csv(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 <<  " 10" <<  std::endl;

        forall_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1);
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

