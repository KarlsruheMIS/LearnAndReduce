#pragma once

#include "graph_access.h"
// #include "dynamic_partition_graph.h"
// #include "parallel_recursive_reduce_algorithm.h"
#include <unordered_set>

template <class Graph>
class solution_check {
    Graph& G;

    public:
    solution_check<Graph>(Graph& graph);
    ~solution_check();

    bool check_overall_solution(NodeWeight weight);
    bool check_IS();
    bool check_weight(NodeWeight weight);
    bool is_free(NodeID node);
    bool is_maximal();
    NodeWeight make_maximal();
    // bool check_dyn_graph(dynamic_partition_graph& G, std::vector<parallel_recursive_reduce_algorithm::IS_status>& status);
    bool check_graph();
};

template <class Graph>
solution_check<Graph>::solution_check(Graph& graph) : G(graph){ }

template <class Graph>
solution_check<Graph>::~solution_check() { }


template <class Graph>
bool solution_check<Graph>::check_overall_solution(NodeWeight weight) {
    if (!check_IS()) {
        return false;
    }
    if (!check_weight(weight)) {
        return false;
    }
    if (!is_maximal()) {
        return false;
    }
    return true;
}

template <class Graph>
bool solution_check<Graph>::check_weight(NodeWeight weight) { 
  ::NodeWeight computed_w = 0;
  forall_nodes(G, node) { 
      if (G.getPartitionIndex(node)==1) computed_w += G.getNodeWeight(node);
  } endfor 

  if (computed_w != weight) {
      std::cout << "ERROR: not correct weight of IS\nrecomputed weight: " << computed_w << "=!="<< weight << "(given weight)"<< std::endl;
      return false;
  }
  return true;
}

template <class Graph>
bool solution_check<Graph>::check_IS() { 
	forall_nodes(G, node) {
		if (G.getPartitionIndex(node) == 1) {
			forall_out_edges(G, edge, node) {
				NodeID neighbor = G.getEdgeTarget(edge);
				if (G.getPartitionIndex(neighbor) == 1) {
                    std::cout << "ERROR: Not an independent set!" << std::endl;;
                    std::cout << " nodes " << neighbor << " and " << node << " both included" << std::endl;
					return false;
				}
			} endfor
		}
	} endfor

    return true;
}

template <class Graph>
bool solution_check<Graph>::is_free(NodeID node) {
    //node itself is in solution
    if (G.getPartitionIndex(node) == 1) {
        return false;
    }

    //check if neighbor node is in solution
    forall_out_edges(G, edge, node) {
        NodeID neighbor = G.getEdgeTarget(edge);
        if (G.getPartitionIndex(neighbor) == 1) {
            return false;
        }
    } endfor
    return true;
}

template <class Graph>
bool solution_check<Graph>::is_maximal() {
    forall_nodes(G, node) {
        if (is_free(node)) {
            std::cout << "ERROR: Independent set not maximal!" << std::endl;
            return false;
        }
    } endfor
    return true;
}

template <class Graph>
NodeWeight solution_check<Graph>::make_maximal() {
    NodeWeight weight = 0;
    forall_nodes(G, node) {
        if (is_free(node)) {
            G.setPartitionIndex(node, 1);
            weight += G.getNodeWeight(node);
        }
    } endfor
    return weight;
}

// template <class Graph>
// bool solution_check<Graph>::check_dyn_graph(dynamic_partition_graph& G, std::vector<parallel_recursive_reduce_algorithm::IS_status>& status)
// {
//     for (NodeID node = 0; node < G.size(); node++) {
//         if (!(status[node] == IS_status::not_set || 
//              status[node] == IS_status::boundary)) {
//             continue;
//         }
//         for (NodeID neighbor : G[node])
//         {
//             if (neighbor == node) {
//                 std::cout << "The graph contains a self-loop." << std::endl;
//                 std::cout << "Node " << node << " has an edge to itself." << std::endl;
//                 return false;
//             }
//             bool directed_edge_exists = false;
//             for (NodeID back_neighbor : G[neighbor]) {
//                 if (back_neighbor == node) {
//                     directed_edge_exists = true;
//                 }
//             }
//             if (!directed_edge_exists) 
//             {
//                 std::cout << "There is a edge with no backwards edge." << std::endl;
//                 std::cout << "The edge " << neighbor << "--" << node << " is missing." << std::endl;
//                 return false;
//             }
//         }
//     }
//     return true;
// }

template <class Graph>
bool solution_check<Graph>::check_graph() {
    long n = G.number_of_nodes();
    long m = G.number_of_edges() / 2;

    std::vector<long> node_starts(n + 1, 0);
    std::vector<long> adjacent_nodes;
    adjacent_nodes.reserve(m * 2);

    long edge_counter = 0;
    long node_counter = 0;

    forall_nodes(G, node) {
        long node_weight = G.getNodeWeight(node);
        if (node_weight < 0) {
            std::cout << "Node " << node + 1 << " has weight < 0." << std::endl;
            return false;
        }

        long node_degree = 0;
        forall_out_edges(G, e, node) {
            long target = G.getEdgeTarget(e);
            if (target >= n || target < 0) {
                std::cout << "Node " << node + 1 << " has an edge to an invalid node " << target + 1 << "." << std::endl;
                return false;
            }
            adjacent_nodes.push_back(target);
            node_degree++;
            edge_counter++;
        } endfor
        node_starts[node + 1] = node_starts[node] + node_degree;
        node_counter++;
    } endfor

    if (edge_counter != m * 2) {
        std::cout << "Edge count mismatch. Expected " << m * 2 << ", found " << edge_counter << "." << std::endl;
        return false;
    }
    
    if (node_counter != n) {
        std::cout << "Node count mismatch. Expected " << n<< ", found " << node_counter << "." << std::endl;
        return false;
    }

    // Check for parallel edges and self-loops
    for (long node = 0; node < n; node++) {
        std::unordered_set<long> seen;
        for (long e = node_starts[node]; e < node_starts[node + 1]; e++) {
            long target = adjacent_nodes[e];
            if (!seen.insert(target).second) {
                std::cout << "Parallel edge or self-loop found at node " << node + 1 << "." << std::endl;
                return false;
            }
            if (target == node) {
                std::cout << "Self-loop detected at node " << node + 1 << "." << std::endl;
                return false;
            }
        }
    }


    // Checking for all forward and backward edges
    for (long node = 0; node < n; node++) {
        for (long e = node_starts[node]; e < node_starts[node + 1]; e++) {
            long target = adjacent_nodes[e];
            bool found = false;
            for (long e_bar = node_starts[target]; e_bar < node_starts[target + 1]; e_bar++) {
                if (adjacent_nodes[e_bar] == node) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "The graph does not contain all forward and backward edges. " << std::endl;
                std::cout << "Node " << node + 1 << " has an edge to node " << target + 1 
                          << ", but there is no edge (" << target + 1 << "," << node + 1 
                          << ") present." << std::endl;
                return false;
            }
        }
    }

    return true; // Graph passed all checks
}

