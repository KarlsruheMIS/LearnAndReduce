/******************************************************************************
 * graph_access.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_ACCESS_MWIS
#define GRAPH_ACCESS_MWIS

#include <bitset>
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

#include "definitions.h"

struct Node {
    EdgeID firstEdge;
    NodeWeight weight;
    PartitionID partitionID;
};

struct Edge {
    NodeID target;
    EdgeWeight weight;
};

class graph_access;

//construction etc. is encapsulated in basicGraph / access to properties etc. is encapsulated in graph_access
class basicGraph {
    friend class graph_access;

public:
    basicGraph() : m_building_graph(false) {
    }

private:
    //methods only to be used by friend class
    EdgeID number_of_edges() {
        return m_edges.size();
    }

    NodeID number_of_nodes() {
        return m_nodes.size()-1;
    }

    inline EdgeID get_first_edge(const NodeID & node) {
        return m_nodes[node].firstEdge;
    }

    inline EdgeID get_first_invalid_edge(const NodeID & node) {
        return m_nodes[node+1].firstEdge;
    }

    // construction of the graph
    void start_construction(NodeID n, EdgeID m) {
        m_building_graph = true;
        node             = 0;
        e                = 0;
        m_last_source    = -1;

        //resizes property arrays
        m_nodes.resize(n+1);
        m_edges.resize(m);

        m_nodes[node].firstEdge = e;
    }

    // Add a new edge from node 'source' to node 'target'.
    // If an edge with source = n has been added, adding
    // edges with source < n will lead to a broken graph.
    EdgeID new_edge(NodeID source, NodeID target) {
        assert(m_building_graph);
        assert(e < m_edges.size());
       
        m_edges[e].target = target;
        EdgeID e_bar = e;
        ++e;

        assert(source+1 < m_nodes.size());
        m_nodes[source+1].firstEdge = e;

        //fill isolated sources at the end
        if ((NodeID)(m_last_source+1) < source) {
            for (NodeID i = source; i>(NodeID)(m_last_source+1); i--) {
                m_nodes[i].firstEdge = m_nodes[m_last_source+1].firstEdge;
            }
        }
        m_last_source = source;
        return e_bar;
    }

    NodeID new_node() {
        assert(m_building_graph);
        return node++;
    }

    void finish_construction() {
        // insert dummy node
        m_nodes.resize(node+1);
        m_edges.resize(e);

        m_building_graph = false;

        //fill isolated sources at the end
        if ((unsigned int)(m_last_source) != node-1) {
                //in that case at least the last node was an isolated node
                for (NodeID i = node; i>(unsigned int)(m_last_source+1); i--) {
                        m_nodes[i].firstEdge = m_nodes[m_last_source+1].firstEdge;
                }
        }
    }

    // %%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    
    // construction properties
    bool m_building_graph;
    int m_last_source;
    NodeID node; //current node that is constructed
    EdgeID e;    //current edge that is constructed
};

//makros - graph access
#define forall_edges(G,e) { for(EdgeID e = 0, end = G.number_of_edges(); e < end; ++e) {
#define forall_nodes(G,n) { for(NodeID n = 0, end = G.number_of_nodes(); n < end; ++n) {
#define forall_out_edges(G,e,n) { for(EdgeID e = G.get_first_edge(n), end = G.get_first_invalid_edge(n); e < end; ++e) {
#define forall_out_edges_starting_at(G,e,n,e_bar) { for(EdgeID e = e_bar, end = G.get_first_invalid_edge(n); e < end; ++e) {
#define forall_blocks(G,p) { for (PartitionID p = 0, end = G.get_partition_count(); p < end; p++) {
#define endfor }}



class graph_access {
        public:
                graph_access() { 
                    m_max_degree_computed = false; m_max_degree = 0; m_max_weight_computed = false; m_max_weight = 0; graphref = new basicGraph(); }
                virtual ~graph_access(){ delete graphref; };

                graph_access(const graph_access&) = delete;

                /* ============================================================= */
                /* build methods */
                /* ============================================================= */
                void start_construction(NodeID nodes, EdgeID edges);
                NodeID new_node();
                EdgeID new_edge(NodeID source, NodeID target);
                void finish_construction();

                /* ============================================================= */
                /* graph access methods */
                /* ============================================================= */
                NodeID number_of_nodes();
                EdgeID number_of_edges();

                EdgeID get_first_edge(NodeID node);
                EdgeID get_first_invalid_edge(NodeID node);

                PartitionID getPartitionIndex(NodeID node);
                void setPartitionIndex(NodeID node, PartitionID partition);

                NodeWeight getNodeWeight(NodeID node);
                void setNodeWeight(NodeID node, NodeWeight weight);

                EdgeWeight getNodeDegree(NodeID node);
                EdgeWeight getWeightedNodeDegree(NodeID node);
                EdgeWeight getMaxDegree();
                EdgeWeight getMaxWeight();

                EdgeWeight getEdgeWeight(EdgeID edge);
                void setEdgeWeight(EdgeID edge, EdgeWeight weight);

                NodeID getEdgeTarget(EdgeID edge);
                void setEdgeTarget(EdgeID edge, NodeID target);

                int* UNSAFE_metis_style_xadj_array();
                int* UNSAFE_metis_style_adjncy_array();

                int* UNSAFE_metis_style_vwgt_array();
                int* UNSAFE_metis_style_adjwgt_array();

                int build_from_metis(NodeID n, EdgeID* xadj, NodeID* adjncy);
                int build_from_metis_weighted(NodeID n, EdgeID* xadj, NodeID* adjncy, NodeWeight* vwgt, EdgeWeight* adjwgt);

                void copy(graph_access & Gcopy);
        protected:
                basicGraph* graphref;     
                bool         m_max_degree_computed;
                bool         m_max_weight_computed;
                EdgeWeight   m_max_degree;
                EdgeWeight   m_max_weight;
};


/* graph build methods */
inline void graph_access::start_construction(NodeID nodes, EdgeID edges) {
        graphref->start_construction(nodes, edges);
}

inline NodeID graph_access::new_node() {
        return graphref->new_node();
}

inline EdgeID graph_access::new_edge(NodeID source, NodeID target) {
        return graphref->new_edge(source, target);
}

inline void graph_access::finish_construction() {
        graphref->finish_construction();
}

/* graph access methods */
inline NodeID graph_access::number_of_nodes() {
        return graphref->number_of_nodes();
}

inline EdgeID graph_access::number_of_edges() {
        return graphref->number_of_edges();
}

inline EdgeID graph_access::get_first_edge(NodeID node) {
#ifdef NDEBUG
        return graphref->m_nodes[node].firstEdge;
#else
        return graphref->m_nodes.at(node).firstEdge;
#endif
}

inline EdgeID graph_access::get_first_invalid_edge(NodeID node) {
        return graphref->m_nodes[node+1].firstEdge;
}

inline void graph_access::setPartitionIndex(NodeID node, PartitionID partitionID) {
#ifdef NDEBUG
        graphref->m_nodes[node].partitionID = partitionID;
#else
        graphref->m_nodes.at(node).partitionID = partitionID;
#endif
}

inline PartitionID graph_access::getPartitionIndex(NodeID node) {
#ifdef NDEBUG
        return graphref->m_nodes[node].partitionID;
#else
        return graphref->m_nodes.at(node).partitionID;
#endif
}

inline NodeWeight graph_access::getNodeWeight(NodeID node){
#ifdef NDEBUG
        return graphref->m_nodes[node].weight;        
#else
        return graphref->m_nodes.at(node).weight;        
#endif
}

inline void graph_access::setNodeWeight(NodeID node, NodeWeight weight){
#ifdef NDEBUG
        graphref->m_nodes[node].weight = weight;        
#else
        graphref->m_nodes.at(node).weight = weight;        
#endif
}

inline EdgeWeight graph_access::getEdgeWeight(EdgeID edge){
#ifdef NDEBUG
        return graphref->m_edges[edge].weight;        
#else
        return graphref->m_edges.at(edge).weight;        
#endif
}

inline void graph_access::setEdgeWeight(EdgeID edge, EdgeWeight weight){
#ifdef NDEBUG
        graphref->m_edges[edge].weight = weight;        
#else
        graphref->m_edges.at(edge).weight = weight;        
#endif
}

inline NodeID graph_access::getEdgeTarget(EdgeID edge){
#ifdef NDEBUG
        return graphref->m_edges[edge].target;        
#else
        return graphref->m_edges.at(edge).target;        
#endif
}

inline void graph_access::setEdgeTarget(EdgeID edge, NodeID target) {
#ifdef NDEBUG
    graphref->m_edges[edge].target = target;
#else
    graphref->m_edges.at(edge).target = target;
#endif
}

inline EdgeWeight graph_access::getNodeDegree(NodeID node) {
        return graphref->m_nodes[node+1].firstEdge-graphref->m_nodes[node].firstEdge;
}

inline EdgeWeight graph_access::getWeightedNodeDegree(NodeID node) {
	EdgeWeight degree = 0;
	for( unsigned e = graphref->m_nodes[node].firstEdge; e < graphref->m_nodes[node+1].firstEdge; ++e) {
		degree += getEdgeWeight(e);
	}
        return degree;
}

inline EdgeWeight graph_access::getMaxWeight() {
        if(!m_max_weight_computed) {
                //compute it
                basicGraph& ref = *graphref;
                forall_nodes(ref, node) {
                        if(getNodeWeight(node) > static_cast<NodeWeight>(m_max_weight)) {
                                m_max_weight = getNodeWeight(node);
                        }
                } endfor
                m_max_weight_computed = true;
        }

        return m_max_weight;
}

inline EdgeWeight graph_access::getMaxDegree() {
        if(!m_max_degree_computed) {
                //compute it
                basicGraph& ref = *graphref;
                forall_nodes(ref, node) {
                        EdgeWeight cur_degree = 0;
                        forall_out_edges(ref, e, node) {
                                cur_degree += getEdgeWeight(e);
                        } endfor
                        if(cur_degree > m_max_degree) {
                                m_max_degree = cur_degree;
                        }
                } endfor
                m_max_degree_computed = true;
        }

        return m_max_degree;
}

inline int* graph_access::UNSAFE_metis_style_xadj_array() {
        int* xadj      = new int[graphref->number_of_nodes()+1];
        basicGraph& ref = *graphref;

        forall_nodes(ref, n) {
                xadj[n] = graphref->m_nodes[n].firstEdge;
        } endfor
        xadj[graphref->number_of_nodes()] = graphref->m_nodes[graphref->number_of_nodes()].firstEdge;
        return xadj;
}


inline int* graph_access::UNSAFE_metis_style_adjncy_array() {
        int* adjncy    = new int[graphref->number_of_edges()];
        basicGraph& ref = *graphref;
        forall_edges(ref, e) {
                adjncy[e] = graphref->m_edges[e].target;
        } endfor 

        return adjncy;
}


inline int* graph_access::UNSAFE_metis_style_vwgt_array() {
        int* vwgt      = new int[graphref->number_of_nodes()];
        basicGraph& ref = *graphref;

        forall_nodes(ref, n) {
                vwgt[n] = (int)graphref->m_nodes[n].weight;
        } endfor
        return vwgt;
}

inline int* graph_access::UNSAFE_metis_style_adjwgt_array() {
        int* adjwgt    = new int[graphref->number_of_edges()];
        basicGraph& ref = *graphref;

        forall_edges(ref, e) {
                adjwgt[e] = (int)graphref->m_edges[e].weight;
        } endfor 

        return adjwgt;
}

inline int graph_access::build_from_metis(NodeID n, EdgeID* xadj, NodeID* adjncy) {
        if (graphref != nullptr) {
            delete graphref;
        }
        graphref = new basicGraph();
        start_construction(n, xadj[n]);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, 1);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, 1);
                }

        }
        
        finish_construction();
        return 0;
}

inline int graph_access::build_from_metis_weighted(NodeID n, EdgeID* xadj, NodeID* adjncy, NodeWeight* vwgt, EdgeWeight* adjwgt) {
        if (graphref != nullptr) {
            delete graphref;
        }
        graphref = new basicGraph();
        start_construction(n, xadj[n]);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, vwgt[i]);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, adjwgt[e]);
                }
        }
        
        finish_construction();
        return 0;
}

inline void graph_access::copy(graph_access & G_bar) {
        G_bar.start_construction(number_of_nodes(), number_of_edges());

        basicGraph& ref = *graphref;
        forall_nodes(ref, node) {
                NodeID shadow_node = G_bar.new_node();
                G_bar.setNodeWeight(shadow_node, getNodeWeight(node));
                forall_out_edges(ref, e, node) {
                        NodeID target                   = getEdgeTarget(e);
                        EdgeID shadow_edge              = G_bar.new_edge(shadow_node, target);
                        G_bar.setEdgeWeight(shadow_edge, getEdgeWeight(e));
                } endfor
        } endfor

        G_bar.finish_construction();
}

#endif /* end of include guard: GRAPH_ACCESS_EFRXO4X2 */
