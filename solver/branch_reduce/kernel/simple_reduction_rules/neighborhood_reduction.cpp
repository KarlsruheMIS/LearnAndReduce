#include "neighborhood_reduction.h"
#include "general_reduction.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    size_t oldn = br_alg->status.remaining_nodes;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    for_each_changed_vertex(br_alg, [&br_alg, &status](NodeID v)
                            {
                                NodeWeight neighbor_weights = 0;
                                for (NodeID u : status.graph[v])
                                {
                                    neighbor_weights += status.weights[u];
                                    if (status.weights[v] < neighbor_weights)
                                    {
                                        return;
                                    }
                                }
                                br_alg->set(v, IS_status::included);
                            });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != status.remaining_nodes) std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;
    return oldn != br_alg->status.remaining_nodes;
}

inline bool neighborhood_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    // if (br_alg->config.disable_neighborhood) return false;
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (status.weights[v] >= neighbor_weights)
    {
        br_alg->set(v, IS_status::included);
    }
    return oldn != br_alg->status.remaining_nodes;
}
bool neighborhood_reduction::is_suited(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    return (status.weights[v] >= static_cast<NodeWeight>(std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&](NodeWeight sum, NodeID neighbor)
                                                                         { 
			if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor]; 
			else return sum; })));
}

inline int neighborhood_reduction::generate_data(branch_and_reduce_algorithm *br_alg, NodeID v, std::vector<NodeID>& label)
{
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    NodeWeight neighbor_weights = get_neighborhood_weight(v, br_alg);
    if (status.weights[v] >= neighbor_weights)
    {
        label.push_back(v);
    }
    return label.size() > 0;
}