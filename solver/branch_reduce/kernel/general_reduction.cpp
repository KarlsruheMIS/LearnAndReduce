
#include "general_reduction.h"

#include "branch_and_reduce_algorithm.h"
#include "bounds.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool general_reduction::is_suited(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    return br_alg->status.node_status[v] == IS_status::not_set;
}
NodeID general_reduction::get_max_weight_neighbor(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    NodeID max_neighbor = status.graph[v][0];
    NodeWeight max_weight = status.weights[max_neighbor];

    for (NodeID neighbor : status.graph[v])
    {
        if (status.weights[neighbor] > max_weight && status.node_status[neighbor] == IS_status::not_set)
        {
            max_weight = status.weights[neighbor];
            max_neighbor = neighbor;
        }
    }

    return max_neighbor;
}
NodeWeight general_reduction::get_neighborhood_weight(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    return std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&status](NodeWeight sum, NodeID neighbor)
                           { 
        if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor];
        else return sum; });
}
void general_reduction::get_neighborhood_set(NodeID v, branch_and_reduce_algorithm *br_alg, fast_set &set)
{
    auto &status = br_alg->status;
    set.clear();
    for (NodeID neighbor : status.graph[v])
    {
        if (status.node_status[neighbor] == IS_status::not_set)
            set.add(neighbor);
    }
}
void general_reduction::get_neighborhood_vector(NodeID v, branch_and_reduce_algorithm *br_alg, std::vector<NodeID> &vec)
{
    auto &status = br_alg->status;
    vec.clear();
    for (NodeID neighbor : status.graph[v])
    {
        if (status.node_status[neighbor] == IS_status::not_set)
            vec.push_back(neighbor);
    }
}
bool general_reduction::try_neighborhood_reduction(NodeID v, branch_and_reduce_algorithm *br_alg, NodeWeight neighbors_weight)
{
#ifdef gen_training_data
    return false;
#endif
    auto &status = br_alg->status;
    if (status.weights[v] >= neighbors_weight)
    {
        br_alg->set(v, IS_status::included);
        return true;
    }
    return false;
}
bool general_reduction::solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution, graph_access &graph, branch_and_reduce_algorithm *br_alg, std::vector<NodeID> &nodes_vec, const fast_set &nodes_set, std::vector<NodeID> &reverse_mapping, bool apply_solution)
{
    if (nodes_vec.size() == 0)
    {
        solution = 0;
        return true;
    }
    if (nodes_vec.size() == 1)
    {
        solution = br_alg->status.weights[nodes_vec[0]];
        return true;
    }
    reverse_mapping.resize(br_alg->status.n);
    br_alg->build_induced_subgraph(graph, nodes_vec, nodes_set, reverse_mapping);
    return solve_graph(solution, graph, br_alg->config, weight_bound, apply_solution);
}
bool general_reduction::solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, graph_access &neighborhood_graph, branch_and_reduce_algorithm *br_alg, NodeID v, bool apply_solution)
{
    br_alg->build_induced_neighborhood_subgraph(neighborhood_graph, v);
    return solve_graph(solution, neighborhood_graph, br_alg->config, weight_bound, apply_solution);
}

bool general_reduction::solve_graph(NodeWeight &solution, graph_access &graph, ReductionConfig &config, NodeWeight weight_bound, bool apply_solution)
{
    if (graph.number_of_nodes() == 0)
    {
        solution = 0;
        return true;
    }
    if (graph.number_of_edges() == 0)
    {
        solution = 0;
        forall_nodes(graph, node)
        {
            if (graph.getNodeWeight(node) > 0)
            {
                graph.setPartitionIndex(node, 1);
                solution += graph.getNodeWeight(node);
            }
        }
        endfor return true;
    }
    auto c = config;
    c.disable_heuristic_exclude = true;
    c.disable_heuristic_include = true;
    c.use_partition_cover = false;
    c.disable_critical_set = true;
    c.disable_heavy_set = true;
    c.disable_blow_up = true;
    c.disable_generalized_fold = true;
    // c.time_limit = graph.number_of_nodes() / 10.0;
    c.time_limit = config.reduction_time_limit*0.1;
    c.max_swaps = 1000;

    branch_and_reduce_algorithm solver(graph, c, true);
    solver.ch.disable_cout();
    if (config.disable_early_termination)
        weight_bound = std::numeric_limits<NodeWeight>::max();

    int lb = greedy_lb(graph);

    if (lb > weight_bound)
    {
        solver.ch.enable_cout();
        return false;
    }

    int ub = greedy_ub(graph);

    if (weight_bound != -1 && ub <= weight_bound)
    {
        // printf("UB hit %d %d %d\n", weight_bound, lb, ub);
        solver.ch.enable_cout();
        return true;
    }

    // timer t;
    bool solved = solver.run_branch_reduce(weight_bound);
    // double t0 = t.elapsed();
    // printf("%lf %d %d %d %d %d %d\n", t0, graph.number_of_nodes(), weight_bound, lb, ub, solver.get_current_is_weight(), solved);

    if (!solved)
    {
        solver.ch.enable_cout();
        return false;
    }
    if (apply_solution)
    {
        solver.apply_branch_reduce_solution(graph);
    }
	solver.ch.enable_cout();
    solution = solver.get_is_weight();
    return true;
}
bool general_reduction::is_reduced(NodeID v, branch_and_reduce_algorithm *br_alg)
{
    return br_alg->status.node_status[v] != IS_status::not_set;
}
