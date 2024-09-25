
#include "general_reduction.h"

#include "reduce_algorithm.h"
#include "../tiny_solver/tiny_solver.h" 
#include "bounds.h"

typedef reduce_algorithm::IS_status IS_status;

bool general_reduction::is_suited(NodeID v, reduce_algorithm *br_alg)
{
    return br_alg->status.node_status[v] == IS_status::not_set;
}
NodeID general_reduction::get_max_weight_neighbor(NodeID v, reduce_algorithm *br_alg)
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
NodeWeight general_reduction::get_neighborhood_weight(NodeID v, reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    return std::accumulate(status.graph[v].begin(), status.graph[v].end(), 0, [&status](NodeWeight sum, NodeID neighbor)
                           { 
        if (status.node_status[neighbor] == IS_status::not_set) return sum + status.weights[neighbor];
        else return sum; });
}
void general_reduction::get_neighborhood_set(NodeID v, reduce_algorithm *br_alg, fast_set &set)
{
    auto &status = br_alg->status;
    set.clear();
    for (NodeID neighbor : status.graph[v])
    {
        if (status.node_status[neighbor] == IS_status::not_set)
            set.add(neighbor);
    }
}
void general_reduction::get_neighborhood_vector(NodeID v, reduce_algorithm *br_alg, std::vector<NodeID> &vec)
{
    auto &status = br_alg->status;
    vec.clear();
    for (NodeID neighbor : status.graph[v])
    {
        if (status.node_status[neighbor] == IS_status::not_set)
            vec.push_back(neighbor);
    }
}
bool general_reduction::try_neighborhood_reduction(NodeID v, reduce_algorithm *br_alg, NodeWeight neighbors_weight)
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
bool general_reduction::solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, std::vector<NodeID> &nodes_vec, const fast_set &nodes_set)
{
    auto &solver = br_alg->subgraph_solver;
    if (nodes_vec.size() == 0)
    {
        solution = 0;
        return true;
    }
    br_alg->tiny_solver_solve_subgraph(weight_bound, nodes_vec, nodes_set);
    solution = solver->independent_set_weight;
    return !(solver->time_limit_exceeded || solver->node_limit_exceeded || solver->weight_limit_exceeded);
}
bool general_reduction::solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, NodeID v)
{
    auto &solver = br_alg->subgraph_solver;
    br_alg->tiny_solver_solve_neighbourhood(weight_bound, v);
    solution = solver->independent_set_weight;
    return !(solver->time_limit_exceeded || solver->node_limit_exceeded || solver->weight_limit_exceeded);
}

bool general_reduction::is_reduced(NodeID v, reduce_algorithm *br_alg)
{
    return br_alg->status.node_status[v] != IS_status::not_set;
}
// used for generating training data
bool general_reduction::solve_induced_subgraph_from_set(NodeWeight weight_bound, NodeWeight &solution,  reduce_algorithm *br_alg, std::vector<NodeID> &nodes_vec, const fast_set &nodes_set, int &label)
{
    if (nodes_vec.size() == 0)
    {
        solution = 0;
        label = 1;
        return true;
    }
    if (nodes_vec.size() == 1)
    {
        solution = br_alg->status.weights[nodes_vec[0]];
        label = 1;
        return true;
    }
    auto &solver = br_alg->subgraph_solver;
    solve_induced_subgraph_from_set(weight_bound, solution, br_alg, nodes_vec, nodes_set);
    if (solver->time_limit_exceeded)
    {
        label = 2;
    } 
    else if (solver->weight_limit_exceeded || solver->node_limit_exceeded)
    {
        label = 0;
    }
    else
    {
        label = 1;
        return true; 
    }
    return false;
}
bool general_reduction::solve_induced_neighborhood_subgraph(NodeWeight weight_bound, NodeWeight &solution, reduce_algorithm *br_alg, NodeID v, int &label)
{
    auto &solver = br_alg->subgraph_solver;
    solve_induced_neighborhood_subgraph(weight_bound, solution, br_alg, v);
    if (solver->time_limit_exceeded)
    {
        label = 2;
    } 
    else if (solver->weight_limit_exceeded || solver->node_limit_exceeded)
    {
        label = 0;
    }
    else
    {
        label = 1;
        return true; 
    }
    return false;
}
