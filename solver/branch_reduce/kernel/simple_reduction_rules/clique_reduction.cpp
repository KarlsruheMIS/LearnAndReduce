#include "clique_reduction.h"
#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;
bool clique_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
// if (br_alg->config.disable_clique) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    auto &status = br_alg->status;
    auto &set_1 = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    auto &isolated = br_alg->buffers[1];
    std::vector<NodeID> non_isolated;

    size_t oldn = status.remaining_nodes;

    for_each_changed_vertex(br_alg, [&](NodeID v)
                            {
                                neighbors.clear();
                                set_1.clear();
                                set_1.add(v);

                                // find potential clique
                                for (NodeID neighbor : status.graph[v])
                                {
                                    neighbors.push_back(neighbor);
                                    set_1.add(neighbor);
                                }

                                // check if clique
                                isolated.clear();
                                isolated.push_back(v);
                                non_isolated.clear();

                                size_t max_isolated_idx = 0;
                                weighted_node max_isolated{v, status.weights[v]};
                                weighted_node max_non_isolated{0, 0};

                                for (auto neighbor : neighbors)
                                {
                                    size_t count = 0;
                                    bool is_isolated = true;

                                    for (NodeID neighbor_2nd : status.graph[neighbor])
                                    {
                                        if (set_1.get(neighbor_2nd))
                                            count++;
                                        else
                                            is_isolated = false;
                                    }
                                    if (count != neighbors.size())
                                        return;

                                    if (is_isolated)
                                    {
                                        isolated.push_back(neighbor);
                                        if (status.weights[neighbor] > max_isolated.weight)
                                        {
                                            max_isolated = {neighbor, status.weights[neighbor]};
                                            max_isolated_idx = isolated.size() - 1;
                                        }
                                    }
                                    else
                                    {
                                        non_isolated.push_back(neighbor);
                                        if (status.weights[neighbor] > max_non_isolated.weight)
                                        {
                                            max_non_isolated = {neighbor, status.weights[neighbor]};
                                        }
                                    }
                                }

                                // one of "isolated" members has highest weight of clique: Add to IS
                                // also handles completely isolated cliques
                                if (max_isolated.weight >= max_non_isolated.weight)
                                {
                                    br_alg->set(max_isolated.node, IS_status::included);
                                    return;
                                }

                                // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
                                isolated[max_isolated_idx] = isolated.back();
                                isolated.pop_back();

                                for (auto neighbor : isolated)
                                {
                                    br_alg->set(neighbor, IS_status::excluded);
                                }

                                for (size_t i = 0; i < non_isolated.size(); i++)
                                {
                                    NodeID neighbor = non_isolated[i];
                                    if (status.weights[neighbor] <= max_isolated.weight)
                                    {
                                        br_alg->set(neighbor, IS_status::excluded);
                                        non_isolated[i] = non_isolated.back();
                                        non_isolated.pop_back();
                                        i--;
                                    }
                                }

                                fold(br_alg, std::move(max_isolated), std::move(non_isolated)); });

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
	// if (oldn != status.remaining_nodes) std::cout << "clique redu -> " << (oldn - status.remaining_nodes) << std::endl;
    return oldn != status.remaining_nodes;
}

inline bool clique_reduction::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    auto &status = br_alg->status;
    auto &set_1 = br_alg->set_1;
    auto &neighbors = br_alg->buffers[0];
    auto &isolated = br_alg->buffers[1];
    std::vector<NodeID> non_isolated;
    size_t oldn = status.remaining_nodes;

    get_neighborhood_set(v, br_alg, set_1);
    get_neighborhood_vector(v, br_alg, neighbors);
    set_1.add(v);

    // check if clique
    // non_isolated.clear();
    isolated.clear();
    isolated.push_back(v);

    size_t max_isolated_idx = 0;
    weighted_node max_isolated{v, status.weights[v]};
    weighted_node max_non_isolated{0, 0};

    for (auto neighbor : neighbors)
    {
        size_t count = 0;
        bool is_isolated = true;

        for (NodeID neighbor_2nd : status.graph[neighbor])
        {
            if (set_1.get(neighbor_2nd))
                count++;
            else
                is_isolated = false;
        }

        if (is_isolated)
        {
            isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_isolated.weight)
            {
                max_isolated = {neighbor, status.weights[neighbor]};
                max_isolated_idx = isolated.size() - 1;
            }
        }
        else
        {
            non_isolated.push_back(neighbor);
            if (status.weights[neighbor] > max_non_isolated.weight)
            {
                max_non_isolated = {neighbor, status.weights[neighbor]};
            }
        }

        if (count != neighbors.size())
            return false;
    }

    // one of "isolated" members has highest weight of clique: Add to IS
    // also handles completely isolated cliques
    if (max_isolated.weight >= max_non_isolated.weight)
    {
        br_alg->set(max_isolated.node, IS_status::included);
        return oldn != br_alg->status.remaining_nodes;
    }

    // remove all nodes from the clique which have a smaller or eqaul weight than "max_isolated" -> we can always pick "max_isolated" over them
    isolated[max_isolated_idx] = isolated.back();
    isolated.pop_back();

    for (auto neighbor : isolated)
    {
        br_alg->set(neighbor, IS_status::excluded);
    }

    for (size_t i = 0; i < non_isolated.size(); i++)
    {
        NodeID neighbor = non_isolated[i];
        if (status.weights[neighbor] <= max_isolated.weight)
        {
            br_alg->set(neighbor, IS_status::excluded);
            non_isolated[i] = non_isolated.back();
            non_isolated.pop_back();
            i--;
        }
    }

    fold(br_alg, std::move(max_isolated), std::move(non_isolated));

    return oldn != status.remaining_nodes;
}
void clique_reduction::fold(branch_and_reduce_algorithm *br_alg, const weighted_node &isolated, std::vector<NodeID> &&non_isolated)
{
    auto &status = br_alg->status;

    br_alg->set(isolated.node, IS_status::folded);
    status.reduction_offset += isolated.weight;

    for (auto node : non_isolated)
    {
        status.weights[node] -= isolated.weight;
        br_alg->add_next_level_neighborhood(node);
    }

    status.folded_stack.push_back(get_reduction_type());
    br_alg->add_next_level_neighborhood(non_isolated);

    restore_vec.emplace_back(isolated, std::move(non_isolated));
}
void clique_reduction::restore(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto &data = restore_vec.back();

    br_alg->unset(data.isolated.node);
    status.reduction_offset -= data.isolated.weight;

    for (auto node : data.non_isolated)
    {
        status.weights[node] += data.isolated.weight;
    }

    restore_vec.pop_back();
}
void clique_reduction::apply(branch_and_reduce_algorithm *br_alg)
{
    auto &status = br_alg->status;
    auto isolated = restore_vec.back().isolated.node;

    bool set_isolated = true;

    for (auto node : restore_vec.back().non_isolated)
    {
        if (status.node_status[node] == IS_status::included)
        {
            set_isolated = false;
            break;
        }
    }

    status.is_weight += restore_vec.back().isolated.weight;

    restore(br_alg);

    if (set_isolated)
    {
        status.node_status[isolated] = IS_status::included;
    }
    else
    {
        status.node_status[isolated] = IS_status::excluded;
    }
}
