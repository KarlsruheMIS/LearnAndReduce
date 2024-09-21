/******************************************************************************
 * reductions.cpp
 *
 * Copyright (C) 2015-2018 Robert Williger
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "critical_set_reduction.h"

#include "reduce_algorithm.h"
#include "flow_graph.h"
#include "push_relabel.h"

typedef reduce_algorithm::IS_status IS_status;

bool critical_set_reduction::reduce(reduce_algorithm *br_alg)
{
    if (br_alg->blowing_up)
        return false;
    if (br_alg->heuristically_reducing)
        return false;
// if (br_alg->config.disable_critical_set) return false;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif

    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    auto &mapping = br_alg->buffers[0];
    auto &inverse_mapping = br_alg->buffers[1];
    mapping.resize(status.n);
    inverse_mapping.resize(status.remaining_nodes);
    NodeID n = 0;
    for (NodeID node = 0; node < br_alg->status.n; node++)
    {
        if (status.node_status[node] == IS_status::not_set)
        {
            mapping[node] = n;
            inverse_mapping[n] = node;
            ++n;
        }
    }

    // build bipartite flow graph
    // node '+ n' shows that we refer to the node in the rest[1] partition
    flow_graph fg;
    fg.start_construction(2 * n + 2);

    const NodeID source = 2 * n;
    const NodeID sink = source + 1;

    for (NodeID id = 0; id < n; id++)
    {
        NodeID node = inverse_mapping[id];
        // add source and target edges
        fg.new_edge(source, id, status.weights[node]);
        fg.new_edge(id + n, sink, status.weights[node]);

        // add edges between node and its neighbors in the other partition
        for (NodeID neighbor : status.graph[node])
        {
            // each outgoing edge has enough capacity to support full flow of the single incoming edge into 'node'
            fg.new_edge(id, mapping[neighbor] + n, status.weights[node]);
        }
    }
    fg.finish_construction();

    // solve max-flow problem
    push_relabel flow_solver;
    std::vector<NodeID> dummy_vec;
    flow_solver.solve_max_flow_min_cut(fg, source, sink, false, dummy_vec);

    auto &max_cs_set = br_alg->double_set;

    max_cs_set.clear();
    // (source, node) edges where flow < capacity indicate that node is in the maximum critical set
    forall_out_edges(fg, edge, source)
    {
        NodeID id = fg.getEdgeTarget(source, edge);
        if (fg.getEdgeFlow(source, edge) < fg.getEdgeCapacity(source, edge))
        {
            max_cs_set.add(inverse_mapping[id]);
        }
    }
    endfor

        // isolated nodes in the maximum critical set form the maximum independent critical set
        for (NodeID id = 0; id < n; id++)
    {
        NodeID node = inverse_mapping[id];
        if (max_cs_set.get(node))
        {
            bool isolated = true;

            for (NodeID neighbor : status.graph[node])
            {
                if (max_cs_set.get(neighbor))
                {
                    isolated = false;
                    break;
                }
            }

            if (isolated)
            {
                // found isolated node
                br_alg->set(node, IS_status::included);
            }
        }
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    // if (oldn != status.remaining_nodes) std::cout << "critical redu -> " << (oldn - status.remaining_nodes) << std::endl;
    return oldn != status.remaining_nodes;
}
void critical_set_reduction::generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index)
{
    auto &status = br_alg->status;
    size_t oldn = status.remaining_nodes;

    auto &mapping = br_alg->buffers[0];
    auto &inverse_mapping = br_alg->buffers[1];
    mapping.resize(status.n);
    inverse_mapping.resize(status.remaining_nodes);
    NodeID n = 0;
    for (NodeID node = 0; node < br_alg->status.n; node++)
    {
        if (status.node_status[node] == IS_status::not_set)
        {
            mapping[node] = n;
            inverse_mapping[n] = node;
            ++n;
        }
    }

    // build bipartite flow graph
    // node '+ n' shows that we refer to the node in the rest[1] partition
    flow_graph fg;
    fg.start_construction(2 * n + 2);

    const NodeID source = 2 * n;
    const NodeID sink = source + 1;

    for (NodeID id = 0; id < n; id++)
    {
        NodeID node = inverse_mapping[id];
        // add source and target edges
        fg.new_edge(source, id, status.weights[node]);
        fg.new_edge(id + n, sink, status.weights[node]);

        // add edges between node and its neighbors in the other partition
        for (NodeID neighbor : status.graph[node])
        {
            // each outgoing edge has enough capacity to support full flow of the single incoming edge into 'node'
            fg.new_edge(id, mapping[neighbor] + n, status.weights[node]);
        }
    }
    fg.finish_construction();

    // solve max-flow problem
    push_relabel flow_solver;
    std::vector<NodeID> dummy_vec;
    flow_solver.solve_max_flow_min_cut(fg, source, sink, false, dummy_vec);

    auto &max_cs_set = br_alg->double_set;

    max_cs_set.clear();
    // (source, node) edges where flow < capacity indicate that node is in the maximum critical set
    forall_out_edges(fg, edge, source)
    {
        NodeID id = fg.getEdgeTarget(source, edge);
        if (fg.getEdgeFlow(source, edge) < fg.getEdgeCapacity(source, edge))
        {
            max_cs_set.add(inverse_mapping[id]);
        }
    }
    endfor

        // isolated nodes in the maximum critical set form the maximum independent critical set
        for (NodeID id = 0; id < n; id++)
    {
        NodeID node = inverse_mapping[id];
        if (max_cs_set.get(node))
        {
            bool isolated = true;

            for (NodeID neighbor : status.graph[node])
            {
                if (max_cs_set.get(neighbor))
                {
                    isolated = false;
                    break;
                }
            }

            if (isolated)
            {
                reduction_data[reduction_index][node] = 1;
            }
        }
    }
}
