#include "heuristic_exclude_reduction.h"

#include "general_reduction.h"
#include "reductions.h"
#include "branch_and_reduce_algorithm.h"

typedef branch_and_reduce_algorithm::IS_status IS_status;

bool heuristic_exclude_reduction::reduce(branch_and_reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_heuristic_exclude) return false;
    if (br_alg->blowing_up)
        return false;
    size_t oldn = br_alg->status.remaining_nodes;
    auto &weights = br_alg->status.weights;
    auto &config = br_alg->config;
#ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
#endif
    if (marker.current.empty())
    {
        has_run = false; // check everything in next round again
        has_filtered_marker = true;
        return false;
    }

    auto &unsafe_to_reduce = br_alg->set_1;
    unsafe_to_reduce.clear();
    br_alg->heuristically_reducing = true;
    if (config.heuristic_style == ReductionConfig::Heuristic_Style::multiple_safe ||
        config.heuristic_style == ReductionConfig::Heuristic_Style::multiple_very_safe)
    {
        // std::sort(marker.current.begin(), marker.current.end(), [&](NodeID a, NodeID b) {
        //     return weights[a] - get_neighborhood_weight(a, br_alg) < weights[b] - get_neighborhood_weight(b, br_alg);
        // });

        for (NodeID v : marker.current)
        {
            if (unsafe_to_reduce.get(v))
                continue;
            if (is_reduced(v, br_alg))
                continue;

            if (br_alg->status.weights[v] > get_neighborhood_weight(v, br_alg))
            {
                br_alg->set(v, IS_status::included);
            }
            else
            {
                br_alg->set(v, IS_status::excluded);
            }
            unsafe_to_reduce.add(v);
            for (NodeID u : br_alg->status.graph[v])
            {
                unsafe_to_reduce.add(u);
                if (config.heuristic_style == ReductionConfig::Heuristic_Style::multiple_very_safe)
                {
                    for (NodeID w : br_alg->status.graph[u])
                    {
                        unsafe_to_reduce.add(w);
                    }
                }
            }
        }
    }
    else if (config.heuristic_style == ReductionConfig::Heuristic_Style::all)
    {
        // assert(config.gnn_filter);
        for (NodeID v : marker.current)
        {
            if (is_reduced(v, br_alg))
                continue;
            if (br_alg->status.weights[v] > get_neighborhood_weight(v, br_alg))
            {
                br_alg->set(v, IS_status::included);
            }
            else
            {
                br_alg->set(v, IS_status::excluded);
            }
        }
    }
    else
    { // only exclude single vertex with smallest score
        // printf("%d\n", marker.current[0]);
        NodeID v = marker.current[0];

        // NodeID v = std::min_element(marker.current.begin(), marker.current.end(), [&](NodeID a, NodeID b) {
        //     return weights[a] - get_neighborhood_weight(a, br_alg) < weights[b] - get_neighborhood_weight(b, br_alg);
        // })[0];
        assert(!is_reduced(v, br_alg));
        br_alg->set(v, IS_status::excluded);
    }

#ifdef REDUCTION_INFO
    reduced_nodes += (oldn - br_alg->status.remaining_nodes);
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    br_alg->status.heuristically_reduced_n += (oldn - br_alg->status.remaining_nodes);
    has_run = false; // check everything in next round again
    has_filtered_marker = true;
    return oldn != br_alg->status.remaining_nodes;
}
