#include "heuristic_include_reduction.h"

#include "general_reduction.h"
#include "reductions.h"
#include "reduce_algorithm.h"

typedef reduce_algorithm::IS_status IS_status;

bool heuristic_include_reduction::reduce(reduce_algorithm *br_alg)
{
    // if (br_alg->config.disable_heuristic_include) return false;
    if (br_alg->blowing_up)
        return false;
    size_t oldn = br_alg->status.remaining_nodes;
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
    auto &weights = br_alg->status.weights;
    br_alg->heuristically_reducing = true;

    auto &unsafe_to_reduce = br_alg->set_1;
    unsafe_to_reduce.clear();
    if (config.heuristic_style != ReductionConfig::Heuristic_Style::single) {
        std::sort(marker.current.begin(), marker.current.end(), [&](NodeID a, NodeID b) {
            return weights[a] - get_neighborhood_weight(a, br_alg) > weights[b] - get_neighborhood_weight(b, br_alg);
        });

        for (NodeID v : marker.current)
        {
            if (unsafe_to_reduce.get(v))
                continue;
            if (is_reduced(v, br_alg))
                continue;

            br_alg->set(v, IS_status::included);
            unsafe_to_reduce.add(v);

            if (config.heuristic_style != ReductionConfig::Heuristic_Style::all)
            {
                for (NodeID u : br_alg->status.graph[v])
                {
                    unsafe_to_reduce.add(u);
                    for (NodeID w : br_alg->status.graph[u])
                    {
                        unsafe_to_reduce.add(w);
                        if (config.heuristic_style == ReductionConfig::Heuristic_Style::multiple_very_safe)
                        {
                            for (NodeID x : br_alg->status.graph[w])
                            {
                                unsafe_to_reduce.add(x);
                            }
                        }
                    }
                }
            }
        }
    }
    else if (config.heuristic_style == ReductionConfig::Heuristic_Style::single)
    {

        // only include the vertex with the highest score  without sorting whole marker
        NodeID v = std::max_element(marker.current.begin(), marker.current.end(), [&](NodeID a, NodeID b)
                                    {
            // return weights[a] - get_neighborhood_weight(a, br_alg) < weights[b] - get_neighborhood_weight(b, br_alg);
            return weights[a] < weights[b]; })[0];
        assert(!is_reduced(v, br_alg));
        br_alg->set(v, IS_status::included);
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
