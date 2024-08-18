#include "struction_reductions.h"

#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "key_functions.h"


typedef branch_and_reduce_algorithm::IS_status IS_status;

template<typename struction_type, reduction_type type, int vertex_increase>
bool iterative_struction<struction_type, type, vertex_increase>::reduce(branch_and_reduce_algorithm* br_alg) {
    #ifdef REDUCTION_INFO
        br_alg->reduction_timer.restart();
        NodeID oldn = br_alg->status.remaining_nodes;
    #endif
    bool applied = false;
    for_each_changed_vertex(br_alg, [&](NodeID v)
                            {
        if (reduce_vertex(br_alg, v))
            applied = true; });

#ifdef REDUCTION_INFO
    reduced_nodes += oldn - br_alg->status.remaining_nodes;
    reduction_time += br_alg->reduction_timer.elapsed();
#endif
    return applied;
}
template <typename struction_type, reduction_type type, int vertex_increase>
inline bool iterative_struction<struction_type, type, vertex_increase>::reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v)
{
    bool applied = false;
    if (br_alg->deg(v) > br_alg->config.struction_degree || !s.reduce(br_alg, v, s.removed_vertices(br_alg, v) + vertex_increase))
        return false;

    br_alg->status.folded_stack.push_back(get_reduction_type());
    applied = true;

    return applied;
}

template class iterative_struction<extended_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<extended_struction<true>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<false>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<false>, reduction_type::struction_plateau, 0>;
template class iterative_struction<original_struction<true>, reduction_type::struction_decrease, -1>;
template class iterative_struction<original_struction<true>, reduction_type::struction_plateau, 0>;

template <typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::reduce(branch_and_reduce_algorithm *br_alg)
{
    this->br_alg = br_alg;
    auto &status = br_alg->status;
    #ifdef REDUCTION_INFO
    br_alg->reduction_timer.restart();
    #endif
    init_blow_up_phase();

    while (!is_done() && clean_up_queue())
    {
        NodeID n = queue.maxElement();
        Gain key = denoise(queue.maxValue());
        size_t set_limit_by_key = f.set_limit(br_alg, n, key);
        size_t plain_set_limit = br_alg->config.set_limit;
        if (f.set_estimate(br_alg, n, key) > plain_set_limit)
            break;

        queue.deleteMax();
        if (s.reduce(br_alg, n, std::min(set_limit_by_key, plain_set_limit)))
        {
            // struction was successfully executed
            if (br_alg->config.backtrack_style != ReductionConfig::Backtrack_Type::IMMEDIATE_EXCLUDE && update_set.add(n))
                update_list.emplace_back(n, key);
            status.folded_stack.push_back(get_reduction_type());
            ++blow_ups;
            //"blow up" actually can reduce remaining nodes if multiple blow ups per phase are enabled.
            br_alg->min_kernel = std::min(status.remaining_nodes, br_alg->min_kernel);
        }
        else if (plain_set_limit <= set_limit_by_key)
        {
            break;
        }
        else
        {
            // reinsert with tighter bound
            update_queue_by_key(n, f.key_by_set_estimate(br_alg, n, set_limit_by_key + 1));
        }
    }
    #ifdef REDUCTION_INFO
    reduction_time += br_alg->reduction_timer.elapsed();
    #endif
    return blow_ups != 0;
}

template <typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::update_queue(NodeID n)
{
    if (br_alg->deg(n) > br_alg->config.struction_degree)
        return;
    update_queue_by_key(n, f.key(br_alg, n));
}

template <typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::clean_up_queue()
{
    auto &status = br_alg->status;
    update_set.resize(status.n);

    // Update candidate set
    for (NodeID n : marker.next)
        update_queue(n);
    marker.clear_next();

    while (queue.size())
    {
        NodeID n = queue.maxElement();
        if (status.node_status[n] == IS_status::not_set && br_alg->deg(n) <= br_alg->config.struction_degree)
            return true;
        if (update_set.add(n))
            update_list.emplace_back(n, denoise(queue.maxValue()));
        queue.deleteMax();
    }
    return false;
}

template <typename key_function, typename struction_type, reduction_type type>
bool blow_up_struction<key_function, struction_type, type>::is_done()
{
    auto &config = br_alg->config;
    auto &status = br_alg->status;
    size_t cur_kernel = status.remaining_nodes;
    size_t max_kernel = config.phase_blow_up_factor * phase_start_kernel;
    return max_kernel < cur_kernel || blow_ups == config.phase_blow_ups;
}

template <typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::restore(branch_and_reduce_algorithm *br_alg)
{
    s.restore(br_alg);
    restored = true;
}

template <typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::reset(branch_and_reduce_algorithm *br_alg, size_t comp_size)
{
    restored = false;
    queue.clear();
}

template <typename key_function, typename struction_type, reduction_type type>
void blow_up_struction<key_function, struction_type, type>::init_blow_up_phase()
{
    auto &status = br_alg->status;

    update_set.resize(status.n);
    blow_ups = 0;
    if (restored)
    {
        for (NodeID n : added_list)
            if (queue.contains(n))
                queue.deleteNode(n);
        for (auto &e : update_list)
            update_queue_by_key(e.first, e.second);
    }
    else
    {
        for_each_changed_vertex(br_alg, [&](NodeID n)
                                { update_queue(n); });
    }
    phase_start_kernel = status.remaining_nodes;
    update_list.clear();
    added_list.clear();
    update_set.clear();

    restored = false;
}

template class blow_up_struction<DegreeKey>;
template class blow_up_struction<IncreaseKey>;
template class blow_up_struction<ApproximateIncreaseKey>;
template class blow_up_struction<RandomKey>;

template <reduction_type type, int new_nodes>
reduction_ptr make_iterative_struction(const ReductionConfig &config, size_t n)
{
    const auto s = config.struction_type;
    if (s == Struction_Type::ORIGINAL)
        return reduction_ptr(new iterative_struction<original_struction<false>, type, new_nodes>(n));
    if (s == Struction_Type::MODIFIED)
        return reduction_ptr(new iterative_struction<original_struction<true>, type, new_nodes>(n));
    if (s == Struction_Type::EXTENDED_REDUCED)
        return reduction_ptr(new iterative_struction<extended_struction<true>, type, new_nodes>(n));
    return reduction_ptr(new iterative_struction<extended_struction<false>, type, new_nodes>(n));
};

reduction_ptr make_decreasing_struction(const ReductionConfig &config, size_t n)
{
    return make_iterative_struction<reduction_type::struction_decrease, -1>(config, n);
};

reduction_ptr make_plateau_struction(const ReductionConfig &config, size_t n)
{
    return make_iterative_struction<reduction_type::struction_plateau, 0>(config, n);
};

reduction_ptr make_increasing_struction(const ReductionConfig &config, size_t n)
{
    if (config.key_type == Key_Type::RANDOM)
        return reduction_ptr(new blow_up_struction<RandomKey>(config, n));
    if (config.key_type == Key_Type::DEGREE)
        return reduction_ptr(new blow_up_struction<DegreeKey>(config, n));
    if (config.key_type == Key_Type::INCREASE)
        return reduction_ptr(new blow_up_struction<IncreaseKey>(config, n));
    return reduction_ptr(new blow_up_struction<ApproximateIncreaseKey>(config, n));
};
