
#ifndef STRUCTION_REDUCTIONS_H
#define STRUCTION_REDUCTIONS_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "extended_struction.h"
#include "original_struction.h"
#include "key_functions.h"
#include "MaxHeap.h"
#include "reduction_config.h"
#include "reductions.h"

// system includes
#include <vector>
#include <memory>
#include <array>

class branch_and_reduce_algorithm;
using Struction_Type = ReductionConfig::Struction_Type;
using Key_Type = ReductionConfig::Key_Type;


template <typename struction_type, reduction_type type, int vertex_increase>
struct iterative_struction : public general_reduction
{
    iterative_struction(size_t n) : general_reduction(n) { has_filtered_marker = false; }
    ~iterative_struction() {}
    virtual iterative_struction *clone() const final { return new iterative_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual std::string get_reduction_name() final { 
        if (type == struction_decrease)
            return "struction_decrease";
        else
            return "struction_plateau";
        }
    virtual std::string get_model_path() final
    {
        if (type == struction_decrease)
            return "models/decreasing_struction.gnn";
        else
            return "";
    }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    virtual bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final { s.restore(br_alg); };
    virtual void apply(branch_and_reduce_algorithm *br_alg) final { s.apply(br_alg); };

private:
    struction_type s;
};

using reduction_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_decrease, -1>;
using plateu_struction = iterative_struction<extended_struction<false>, reduction_type ::struction_plateau, 0>;

template <typename key_function = ApproximateIncreaseKey, typename struction_type = extended_struction<false>, reduction_type type = reduction_type ::struction_blow /**/>
struct blow_up_struction : public general_reduction
{
    blow_up_struction(const ReductionConfig &config, size_t n) : general_reduction(n), update_set(n), f(config), generator(config.seed),
                                                                 distribution(0, 1) {}
    ~blow_up_struction() {}
    virtual blow_up_struction *clone() const final { return new blow_up_struction(*this); }

    virtual reduction_type get_reduction_type() const final { return type; }
    virtual std::string get_reduction_name() final { return "blow_up_struction"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    // virtual bool reduce_vertex(branch_and_reduce_algorithm* br_alg, NodeID v) final;
    virtual void restore(branch_and_reduce_algorithm *br_alg) final;
    virtual void apply(branch_and_reduce_algorithm *br_alg) final { s.apply(br_alg); };
    virtual void reset(branch_and_reduce_algorithm *br_alg, size_t comp_size) final;

private:
    size_t phase_start_kernel;

    fast_set update_set;
    std::vector<std::pair<NodeID, Gain>> update_list;
    std::vector<NodeID> added_list;

    key_function f;
    struction_type s;
    bool restored = false;
    size_t blow_ups;
    branch_and_reduce_algorithm *br_alg;

    MaxHeap<float> queue;

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution;

    void init_blow_up_phase();

    bool clean_up_queue();
    bool is_done();

    Gain denoise(float key)
    {
        return static_cast<Gain>(key);
    }

    float apply_noise(Gain key)
    {
        return key + distribution(generator);
    }

    void update_queue_by_key(NodeID n, Gain key)
    {
        if (queue.contains(n))
        {
            if (blow_ups != 0 && update_set.add(n))
                update_list.emplace_back(n, denoise(queue.getKey(n)));
            queue.changeKey(n, apply_noise(key));
        }
        else
        {
            if (blow_ups != 0 && update_set.add(n))
                added_list.push_back(n);
            queue.insert(n, apply_noise(key));
        }
    }

    void update_queue(NodeID n);
};


template <reduction_type type, int new_nodes>
reduction_ptr make_iterative_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_decreasing_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_plateau_struction(const ReductionConfig &config, size_t n);

reduction_ptr make_increasing_struction(const ReductionConfig &config, size_t n);

#endif // STRUCTION_REDUCTIONS_H
