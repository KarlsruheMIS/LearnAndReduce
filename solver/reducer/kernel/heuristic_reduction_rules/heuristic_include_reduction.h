
#ifndef heuristic_INCLUDE_REDUCTION_H
#define heuristic_INCLUDE_REDUCTION_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "fast_set.h"
#include "reduction_config.h"


class reduce_algorithm;

struct heuristic_include_reduction : public general_reduction
{
    heuristic_include_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true; }
    ~heuristic_include_reduction() {has_filtered_marker = true;}
    virtual heuristic_include_reduction *clone() const final { return new heuristic_include_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::heuristic_include; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    virtual std::string get_reduction_name() final { return "heuristic_include"; }
    virtual std::string get_model_path() final { return ""; }

};

#endif // HEURISTIC_INCLUDE_REDUCTION_H 
