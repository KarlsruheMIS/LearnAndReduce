
#ifndef _REDUCTIONS_H
#define _REDUCTIONS_H

#define REDUCTION_INFO

// local includes
#include "definitions.h"
#include "general_reduction.h"

constexpr size_t REDUCTION_NUM = 29;

struct reduction_ptr
{
    general_reduction *reduction = nullptr;

    reduction_ptr() = default;

    ~reduction_ptr()
    {
        release();
    }

    reduction_ptr(general_reduction *reduction) : reduction(reduction){};

    reduction_ptr(const reduction_ptr &other) : reduction(other.reduction->clone()){};

    reduction_ptr &operator=(const reduction_ptr &other)
    {
        release();
        reduction = other.reduction->clone();
        return *this;
    };

    reduction_ptr(reduction_ptr &&other) : reduction(std::move(other.reduction))
    {
        other.reduction = nullptr;
    };

    reduction_ptr &operator=(reduction_ptr &&other)
    {
        reduction = std::move(other.reduction);
        other.reduction = nullptr;
        return *this;
    };

    general_reduction *operator->() const
    {
        return reduction;
    }

    void release()
    {
        if (reduction)
        {
            delete reduction;
            reduction = nullptr;
        }
    };
};

template <class Last>
void make_reduction_vector_helper(std::vector<reduction_ptr> &vec, size_t n)
{
    vec.emplace_back(new Last(n));
};

template <class First, class Second, class... Redus>
void make_reduction_vector_helper(std::vector<reduction_ptr> &vec, size_t n)
{
    vec.emplace_back(new First(n));
    make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template <class... Redus>
std::vector<reduction_ptr> make_reduction_vector(size_t n)
{
    std::vector<reduction_ptr> vec;
    make_reduction_vector_helper<Redus...>(vec, n);
    return vec;
};

#endif // REDUCTIONS_H
