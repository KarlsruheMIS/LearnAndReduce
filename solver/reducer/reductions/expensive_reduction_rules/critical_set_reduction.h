/******************************************************************************
 * critical_set_reduction.h
 *
 * Copyright (C) 2015-2018 Robert Williger
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *****************************************************************************/

#ifndef CRITICAL_SET_H
#define CRITICAL_SET_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "reductions.h"

struct critical_set_reduction : public general_reduction
{
    critical_set_reduction(size_t n) : general_reduction(n) { has_filtered_marker = true;}
    ~critical_set_reduction() {}
    virtual critical_set_reduction *clone() const final { return new critical_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::critical_set; }
    virtual std::string get_reduction_name() final { return "critical_set"; }
    virtual std::string get_model_path() final { return get_models_folder_path() + "critical_set.lr_gcn"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);
};

#endif // CRITICAL_SET_H
