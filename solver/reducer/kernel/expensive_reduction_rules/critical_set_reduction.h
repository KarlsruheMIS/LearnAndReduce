/******************************************************************************
 * reductions.h
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

#ifndef CRITICAL_SET_H
#define CRITICAL_SET_H

// local includes
#include "definitions.h"
#include "general_reduction.h"
#include "reductions.h"

struct critical_set_reduction : public general_reduction
{
    critical_set_reduction(size_t n) : general_reduction(n) {}
    ~critical_set_reduction() {}
    virtual critical_set_reduction *clone() const final { return new critical_set_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::critical_set; }
    virtual std::string get_reduction_name() final { return "critical_set"; }
    virtual std::string get_model_path() final { return "models/critical_set.lr_gcn"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);
};

#endif // CRITICAL_SET_H
