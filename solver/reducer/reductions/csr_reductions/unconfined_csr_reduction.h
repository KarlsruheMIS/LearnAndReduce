
#ifndef UNCONFINED_CSR_REDUCTION_H
#define UNCONFINED_CSR_REDUCTION_H
// local includes
#include "definitions.h"
#include "fast_set.h"
#include "general_reduction.h"
#include "reduction_config.h"

#ifdef __cplusplus
extern "C"
{
#endif
#include "csr_graph.h"
#ifdef __cplusplus
}
#endif

class reduce_algorithm;

struct unconfined_csr_reduction : public general_reduction
{
    unconfined_csr_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = true;
    }
    ~unconfined_csr_reduction() {}
    virtual unconfined_csr_reduction *clone() const final { return new unconfined_csr_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::unconfined_csr; }
    virtual std::string get_reduction_name() final { return "unconfined"; }
    virtual std::string get_model_path() final { return "models/unconfined.lr_gcn"; }
    virtual bool reduce(reduce_algorithm *br_alg) final;
    bool reduce_vertex(reduce_algorithm *br_alg, NodeID v, csr_graph &g, int *A, void *R);
    void generate_global_data(reduce_algorithm *br_alg, std::vector<std::vector<int>> &reduction_data, int reduction_index);
};
#endif // UNCONFINED_CSR_REDUCTION_H
