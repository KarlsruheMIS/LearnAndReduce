
#ifndef UNCONFINED_CSR_REDUCTION_H
#define UNCONFINED_CSR_REDUCTION_H
// local includes
#include "definitions.h"
#include "fast_set.h"
#include "general_reduction.h"
#include "reduction_config.h"

#ifdef __cplusplus
extern "C" {
#endif
    #include "include/graph.h"
#ifdef __cplusplus
}
#endif

class branch_and_reduce_algorithm;

struct unconfined_csr_reduction : public general_reduction
{
    unconfined_csr_reduction(size_t n) : general_reduction(n)
    {
        has_filtered_marker = false;
    }
    ~unconfined_csr_reduction() { }
    virtual unconfined_csr_reduction *clone() const final { return new unconfined_csr_reduction(*this); }

    virtual reduction_type get_reduction_type() const final { return reduction_type::unconfined_csr; }
    virtual std::string get_reduction_name() final { return "unconfined_csr"; }
    virtual bool reduce(branch_and_reduce_algorithm *br_alg) final;
    bool reduce_vertex(branch_and_reduce_algorithm *br_alg, NodeID v, graph &g, int *A, void *R);
    bool generate_global_data(branch_and_reduce_algorithm *br_alg, std::vector<NodeID> &label);

    graph build_graph(branch_and_reduce_algorithm *br_alg);

};
#endif // UNCONFINED_CSR_REDUCTION_H
