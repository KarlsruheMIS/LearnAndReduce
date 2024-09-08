#pragma once
#include "definitions.h"

void *reduction_init(NodeID N, EdgeID M);

void reduction_free(void *R);

int reduction_neighborhood_csr(void *R, NodeID N, const EdgeID *V, const NodeID *E,
                               const NodeWeight *W, const int *A, NodeID u);

int reduction_unconfined_csr(void *R, NodeID N, const EdgeID *V, const NodeID *E,
                             const NodeWeight *W, const int *A, NodeID u);