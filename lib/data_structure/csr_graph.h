#pragma once
#include <stdio.h>
#include "definitions.h"

typedef struct
{
    NodeID N;
    EdgeID *V;
    NodeID *E;
    NodeWeight *W;
} csr_graph;

csr_graph graph_parse(FILE *f);

void graph_free(csr_graph g);

int graph_validate(NodeID N, const EdgeID *V, const NodeID *E);

csr_graph graph_subgraph(csr_graph g, NodeID *mask, NodeID *reverse_map);