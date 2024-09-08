#pragma once
#include <stdio.h>
#include "definitions.h"

typedef struct
{
    NodeID N;
    EdgeID *V;
    NodeID *E;
    NodeWeight *W;
} graph;

graph graph_parse(FILE *f);

void graph_free(graph g);

int graph_validate(NodeID N, const EdgeID *V, const NodeID *E);

graph graph_subgraph(graph g, int *mask, int *rm);