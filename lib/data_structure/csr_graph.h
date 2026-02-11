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
