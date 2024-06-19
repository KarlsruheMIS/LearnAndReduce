#pragma once

#include <memory>
#include <vector>

#include "definitions.h"
#include "reduction_config.h"
#include "partition_config.h"
#include "graph_access.h"

class partition_cover {
public:
    partition_cover(PartitionID number_of_blocks);
    ~partition_cover();

    void create_partition(graph_access &G, ReductionConfig &config);


    NodeWeight solve_partition(graph_access &G, ReductionConfig &config);

private:
    // Graph data structures used for the KaHIP-library calls.
    int k;
};