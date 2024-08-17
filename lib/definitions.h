/******************************************************************************
 * definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DEFINITIONS_H_CHR
#define DEFINITIONS_H_CHR

#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include "limits.h"
#include "stdio.h"


/**********************************************
 * Constants
 * ********************************************/
//Types needed for the graph ds
typedef unsigned int 	NodeID;
typedef unsigned int 	PartitionID;
typedef int 		EdgeWeight;

#ifdef MODE64BITNODEWEIGHTS
typedef unsigned long long 	NodeWeight;
#else
typedef unsigned int 	NodeWeight;
#endif
#ifdef MODE64BITEDGES
typedef uint64_t 	EdgeID;
#else
typedef unsigned int 	EdgeID;
#endif
typedef long FlowType;
typedef NodeWeight Gain;
typedef std::vector<NodeID> NodePermutationMap;


#endif

