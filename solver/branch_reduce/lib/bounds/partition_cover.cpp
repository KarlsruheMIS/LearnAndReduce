#include "partition_cover.h"

#include "random_functions.h"
#include "graph_partitioner.h"
#include "reduction_config.h"
#include "partition_config.h"
#include "graph_extractor.h"
#include "definitions.h"
#include "branch_and_reduce_algorithm.h"
#include "../app/balance_configuration.h"
#include "../app/configuration.h"
#include "solution_check.h"

partition_cover::partition_cover(PartitionID number_of_blocks) : k(number_of_blocks)
{
}

partition_cover::~partition_cover() {}

void partition_cover::create_partition(graph_access &G, ReductionConfig &config)
{
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = k;
        cfg.fast(partition_config);
        // cfg.strong(partition_config);
        // cfg.strongsocial(partition_config);
        // cfg.fastsocial(partition_config);
        // cfg.eco(partition_config);

        G.set_partition_count(partition_config.k); 
 
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
       
        // if(vwgt != NULL) {
                forall_nodes(G, node) {
                        // G.setNodeWeight(node, 1);
                        G.setNodeWeight(node, G.getNodeDegree(node)+1);
                } endfor
        // }

        // if(adjcwgt != NULL) {
        //         forall_edges(G, e) {
        //                 G.setEdgeWeight(e, adjcwgt[e]);
        //         } endfor 
        // }

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);
        partition_config.largest_graph_weight = G.number_of_nodes();
        
        double epsilon = partition_config.imbalance/100;
        partition_config.upper_bound_partition = ceil((1+epsilon)*partition_config.largest_graph_weight/(double)partition_config.k);
        partition_config.graph_allready_partitioned = false;

        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);
}


NodeWeight partition_cover::solve_partition(graph_access &G, ReductionConfig &config)
{
        NodeWeight global_partition_solution_weight = 0;
        graph_extractor extractor;
        std::vector<NodeID> component_mapping;
        graph_access component_graph;
        ReductionConfig component_config = config;
        component_config.use_partition_cover = false;

        for (size_t i = 0; i < k; i++)
        {
                component_mapping.clear();
                extractor.extract_block(G, component_graph, i, component_mapping);
                // #ifdef DEBUG
                // solution_check<graph_access> sc(component_graph);
                // assert(sc.check_graph());
                // #endif

                branch_and_reduce_algorithm solver(component_graph, component_config, true);
                solver.ch.disable_cout();
                if (solver.run_branch_reduce()) {
                        global_partition_solution_weight += solver.get_current_is_weight();
                        solver.ch.enable_cout();
                } else {
                        solver.ch.enable_cout();
                        std::cout << "component not solved optimally" << std::endl;
                        global_partition_solution_weight = std::numeric_limits<NodeWeight>::max();
                        break;
                }
        }
        return global_partition_solution_weight;
}

