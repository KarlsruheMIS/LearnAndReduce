#!/bin/bash

parameters=" --weight_source='uniform' "
program="kernelization"
instance="parameter_tuning"
run_parallel="run_parallel_single_reductions.txt"
address="results/single_reductions/${instance}"

path_to_program="../deploy/${program}"
path_to_instance="/home/ernestineg/test_instances/"
database="/home/ernestineg/test_instances/${instance}"
mkdir -p $address

create_run() {
    id=$1
    id_instructions=$2
    address2=${address}/${id}
    identify=${instance}/${program}/${id}
    mkdir -p ${address2}
    
    instruction="${path_to_program} ${id_instructions}"
    echo " ${instruction} > ${address2}/result.txt ${solve_kernel} && sed -i '\|${identify}|d' ./${run_parallel} || echo '${address2}' " >> ${run_parallel}
}

# no struction since they will always be disabled
reductions=("neighborhood" "fold1" "fold2" "clique" "twin" "extended_twin" "domination" "extended_domination" "extended_domination_reverse"
            "clique_neighborhood_fast" "funnel" "funnel_fold" "unconfined" "basic_se" "extended_se"
            "critical_set" "cut_vertex"
            "generalized_fold" "heavy_set" "heavy_set3")
n_reductions=${#reductions[@]}
rm -f ${run_parallel}


    t=36000
            #########################################################
            ############### Experiments for soa  ####################
            #########################################################

    for seed in 0 1; do # seeds
        for graph in `ls ${database} | grep "graph" | awk -F"." '{print $1}'`; do
            sg=128
            config="all_reductions_cyclicFast"
            red_style="early_struction"
            for ((k=0; k<${n_reductions}; k++ )); do
                disable=" --disable_blow_up --disable_plateau_struction --disable_decreasing_struction "
                for ((j=0; j<${n_reductions}; j++ )); do
                    if [ $k -ne $j ]; then
                        disable+=" --disable_${reductions[j]} "
                    fi
                done
                overall="${disable} --reduction_config=${config} --reduction_style=${red_style} --print_reduction_info --initial_filter --subgraph_node_limit=${sg} --time_limit=${t} --seed=${seed} ${database}/${graph}.graph"

                ########## exact with gnn filtering
                # id=${reductions[k]}/${graph}/sg${sg}/seed${seed}
                # create_run "$id" " --gnn_filter ${overall}"

                ########## exact without gnn filtering
                id=${reductions[k]}/${graph}/sg${sg}/seed${seed}
                create_run "$id" "${overall}"
            done
    done
done
