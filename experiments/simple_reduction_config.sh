#!/bin/bash

run="simple_config"
program="kernelization"
instance="complete_set"
run_parallel="run_parallel_${run}.txt"
address="results/${run}/${instance}"
database="/home/ernestineg/test_instances/${instance}"

path_to_program="../deploy/${program}"
result_file="results/results_${run}_${instance}.csv"
error_file="results/error_${run}_${instance}.txt"
mkdir -p $address

create_run() {
    id=$1
    id_instructions=$2
    address2=${address}/${id}
    identify=${instance}/${program}/${id}
    mkdir -p ${address2}
    first=`echo ${id} | awk -F  "/" '{print $1}'`
    
    instruction="echo -n '${first},' > ${address2}/result.txt && ${path_to_program} ${id_instructions}"
    echo " ${instruction} >> ${address2}/result.txt || { echo '${id} ${address2}' >> ${error_file}; echo 'error,' >> ${result_file}; } &&  cat ${address2}/result.txt >> ${result_file}; " >> ${run_parallel}
}
rm -f ${run_parallel}

echo "id,graph,redstyle,redconfig,timelimit,seed,gnnfilter,n,m,kn,km,offset,time" > ${result_file}

    t=36000
    sg=128
            #########################################################
            ############### Experiments for soa  ####################
            #########################################################

for seed in 0; do # seeds
    for graph in `ls ${database} | grep "graph" | awk -F"." '{print $1}'`; do
        for config in  all_decreasing; do
            overall=" --kernel=${kernel_path}/${graph} --disable_generalized_fold --disable_unconfined --disable_funnel --disable_funnel_fold --reduction_config=${config} --disable_blow_up --disable_plateau_struction --disable_decreasing_struction --disable_cut_vertex --disable_heavy_set --disable_heavy_set3 --disable_extended_twin --disable_extended_domination --disable_extended_domination_reverse --disable_critical_set --time_limit=${t} --seed=${seed} ${database}/${graph}.graph"

            id="simple/${graph}"
            create_run "$id" "${overall}"
        done
    done
done
