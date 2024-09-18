#!/bin/bash

# Initialize the CSV file with headers
echo "id,config,graph,kernelV,seed,offset,reduction_time,time,weight,optimal,maxComponent" > results.csv

# Find all result.txt files and process them
find . -iname result.txt | while read i; do

    # Extract data from result.txt in a single read
    data=`cat $i | grep -E 'Config:|Filename:|Seed:|Kernel nodes:|Offset:|MaxComponent:|optimal:|Time|Final'`

    id=$(echo "$i" | awk -F "/" '{print $2}')
    graph=$(echo "$data" | grep 'Filename:' | awk '{print $NF}')
    config=$(echo "$data" | grep 'Config:' | awk '{print $NF}')
    seed=$(echo "$data" | grep 'Seed:' | awk '{print $NF}')
    kernelV=$(echo "$data" | grep 'Kernel nodes:' | awk '{print $NF}')
    offset=$(echo "$data" | grep 'Offset:' | awk '{print $NF}')
    reduction_time=$(echo "$data" | grep 'Time:' | awk 'END{print $NF}')
    maxComponent=$(echo "$data" | grep 'MaxComponent:' | awk '{print $NF}')
    optimal=$(echo "$data" | grep 'optimal:' | awk '{print $NF}')
    time=$(echo "$data" | grep 'Time found:' | awk '{print $NF}')
    weight=$(echo "$data" | grep 'Final Weight:' | awk '{print $NF}')

    # Append the results to the CSV file
    echo "${id},${config},${graph},${kernelV},${seed},${offset},${reduction_time},${time},${weight},${optimal},${maxComponent}" >> results.csv
done

