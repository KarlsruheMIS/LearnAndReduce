## Description ##
This is the project LearnAndReduce. Given a graph G=(V,E), the goal of the maximum independent set problem is to compute a maximum cardinality set of vertices I, such that no vertices in the set are adjacent to one another. Such a set is called a maximum independent set. The problem is NP-hard and particularly difficult to solve in large sparse graphs. 
This project provides a GNN guided Preprocessing to reduce input instances for this problem.

## Installation ##
As a first step, compile the source by running *compile_withcmake.sh*. The binaries can then be found in the folder *deploy*.  To compile the programs you need g++, OpenMP and cmake installed. 

The framework contains a graph checking tool to make life a little bit easier:
* graphchecker -- check if the graph file you gave to algorithm is in the correct format

## Usage LearnAndReduce ##
`kernelization FILE [options]`.    


### Options ###
This is a brief overview of the most important options.
For a full description, please take a look at the user guide.

`FILE`
Path to graph file that you want the reduce.

`--help`
Print help.

`--console_log`
Write the log to the console.

`--verbose`
Print detailed information.

`--seed=<int>`
Seed to use for the random number generator.

`--reduction_config=<string>`
Config to use for the reduction configuration [all_reductions_cyclicFast|all_reductions_cyclicStrong|no_gnn_reductions_cyclcicFast|no_gnn_reductions_cyclicStrong].

`--gnn_filter=<string>`
Config to use for the gnn filtering [never(default)|initial|initial_tight|always].

`--time_limit=<double>`
Time limit until the algorithm terminates.


## Usage generate training data ##
`generate_training_data FILE`.    

The program reads a Metis file, generated the training data for the full graph and the reduced instance with training data on the reduced instace stored in training_data/csv.

## Usage Graph Checker ##
`graphchecker FILE`.    

The program reads a Metis file and checks the file for correctness.
