# Learn and Reduce #

This is the project LearnAndReduce. Given a graph G=(V,E,w), the goal is to compute a maximum weight independent set which is NP-hard. This project provides a GNN guided preprocessing to reduce input instances for this problem.

## Installation ##
As a first step, compile the source by running `compile_all.sh`. The binaries can then be found in the folder **deploy**. To compile the programs you need **g++** and **cmake** installed.

The framework contains a graph checking tool to make life a little bit easier:
* graphchecker -- check if the graph file you gave to algorithm is in the correct format

The instances used for our evaluation can be downloaded from Dropbox at this [link](https://www.dropbox.com/scl/fi/kbpttzi2woiqfhwvgjadi/LearnAndReduceInstances.zip?rlkey=ijl6uz9indkihxc7luv92mzyd&st=bkyu8vea&dl=0).

## Usage LearnAndReduce ##
`reduce FILE [options]`.

### Options ###
This is a brief overview of the most important options. For a full overview, use the ```--help``` option.

`FILE`
Path to graph file that you want the reduce.

`--help`
Print help.

`--verbose`
Print detailed information.

`--seed=<int>`
Seed to use for the random number generator.

`--reduction_config=<string>`
Config to use for the reduction configuration [all_reductions (default) | no_gnn_reductions ].

`--gnn_filter=<string>`
Config to use for the gnn filtering [never | initial | initial_tight (default) | always].

`--time_limit=<double>`
Time limit until the algorithm terminates.

The exact command used for the results in Table 5 of the paper is: \
`./deploy/reduce [instance] --reduction_config=all_reductions --gnn_filter=initial_tight --cylcicFast`

## Usage generate training data ##
`generate_training_data FILE`.    

The program reads a Metis file, generated the training data for the full graph and the reduced instance with training data on the reduced instace stored in training_data/csv.

## Usage Graph Checker ##
`graphchecker FILE`.    

The program reads a Metis file and checks the file for correctness.
