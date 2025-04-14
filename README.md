# LearnAndReduce 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


This is the project LearnAndReduce. Given a graph G=(V,E,w), the goal is to compute a maximum weight independent set which is NP-hard. This project provides an exact GNN guided preprocessing to reduce input instances for this problem.

## Installation 
As a first step, compile the source by running 
```
./compile_all.sh
```
The binaries can then be found in the folder **deploy**. To compile the programs you need **g++** and **cmake** installed.

The framework contains a graph checking tool to make life a little bit easier:
* graphchecker -- check if your graph file is in the correct format

## Usage LearnAndReduce 
```
./reduce FILE [options]
```

### Options 

| Option | Decription | Default | Mandatory
|-|-|-|-
|`FILE`| Path to graph file that you want the reduce. || &check;
|`--help`| Print help. ||
|`--verbose`|Print detailed information. ||
|`--seed=<int>` |Set seed. | 0 ||
|`--reduction_config=<string>` |Choose reduction configuration: all_reductions, no_gnn_reductions. | all_reductions | |
|`--gnn_filter=<string>` |Choose gnn filtering: never, initial, initial_tight, always. | initial_tight||
|`--time_limit=<double>` |Set time limit (in seconds). | 1000 s||
|`--cyclicFast` | Set CyclicFast configuration. | &check; ||
|`--cyclicStrong` | Set CyclicStrong configuration. | ||
|`--kernel=<string>` | Path to store reduced instance. | ||


### Example 
An example to use the different options is:
```
./deploy/reduce [instance] --gnn_filter=always --cyclicStrong --time_limit=3000
```

### Output

The output of the program without the **-verbose** option is a single line on the format
```
instance_name,struction_config,reduction_config,seed,gnn_filter,#vertices,#edges,#reduced_instance_vertices,#reduced_instance_edges,offset,reduction_time
```


### Input Format

LearnAndReduce expects graphs on the METIS graph format. A graph with **N** vertices is stored using **N + 1** lines. The first line lists the number of vertices, the number of edges, and the weight type. For weighted instances, the first line should use 10 as the weight type to indicate vertex weights. Each subsequent line first gives the weight and then lists the neighbors of that node in **sorted** order.

Here is an example of a graph with 3 vertices of weight 15, 15, and 20, where the weight 20 vertex is connected to the two 15-weight vertices.

```
3 2 10
15 3
15 3
20 1 2
```
Notice that vertices are 1-indexed, and edges appear in the neighborhoods of both endpoints.
You can use the provided graphchecker tool to check if the format of your file is correct.
```
./graphchecker FILE
```   

## Usage Generate Training Data 
The data used for training our GNN models is availiable [here](https://zenodo.org/records/15210077). You can also create your own training data with our tool.
```
./generate_training_data FILE
``` 

The program reads a metis file, generated the training data for the full graph and the reduced instance with training data on the reduced instace stored in training_data/csv.

# Reproducing Results of the Paper

An archive version of the paper is available on arXiv combined with the new local search [CHILS](https://github.com/KennethLangedal/CHILS). The LearnAndReduce part of this combined paper was accepted at ACDA 2025 and will be published in the conference proceedings later in 2025.
```
@article{grossmann2024accelerating,
  title    = {Accelerating Reductions Using Graph Neural Networks and a New Concurrent Local Search for the Maximum Weight Independent Set Problem},
  author   = {Gro{\ss}mann, Ernestine and Langedal, Kenneth and Schulz, Christian},
  journal  = {arXiv preprint arXiv:2412.14198},
  year     = {2024}
}
```

The default setting gives the results in Table 5 of the paper:
```
./deploy/reduce [instance] 
```
