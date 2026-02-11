# LearnAndReduce 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/786850375.svg)](https://doi.org/10.5281/zenodo.15229974)


This is the project LearnAndReduce. Given a graph G=(V,E,w), the goal is to compute a maximum weight independent set which is NP-hard. This project provides an exact GNN guided preprocessing to reduce input instances for this problem.
## Dependencies
For the local search CHILS, the get_dep.sh script clones and builds the necessary library.
```
./get_dep.sh
```
## Installation 
Once the CHILS library is build, compile the project by running 
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

| Option | Decription | Default
|-|-|-
|`--help`| Print help. |
|`--verbose`|Print detailed information. |
|`--seed=<int>` |Set seed. | 0 |
|`--reduction_config=<string>` |Choose reduction configuration: all_reductions, no_gnn_reductions. | all_reductions |
|`--gnn_filter=<string>` |Choose gnn filtering: never, initial, initial_tight, always. | initial_tight|
|`--time_limit=<double>` |Set time limit (in seconds). | 1000 s|
|`--cyclicFast` | Set CyclicFast configuration. | &check; |
|`--cyclicStrong` | Set CyclicStrong configuration. | |
|`--kernel=<string>` | Path to store reduced instance. | |
|`--solution_from_file` | Option to lift your own solution on the reduced instance. After the reduction, you are asked to provide a path to that solution. | |
|`--chils_time_limit=<double>` | Set time limit (in seconds) for running the local search CHILS to compute a solution on the reduced instance. | 600 s|
|`--chils_n_solutions=<int>` | Set number of solutions used for running the local search CHILS to compute a solution on the reduced instance. | 16 |


### Example 
An example to use the different options is:
```
./deploy/reduce [instance] --gnn_filter=always --cyclicStrong --time_limit=3000
```

### Output

The output of the program without the **--verbose** option is a single line on the format
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

### Solution Format
The solution file for a graph with $n$ vertices contains $n$ lines, one for each vertex. If a vertex is in the solution, the line contains a 1; if it is not, it contains a 0.


## Usage Generate Training Data 
The data used for training our GNN models is availiable [here](https://zenodo.org/records/15210077). You can also create your own training data with our tool.
```
./generate_training_data FILE
``` 

The program reads a metis file and then generates the training data for the full and reduced graphs. The training data is stored in training_data/csv.

# Reproducing Results of the Paper

When using or comparing against this method, please cite the following paper.
```
@inproceedings{grossmann2025accelerating,
  title={Accelerating Reductions Using Graph Neural Networks for the Maximum Weight Independent Set Problem},
  author={Gro{\ss}mann, Ernestine and Langedal, Kenneth and Schulz, Christian},
  booktitle={2025 Proceedings of the Conference on Applied and Computational Discrete Algorithms (ACDA)},
  pages={155--168},
  year={2025},
  doi={https://doi.org/10.1137/1.9781611978759.12},
  organization={SIAM}
}

```

The default setting gives the results in Table 5 of the paper:
```
./deploy/reduce [instance] 
```
