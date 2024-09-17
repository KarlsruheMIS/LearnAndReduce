/**
 * config.h
 * Purpose: Configuration used for different algorithms.
 *
 *****************************************************************************/

#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <string>
#include <cctype>
/* #include <algorithm> */

#include "definitions.h"

// Configuration for the calculation 
struct Config {
    enum Weight_Source {FILE, HYBRID, SMALL_UNIFORM, UNIFORM, GEOMETRIC, UNIT};
    enum Key_Type {RANDOM, DEGREE, INCREASE, APPROXIMATE_INCREASE};

    // Name of the graph file.
    std::string graph_filename;
    // Name of the output file.
    std::string output_filename;
    // Seed for the RNG.
    int seed;
    // Time limit for the total algorithm
    double time_limit;
    // Force parameter for the mutation.
    unsigned int force_k;
    // Number of candidates for forced insertion.
    unsigned int force_cand;
    // Threshold for local search 
    unsigned int local_search_threshold;
    // Number of swaps in ils pruning 
    unsigned int max_swaps = 1000;
    // Time limit for the ils 
    /* double ils_time_limit; */
    // Print detailed information 
    bool verbose = false;
    // Write the log into a file
    bool print_log;
    // Write the inpendent set into a file
    bool write_graph;
    // Write the log into the console
    bool console_log;
    // Number of iterations for the ILS.
    unsigned int ils_iterations;
    // Optimize candidates for ILS.
    bool optimize_candidates;
    // Start Greedy adaptive (initial solution)
    bool start_greedy_adaptive;
    // Sort free nodes in local search
    bool sort_freenodes;
    // Check graph sortedness
    bool check_sorted;
    // Source for weights
    Weight_Source weight_source;

    void setWeightSource(const std::string & source) {
        if (strCompare(source, "file")) {
            weight_source = Weight_Source::FILE;
        } else if (strCompare(source, "hybrid")) {
            weight_source = Weight_Source::HYBRID;
        } else if (strCompare(source, "unit")) {
            weight_source = Weight_Source::UNIT;
        } else if (strCompare(source, "uniform")) {
            weight_source = Weight_Source::UNIFORM;
        } else if (strCompare(source, "small_uniform")) {
            weight_source = Weight_Source::SMALL_UNIFORM;
        } else if (strCompare(source, "geometric")) {
            weight_source = Weight_Source::GEOMETRIC;
        }
    }


    protected:

    bool strCompare(const std::string & str1, const std::string & str2) {
        return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](unsigned char c1, unsigned char c2){ return std::toupper(c1) == std::toupper(c2); });
    }
};

#endif
