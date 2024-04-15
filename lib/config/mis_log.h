/******************************************************************************
 * mis_log.h
 *
 *****************************************************************************/

#ifndef _MIS_LOG_H_
#define _MIS_LOG_H_

#include <sstream>
#include <string>

#include "definitions.h"
#include "timer.h"
#include "config.h"
#include "data_structure/graph_access.h"

class mis_log {
    public:
        /**
         * Get the singleton logger instance.
         * 
         * @return Instance of the logger.
         */
        static mis_log *instance() {
            static mis_log inst;
            return &inst;
        };

        /**
         * Set the config.
         *
         * @param config Config for the evolutionary algorithm.
         */
        void set_config(Config & config);

        /**
         * Set the graph.
         *
         * @param G Graph representation.
         */
        void set_graph(graph_access & G, std::string & name);

        /**
         * Write the log to a file.
         */
        void write_log();

        /**
         * Add a newline to the log.
         */
        void print_newline();

        /**
         * Print information about the graph.
         */
        void print_graph();

        /**
         * Print the current config.
         */
        virtual void print_config();

        /**
         * Print banner for the initialization phase.
         */
        void print_init_title();

        /**
         * Print the final results.
         */
        void print_results();

        /**
         * Print a title.
         */
        virtual void print_title();

        /**
         * Restart the timer for the total time excluding IO etc.
         */
        void restart_total_timer();

        /**
         * Get the timer for the evolutionary algorithm.
         */
        double get_total_timer();

        /**
         * Update the weight of the best solution.
         *
         * @param mis_config Config for the logger.
         * @param weight Candidate to replace the best solution weight.
         */
        void set_best_weight(Config & mis_config, NodeWeight weight);

        /**
         * Reset the size of the best solution.
         */
        void reset_best_weight();

        int get_number_of_nodes() {return number_of_nodes;}
        NodeWeight get_best_weight();

    protected:
        // General information
        timer total_timer;
        std::stringstream filebuffer_string;
        Config log_config;

        // Graph informations
        std::string graph_name;
        NodeID number_of_nodes;
        EdgeID number_of_edges;
        unsigned int arc_scans;
        double avg_degree;
        double density;
        

        // Results information
        double total_time_taken;
        double time_taken_best;
        NodeWeight optimum_weight;
        NodeWeight best_solution_weight;

        /**
         * Default Constructor.
         */
        mis_log();

        /**
         * Default Destructor.
         */
        virtual ~mis_log();

};

#endif
