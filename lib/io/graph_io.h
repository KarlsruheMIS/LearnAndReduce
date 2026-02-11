/******************************************************************************
 * graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPHIO_H_
#define GRAPHIO_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "definitions.h"
#include "data_structure/graph_access.h"

class graph_io {
        public:
                graph_io();
                virtual ~graph_io () ;

                static
                int readGraphWeighted(graph_access & G, const std::string & filename);

                static
                int writeGraphWeighted(graph_access & G, const std::string & filename);

                static
                int writeGraph(graph_access & G, const std::string & filename);

                static
                void writeMWIS(graph_access& G, const std::string & filename);

                template<typename vectortype>
                static void writeVector(std::vector<vectortype> & vec, const std::string & filename);

                template<typename vectortype>
                static bool readVector(std::vector<vectortype> & vec, const std::string & filename);


};

template<typename vectortype>
void graph_io::writeVector(std::vector<vectortype> & vec, const std::string & filename) {
        std::ofstream f(filename.c_str());
        for( unsigned i = 0; i < vec.size(); ++i) {
                f << vec[i] <<  std::endl;
        }

        f.close();
}

template<typename vectortype>
bool graph_io::readVector(std::vector<vectortype> & vec, const std::string & filename) {

        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening file " << filename << std::endl;
                return false;
        }

        unsigned pos = 0;
        std::getline(in, line);
        while( !in.eof() ) {
                if (line[0] == '%') { //Comment
                        continue;
                }

                vectortype value = (vectortype) atof(line.c_str());
                vec[pos++] = value;
                std::getline(in, line);
        }

        in.close();
        return true;
}

#endif /*GRAPHIO_H_*/
