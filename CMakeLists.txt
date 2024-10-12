cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(MMWIS C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
# if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

#probably set by Release, but doing anyway.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DNDEBUG)
    add_compile_options(-O3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-omit-frame-pointer -march=native -DNDEBUG -Wall -funroll-loops -fno-stack-limit -fPIC -O3 -g")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DEBUG)
    add_compile_options(-g) 
    # add_compile_options(-fsanitize=address)
endif()

 if(CMAKE_BUILD_TYPE STREQUAL "Profile")
    add_compile_options(-g) 
    add_compile_options(-pg) 
    set(DCMAKE_CXX_FLAGS -pg)
    set(DCMAKE_EXE_LINKER_FLAGS -pg)
    set(DCMAKE_SHARED_LINKER_FLAGS -pg)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-omit-frame-pointer -march=native -DNDEBUG -Wall -funroll-loops -fno-stack-limit -fPIC -O3 -g")
 endif()


add_definitions(-w)
# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
#CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
#if(COMPILER_SUPPORTS_WALL)
  #add_definitions(-Wall)
#endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSED)
if(COMPILER_SUPPORTS_NOUNUSED)
add_definitions(-Wno-unused-value)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-unused-value COMPILER_SUPPORTS_NOUNUSEDRES)
if(COMPILER_SUPPORTS_NOUNUSEDRES)
add_definitions(-Wno-unused-result)
endif()

# check dependencies
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit options
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()
option(64BITMODENODEWEIGHTS "64 bit mode (node weights)" ON)
if(64BITMODENODEWEIGHTS)
  add_definitions("-DMODE64BITNODEWEIGHTS")
  add_definitions("-DPOINTER64=1")
endif()

set(TSOLVER ${CMAKE_CURRENT_SOURCE_DIR}/solver/tiny_solver)
set(BSOLVER ${CMAKE_CURRENT_SOURCE_DIR}/solver/reducer)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/app
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/io
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/config
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/tools
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure/priority_queues
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/basics 
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.2.2
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TSOLVER}
    ${BSOLVER}
    ${BSOLVER}/app 
    ${BSOLVER}/lib
    ${BSOLVER}/lib/bounds
    ${BSOLVER}/tools
    ${BSOLVER}/kernel
    ${BSOLVER}/kernel/simple_reduction_rules
    ${BSOLVER}/kernel/expensive_reduction_rules
    ${BSOLVER}/kernel/struction_reduction_rules
    ${BSOLVER}/kernel/heuristic_reduction_rules
    ${BSOLVER}/kernel/csr_reductions
)


set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure/csr_graph.c
    ${TSOLVER}/tiny_solver.c
    ${BSOLVER}/kernel/csr_reductions/src/reductions.c
    PROPERTIES LANGUAGE C
)

set(GENERAL
    lib/io/graph_io.cpp
    lib/algorithms/push_relabel.cpp
    lib/data_structure/operation_log.cpp
    lib/data_structure/priority_queues/bucket_array.cpp
    lib/data_structure/csr_graph.c
    lib/tools/random_functions.cpp
    lib/config/mis_log.cpp
    lib/config/reduction_log.cpp
    lib/basics/cout_handler.cpp
    extern/argtable3-3.2.2/argtable3.c
    ${TSOLVER}/tiny_solver.c
    ${BSOLVER}/tools/struction_log.cpp
    ${BSOLVER}/lib/reduce_algorithm.cpp
    ${BSOLVER}/kernel/general_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/fold1_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/fold2_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/neighborhood_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/clique_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/clique_neighborhood_fast_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/twin_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/domination_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/extended_domination_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/extended_domination_reverse_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/extended_twin_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/single_edge_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/extended_single_edge_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/funnel_reduction.cpp
    ${BSOLVER}/kernel/simple_reduction_rules/funnel_fold_reduction.cpp
    ${BSOLVER}/kernel/struction_reduction_rules/mwis_finder.cpp
    ${BSOLVER}/kernel/struction_reduction_rules/original_struction.cpp
    ${BSOLVER}/kernel/struction_reduction_rules/extended_struction.cpp
    ${BSOLVER}/kernel/struction_reduction_rules/key_functions.cpp
    ${BSOLVER}/kernel/struction_reduction_rules/struction_reductions.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/clique_neighborhood_reduction.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/critical_set_reduction.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/cut_vertex_reduction.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/heavy_set_reduction.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/heavy_set3_reduction.cpp
    ${BSOLVER}/kernel/expensive_reduction_rules/generalized_fold_reduction.cpp
    ${BSOLVER}/kernel/heuristic_reduction_rules/heuristic_include_reduction.cpp
    ${BSOLVER}/kernel/heuristic_reduction_rules/heuristic_exclude_reduction.cpp
    ${BSOLVER}/kernel/csr_reductions/src/reductions.c
    ${BSOLVER}/kernel/csr_reductions/unconfined_csr_reduction.cpp
    ${BSOLVER}/kernel/LRConv.cpp
    ${BSOLVER}/kernel/lr_gcn.cpp
    ${BSOLVER}/lib/bounds.cpp
)
add_library(lib OBJECT ${GENERAL})
target_link_libraries(lib PRIVATE OpenMP::OpenMP_C)

add_executable(graphchecker app/graphchecker.cpp)
target_link_libraries(graphchecker ${OpenMP_CXX_LIBRARIES})

add_executable(metis_to_cosmo "app/metis_to_cosmo.cpp")
target_link_libraries(metis_to_cosmo PRIVATE lib)

add_executable(generate_full_graph_training_data "app/generate_full_graph_training_data.cpp")
target_link_libraries(generate_full_graph_training_data PRIVATE lib)

add_executable(kernelization "app/kernelization.cpp")
target_link_libraries(kernelization PRIVATE lib)
