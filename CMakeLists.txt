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
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/app
    ${CMAKE_CURRENT_SOURCE_DIR}/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/config
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure/priority_queues
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/basics 
    ${CMAKE_CURRENT_SOURCE_DIR}/solver
    ${CMAKE_CURRENT_SOURCE_DIR}/solver/initial
    ${CMAKE_CURRENT_SOURCE_DIR}/solver/ils
    ${CMAKE_CURRENT_SOURCE_DIR}/solver/hils
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/io
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/tools
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/partition
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/lib/data_structure
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaHIP/extern/argtable3-3.2.2
)

add_subdirectory(extern/KaHIP EXCLUDE_FROM_ALL)


set(GENERAL
    lib/data_structure/operation_log.cpp
    lib/data_structure/priority_queues/bucket_array.cpp
    lib/data_structure/mis_permutation.cpp
    lib/data_structure/candidate_list.cpp
    lib/config/mis_log.cpp
    lib/config/reduction_log.cpp
    lib/basics/cout_handler.cpp
    )
add_library(lib_general OBJECT ${GENERAL})

set(INITIAL_SOLVER
    solver/initial/initial_mis.cpp
    solver/initial/greedy_mis.cpp
    solver/initial/greedy_mwis.cpp
    solver/initial/random_mis.cpp
    solver/initial/greedy_mwis.cpp
    solver/initial/greedy_mis.cpp
    solver/initial/greedy_weighted_vertex.cpp
    solver/initial/greedy_vertex.cpp
    )
add_library(lib_initial_solver OBJECT ${INITIAL_SOLVER})


# build branch and reduce solvers 
add_subdirectory(solver/branch_reduce)

set(LOCAL_SEARCH
    solver/ils/ils.cpp
    solver/ils/local_search.cpp
    solver/hils/Solution.cpp
    )
add_library(lib_local_search OBJECT ${LOCAL_SEARCH})
target_include_directories(lib_local_search PUBLIC ${REDUCTION_INCLUDES})
target_link_libraries(lib_local_search PUBLIC lib_reduction)

add_library(lib STATIC $<TARGET_OBJECTS:libkaffpa>
                       $<TARGET_OBJECTS:libmapping>
                       $<TARGET_OBJECTS:libnodeordering>
                       $<TARGET_OBJECTS:lib_general>
                       $<TARGET_OBJECTS:lib_initial_solver>
                       $<TARGET_OBJECTS:lib_local_search>
                        )


add_executable(metis_to_cosmo "app/metis_to_cosmo.cpp" $<TARGET_OBJECTS:lib_general>)
target_link_libraries(metis_to_cosmo PRIVATE lib)
target_include_directories(metis_to_cosmo PUBLIC ${GENERAL_INCLUDES})
