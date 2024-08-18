
set(REDUCTION_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/app 
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/bounds
        ${CMAKE_CURRENT_SOURCE_DIR}/tools
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/heuristic_reduction_rules
)
set(REDUCTION_INCLUDES ${REDUCTION_INCLUDES} PARENT_SCOPE)

set(LIB_REDUCTION
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/struction_log.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/general_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/fold1_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/fold2_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/neighborhood_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/clique_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/clique_neighborhood_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/clique_neighborhood_fast_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/twin_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/domination_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/extended_domination_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/single_edge_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/extended_single_edge_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/funnel_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/simple_reduction_rules/funnel_fold_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules/mwis_finder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules/original_struction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules/extended_struction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules/key_functions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/struction_reduction_rules/struction_reductions.cpp
        # ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/high_degree_reduction.cpp
        # ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/bound_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/critical_set_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/cut_vertex_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/heavy_set_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/heavy_set3_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/expensive_reduction_rules/generalized_fold_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/heuristic_reduction_rules/heuristic_include_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/heuristic_reduction_rules/heuristic_exclude_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/LRConv.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/bounds.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/branch_and_reduce_algorithm.cpp
        # ${CMAKE_CURRENT_SOURCE_DIR}/lib/bounds/partition_cover.cpp
)
add_library(lib_reduction OBJECT ${LIB_REDUCTION})
target_include_directories(lib_reduction PUBLIC ${REDUCTION_INCLUDES})

add_executable(generate_training_data "app/generate_training_data.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(generate_training_data PRIVATE lib)
target_include_directories(generate_training_data PUBLIC ${REDUCTION_INCLUDES})
add_executable(branch_reduce_convergence "app/branch_reduce.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(branch_reduce_convergence PRIVATE lib)
target_include_directories(branch_reduce_convergence PUBLIC ${REDUCTION_INCLUDES})
add_executable(kernelization "app/kernelization.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(kernelization PRIVATE lib)
target_include_directories(kernelization PUBLIC ${REDUCTION_INCLUDES})

install(TARGETS generate_training_data DESTINATION bin)
install(TARGETS branch_reduce_convergence DESTINATION bin)
install(TARGETS kernelization DESTINATION bin)

