
set(REDUCTION_INCLUDES
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/app 
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/bounds
        ${CMAKE_CURRENT_SOURCE_DIR}/tools
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel
)
set(REDUCTION_INCLUDES ${REDUCTION_INCLUDES} PARENT_SCOPE)

set(LIB_REDUCTION
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/struction_log.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/reductions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/mwis_finder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/original_struction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/extended_struction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/high_degree_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/bound_reduction.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/key_functions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/kernel/LRConv.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/branch_and_reduce_algorithm.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/bounds/partition_cover.cpp
)
add_library(lib_reduction OBJECT ${LIB_REDUCTION})
target_include_directories(lib_reduction PUBLIC ${REDUCTION_INCLUDES})

add_executable(generate_training_data "app/generate_training_data.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(generate_training_data PRIVATE lib)
target_include_directories(generate_training_data PUBLIC ${REDUCTION_INCLUDES})
add_executable(branch_reduce_convergence "app/branch_reduce.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(branch_reduce_convergence PRIVATE lib)
target_include_directories(branch_reduce_convergence PUBLIC ${REDUCTION_INCLUDES})
add_executable(cyclic_blow_up_convergency "app/cyclic_blow_up_convergency.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(cyclic_blow_up_convergency PRIVATE lib)
target_include_directories(cyclic_blow_up_convergency PUBLIC ${REDUCTION_INCLUDES})
add_executable(struction_checker "app/struction_checker.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(struction_checker PRIVATE lib)
target_include_directories(struction_checker PUBLIC ${REDUCTION_INCLUDES})
add_executable(kernelization "app/kernelization.cpp" $<TARGET_OBJECTS:lib_reduction>)
target_link_libraries(kernelization PRIVATE lib)
target_include_directories(kernelization PUBLIC ${REDUCTION_INCLUDES})

target_compile_definitions(cyclic_blow_up_convergency PRIVATE -D OUTPUT_GRAPH_CONVERGENCE)
target_compile_definitions(branch_reduce_convergence PRIVATE -D OUTPUT_WEIGHT_CONVERGENCE)

install(TARGETS generate_training_data DESTINATION bin)
install(TARGETS branch_reduce_convergence DESTINATION bin)
install(TARGETS cyclic_blow_up_convergency DESTINATION bin)
install(TARGETS struction_checker DESTINATION bin)

