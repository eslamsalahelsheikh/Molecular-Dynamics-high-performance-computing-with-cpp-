configure_file(../../clusters/cluster_923.xyz .)
configure_file(../../clusters/cluster_3871.xyz .)

set(milestone7_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )


add_executable(milestone07 ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${milestone7_SRCS})
include_directories(milestone07 ${CMAKE_SOURCE_DIR})
target_link_libraries(milestone07 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})

