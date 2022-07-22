add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
configure_file(lj54.xyz .)

set(milestone4_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )


add_executable(milestone04 ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${milestone4_SRCS})
include_directories(milestone04 ${CMAKE_SOURCE_DIR})
target_link_libraries(milestone04 PUBLIC Eigen3::Eigen gtest gtest_main ${MPI_CXX_LIBRARIES})
