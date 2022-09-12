configure_file(./cluster_923.xyz .)

set(milestone8_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )


add_executable(milestone08 ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${milestone8_SRCS})
include_directories(milestone08 ${CMAKE_SOURCE_DIR})
target_link_libraries(milestone08 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES} ${MPI_INCLUDE_PATH})
