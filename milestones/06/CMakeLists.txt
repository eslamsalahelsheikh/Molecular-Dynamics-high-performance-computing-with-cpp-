
set(milestone6_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )


add_executable(milestone06 ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${milestone6_SRCS})
include_directories(milestone06 ${CMAKE_SOURCE_DIR})
target_link_libraries(milestone06 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})

