
set(milestone5_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        )


add_executable(milestone05 ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${milestone5_SRCS})
include_directories(milestone05 ${CMAKE_SOURCE_DIR})
target_link_libraries(milestone05 PUBLIC Eigen3::Eigen ${MPI_CXX_LIBRARIES})

