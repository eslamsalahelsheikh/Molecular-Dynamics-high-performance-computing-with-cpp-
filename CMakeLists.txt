cmake_minimum_required(VERSION 3.11)

project(MolecularDynamics)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

find_package(MPI)
if(NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is needed!")
endif(NOT MPI_FOUND)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
add_definitions(-DCOSMO_MPI)
if(NOT MPI_CXX_LIBRARIES)
    message(FATAL_ERROR "MPI found but the MPI C++ libraries not found!")
endif(NOT MPI_CXX_LIBRARIES)

include(tools)
add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(MolecularDynamics_HDRS
        # Add header files here
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/verlet.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/types.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/atoms.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/lj_direct_summation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/xyz.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/energy.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/neighbors.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/gupta.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/simulation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/simulation_data.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/utils.h
        ${CMAKE_CURRENT_SOURCE_DIR}/include/Molecular-Dynamics/vector.h
        )


set(MolecularDynamics_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/src/verlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/lj_direct_summation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/energy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gupta.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation_data.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ih.cpp
        )

add_subdirectory(milestones)
add_subdirectory(tests)

