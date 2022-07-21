cmake_minimum_required(VERSION 3.11)

project(MolecularDynamics)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

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
        )

set(MolecularDynamics_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/src/verlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/lj_direct_summation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/energy.cpp)

set(MolecularDynamics_executable_SRCS
        # Add executable source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_subdirectory(tests)

# Simulation code
add_executable(moleculardynamics ${MolecularDynamics_SRCS} ${MolecularDynamics_HDRS} ${MolecularDynamics_executable_SRCS})
target_link_libraries(moleculardynamics PUBLIC Eigen3::Eigen)