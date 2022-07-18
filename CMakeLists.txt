cmake_minimum_required(VERSION 3.11)

project(molecular_dynamics)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)

set(MYPROJECT_HDRS
        # Add header files here
        include/Molecular-Dynamics/verlet.h
        include/Molecular-Dynamics/types.h
        include/Molecular-Dynamics/atoms.h
        include/Molecular-Dynamics/lj_direct_summation.h include/Molecular-Dynamics/xyz.h include/Molecular-Dynamics/energy.h)

set(MYPROJECT_SRCS
        # Add source files here
        src/main.cpp
        src/verlet.cpp
        tests/test_lj_direct_summation.cpp
        src/lj_direct_summation.cpp src/xyz.cpp src/energy.cpp)

set(MYPROJECT_executable_SRCS
        # Add executable source files here
        src/main.cpp)


# Simulation code
add_executable(molecular_dynamics ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${MYPROJECT_executable_SRCS})
target_link_libraries(molecular_dynamics PUBLIC Eigen3::Eigen)
