# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt" AND EXISTS "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitinfo.txt" AND
  "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt" IS_NEWER_THAN "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external/Eigen3"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external/Eigen3'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git" 
            clone --no-checkout --config "advice.detachedHead=false" "https://gitlab.com/libeigen/eigen.git" "Eigen3"
    WORKING_DIRECTORY "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://gitlab.com/libeigen/eigen.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git" 
          checkout "3.4.0" --
  WORKING_DIRECTORY "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external/Eigen3"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '3.4.0'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external/Eigen3"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/external/Eigen3'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitinfo.txt" "/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/eslam/Desktop/Molecular-Dynamics/cmake-build-debug/third-party/Eigen3-download/Eigen3-prefix/src/Eigen3-stamp/Eigen3-gitclone-lastrun.txt'")
endif()
