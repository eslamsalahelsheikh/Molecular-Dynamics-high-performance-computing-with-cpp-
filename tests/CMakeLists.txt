include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(myproject_tests_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/test_hello_world.cpp
        )

# Tests
include_directories(myproject_tests ${CMAKE_SOURCE_DIR})
add_executable(myproject_tests ${MYPROJECT_SRCS} ${MYPROJECT_HDRS} ${myproject_tests_SRCS})
target_link_libraries(myproject_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(myproject_tests)

add_test(NAME myproject_tests COMMAND myproject_tests)
