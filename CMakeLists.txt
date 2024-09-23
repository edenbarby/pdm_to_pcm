cmake_minimum_required(VERSION 3.14)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

set(PROJECT_LIB pdm_to_pcm)
set(PROJECT_TEST pdm_to_pcm_test)

project(${PROJECT_LIB})

add_executable(
    ${PROJECT_TEST}
    src/PdmToPcm.cpp
    test/PdmToPcmTests.cpp
    test/ChebyshevPolynomialTests.cpp
    test/KaiserWindowTests.cpp
)

target_include_directories(
    ${PROJECT_TEST}
    PRIVATE
    src
    inc
    gcem/include # gcem is a header only library so instead of linking we just include the headers
)

target_link_libraries(
    ${PROJECT_TEST}
    GTest::gtest_main
)

add_subdirectory(googletest)
include(GoogleTest)
enable_testing()
gtest_discover_tests(${PROJECT_TEST})