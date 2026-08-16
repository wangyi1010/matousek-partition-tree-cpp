cmake_minimum_required(VERSION 3.20)
project(matousek_partition_tree VERSION 1.0.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)

option(MPT_BUILD_TESTS "Build property tests" ON)
option(MPT_BUILD_BENCHMARKS "Build benchmark programs" ON)

find_path(BOOST_INCLUDE_DIR boost/rational.hpp
    HINTS /opt/homebrew/include /usr/local/include
    DOC "Directory containing Boost headers")
if(NOT BOOST_INCLUDE_DIR)
    message(FATAL_ERROR "Boost headers are required (boost/rational.hpp and boost/multiprecision/cpp_int.hpp)")
endif()

add_library(matousek_partition_tree src/core.cpp)
add_library(MatousekPartitionTree::matousek_partition_tree ALIAS matousek_partition_tree)
target_compile_features(matousek_partition_tree PUBLIC cxx_std_20)
target_include_directories(matousek_partition_tree
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_include_directories(matousek_partition_tree SYSTEM PUBLIC
    $<BUILD_INTERFACE:${BOOST_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${BOOST_INCLUDE_DIR}>)

if(MSVC)
    target_compile_options(matousek_partition_tree PRIVATE /W4 /permissive-)
else()
    target_compile_options(matousek_partition_tree PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(matousek-demo apps/demo.cpp)
target_link_libraries(matousek-demo PRIVATE MatousekPartitionTree::matousek_partition_tree)

add_executable(matousek-visualize apps/visualize.cpp)
target_link_libraries(matousek-visualize PRIVATE MatousekPartitionTree::matousek_partition_tree)

if(MPT_BUILD_BENCHMARKS)
    add_executable(measure-crossings benchmarks/measure_crossings.cpp)
    target_link_libraries(measure-crossings PRIVATE MatousekPartitionTree::matousek_partition_tree)
    add_executable(measure-performance benchmarks/measure_performance.cpp)
    target_link_libraries(measure-performance PRIVATE MatousekPartitionTree::matousek_partition_tree)
endif()

include(GNUInstallDirs)
install(TARGETS matousek_partition_tree EXPORT MatousekPartitionTreeTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT MatousekPartitionTreeTargets
    NAMESPACE MatousekPartitionTree::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatousekPartitionTree)
configure_package_config_file(
    cmake/MatousekPartitionTreeConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/MatousekPartitionTreeConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatousekPartitionTree)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/MatousekPartitionTreeConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/MatousekPartitionTreeConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/MatousekPartitionTreeConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MatousekPartitionTree)

if(MPT_BUILD_TESTS)
    enable_testing()
    add_executable(matousek-tests tests/test_matousek_partition_tree.cpp)
    target_link_libraries(matousek-tests PRIVATE MatousekPartitionTree::matousek_partition_tree)
    add_test(NAME property-tests COMMAND matousek-tests)
endif()
