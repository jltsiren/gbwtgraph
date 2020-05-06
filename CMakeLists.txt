# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(gbwtgraph)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -g")

# We use OpenMP for parallelism
# Use openmp for parallelism, but it's configured differently on OSX
find_package(OpenMP)
if (OPENMP_FOUND)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  endif()
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
message("${sdsl-lite_INCLUDE}")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")
set(googletest_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/googletest/googletest/")
ExternalProject_Get_property(sdsl-lite SOURCE_DIR)
set(googletest_INCLUDE "${SOURCE_DIR}/external/googletest/googletest/include")

# mmmultimap (memory mapped multimap)
ExternalProject_Add(mmmulti
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmmulti SOURCE_DIR)
set(mmmulti_INCLUDE "${SOURCE_DIR}/src")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/ips4o"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/BBHash"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# GBWT
ExternalProject_Add(gbwt-src
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gbwt"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gbwt-src SOURCE_DIR)
set(gbwt_SOURCE "${SOURCE_DIR}")
set(gbwt_INCLUDE "${SOURCE_DIR}/include")
set(gbwt_LIB "${SOURCE_DIR}/lib")

add_library(gbwt STATIC
  ${gbwt_SOURCE}/algorithms.cpp
  ${gbwt_SOURCE}/bwtmerge.cpp
  ${gbwt_SOURCE}/cached_gbwt.cpp
  ${gbwt_SOURCE}/dynamic_gbwt.cpp
  ${gbwt_SOURCE}/files.cpp
  ${gbwt_SOURCE}/gbwt.cpp
  ${gbwt_SOURCE}/internal.cpp
  ${gbwt_SOURCE}/metadata.cpp
  ${gbwt_SOURCE}/support.cpp
  ${gbwt_SOURCE}/utils.cpp
  ${gbwt_SOURCE}/variants.cpp)
add_dependencies(gbwt sdsl-lite)
add_dependencies(gbwt gbwt-src)
target_include_directories(gbwt PUBLIC
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${gbwt_INCLUDE}")

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libhandlegraph"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/include")
set(handlegraph_LIB "${INSTALL_DIR}/lib")

set(CMAKE_BUILD_TYPE Debug)

add_library(gbwtgraph STATIC
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp
  ${CMAKE_SOURCE_DIR}/path_cover.cpp)
add_dependencies(gbwtgraph sdsl-lite)
add_dependencies(gbwtgraph mmmulti)
add_dependencies(gbwtgraph ips4o)
add_dependencies(gbwtgraph bbhash)
add_dependencies(gbwtgraph gbwt)
add_dependencies(gbwtgraph handlegraph)
target_include_directories(gbwtgraph PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmulti_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}")

set(gbwtgraph_INCLUDE
  "${CMAKE_SOURCE_DIR}/include"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmulti_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}")
set(gbwtgraph_LIBS
  "${LIBRARY_OUTPUT_PATH}/libgbwt.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "-latomic"
  z)

# set up our target executable and specify its dependencies and includes
add_executable(gfa2gbwt
  ${CMAKE_SOURCE_DIR}/gfa2gbwt.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(gfa2gbwt gbwtgraph)
target_include_directories(gfa2gbwt PUBLIC ${gbwtgraph_INCLUDE})
target_link_libraries(gfa2gbwt ${gbwtgraph_LIBS})

set(gbwtgraph_tests_INCLUDE
  ${gbwtgraph_INCLUDE}
  ${googletest_INCLUDE})

set(gbwtgraph_tests_LIBS
  ${gbwtgraph_LIBS}
  ${googletest_LIB}/libgtest.a
  ${googletest_LIB}/libgtest_main.a)

add_executable(test_algorithms
  ${CMAKE_SOURCE_DIR}/tests/test_algorithms.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_algorithms gbwtgraph)
target_include_directories(test_algorithms PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_algorithms ${gbwtgraph_tests_LIBS})

add_executable(test_cached_gbwtgraph
  ${CMAKE_SOURCE_DIR}/tests/test_cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_cached_gbwtgraph gbwtgraph)
target_include_directories(test_cached_gbwtgraph PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_cached_gbwtgraph ${gbwtgraph_tests_LIBS})

add_executable(test_gbwtgraph
  ${CMAKE_SOURCE_DIR}/tests/test_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_gbwtgraph gbwtgraph)
target_include_directories(test_gbwtgraph PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_gbwtgraph ${gbwtgraph_tests_LIBS})

add_executable(test_gfa
  ${CMAKE_SOURCE_DIR}/tests/test_gfa.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_gfa gbwtgraph)
target_include_directories(test_gfa PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_gfa ${gbwtgraph_tests_LIBS})

add_executable(test_minimizer
  ${CMAKE_SOURCE_DIR}/tests/test_minimizer.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_minimizer gbwtgraph)
target_include_directories(test_minimizer PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_minimizer ${gbwtgraph_tests_LIBS})

add_executable(test_path_cover
  ${CMAKE_SOURCE_DIR}/tests/test_path_cover.cpp
  ${CMAKE_SOURCE_DIR}/algorithms.cpp
  ${CMAKE_SOURCE_DIR}/cached_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp
  ${CMAKE_SOURCE_DIR}/path_cover.cpp)
add_dependencies(test_path_cover gbwtgraph)
target_include_directories(test_path_cover PUBLIC ${gbwtgraph_tests_INCLUDE})
target_link_libraries(test_path_cover ${gbwtgraph_tests_LIBS})

if (BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

install(TARGETS gfa2gbwt DESTINATION bin)
