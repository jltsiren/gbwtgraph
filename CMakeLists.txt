# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(seqwish)
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
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
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
ExternalProject_Add(mmmultimap
  GIT_REPOSITORY "https://github.com/ekg/mmmultimap.git"
  GIT_TAG "09b6618837deb526de169e71256a0cb5cdf06dc2"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(mmmultimap SOURCE_DIR)
set(mmmultimap_INCLUDE "${SOURCE_DIR}/src")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  GIT_REPOSITORY "https://github.com/SaschaWitt/ips4o.git"
  GIT_TAG "bff3ccf0bf349497f2bb10f825d160b792236367"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  GIT_REPOSITORY "https://github.com/vgteam/BBHash.git"
  GIT_TAG "ccc22805eae48aca2a01a7ff5ef8e002857020d3"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# GBWT
ExternalProject_Add(gbwt
  GIT_REPOSITORY "https://github.com/jltsiren/gbwt.git"
  GIT_TAG "debb9f58c1fa4360aa11984011781b1b8d1e7063"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(gbwt SOURCE_DIR)
set(gbwt_INCLUDE "${SOURCE_DIR}/include")
set(gbwt_LIB "${SOURCE_DIR}/lib")

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  GIT_REPOSITORY "https://github.com/vgteam/libhandlegraph.git"
  GIT_TAG "962018b84d8525ca93e79efd5f3726a0683bcfdf"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/src/handlegraph/src/include")
set(handlegraph_LIB "${INSTALL_DIR}/src/handlegraph-build")

set(CMAKE_BUILD_TYPE Debug)

add_library(gbwtgraph STATIC
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(gbwtgraph sdsl-lite)
add_dependencies(gbwtgraph mmmultimap)
add_dependencies(gbwtgraph ips4o)
add_dependencies(gbwtgraph bbhash)
add_dependencies(gbwtgraph gbwt)
add_dependencies(gbwtgraph handlegraph)
target_include_directories(gbwtgraph PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}")
#target_link_libraries(gbwtgraph
#  "${handlegraph_LIB}/libhandlegraph.a"
#  "${sdsl-lite_LIB}/libsdsl.a"
#  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
#  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
#  "-latomic"
#  z)

# set up our target executable and specify its dependencies and includes
add_executable(gfa2gbwt
  ${CMAKE_SOURCE_DIR}/gfa2gbwt.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(gfa2gbwt gbwtgraph)
target_include_directories(gfa2gbwt PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}")
target_link_libraries(gfa2gbwt
  "${gbwt_LIB}/libgbwt.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "-latomic"
  z)

add_executable(test_gbwtgraph
  ${CMAKE_SOURCE_DIR}/tests/test_gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_gbwtgraph gbwtgraph)
target_include_directories(test_gbwtgraph PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}"
  "${googletest_INCLUDE}")
target_link_libraries(test_gbwtgraph
  "${gbwt_LIB}/libgbwt.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${googletest_LIB}/libgtest.a"
  "${googletest_LIB}/libgtest_main.a"
  "-latomic"
  z)

add_executable(test_minimizer
  ${CMAKE_SOURCE_DIR}/tests/test_minimizer.cpp
  ${CMAKE_SOURCE_DIR}/gbwtgraph.cpp
  ${CMAKE_SOURCE_DIR}/gfa.cpp
  ${CMAKE_SOURCE_DIR}/minimizer.cpp
  ${CMAKE_SOURCE_DIR}/utils.cpp)
add_dependencies(test_minimizer gbwtgraph)
target_include_directories(test_minimizer PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${bsort_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gbwt_INCLUDE}"
  "${handlegraph_INCLUDE}"
  "${googletest_INCLUDE}")
target_link_libraries(test_minimizer
  "${gbwt_LIB}/libgbwt.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${googletest_LIB}/libgtest.a"
  "${googletest_LIB}/libgtest_main.a"
  "-latomic"
  z)

if (BUILD_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
