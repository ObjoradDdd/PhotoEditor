cmake_minimum_required(VERSION 3.10)
project(fe)

set(CMAKE_CXX_STANDARD 20)

# Add AVX and AVX2 compiler flags
if(MSVC)
    add_compile_options(/arch:AVX2)
else()
    add_compile_options(-mavx -mavx2)
endif()
add_executable(fe main.cpp)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES main.cpp)
add_executable(openmp_test_clion ${SOURCE_FILES})

