cmake_minimum_required(VERSION 2.8.12)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/ )

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()


# ==============
project(MonteCarlo CXX)
# ==============

option(WITH_MPI "Compile with MPI" OFF)
option(WITH_OMP "Compile with OpenMP" OFF)

add_executable(MonteCarlo
Message.cpp
Point.cpp
Database.cpp
Main_Monte.cpp)
