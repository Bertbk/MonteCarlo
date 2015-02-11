cmake_minimum_required(VERSION 2.8.12)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/ )

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()


# ==============
project(MonteCarlo CXX)
# ==============

option(HAVE_MPI "Compile with MPI" OFF)
option(HAVE_OMP "Compile with OpenMP" OFF)


if(HAVE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_FOUND)
  set(CMAKE_CXX_COMPILER mpicxx)
   add_definitions(-DHAVE_MPI)
  endif(MPI_FOUND)
endif(HAVE_MPI)


add_subdirectory(Src)


