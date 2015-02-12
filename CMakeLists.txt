cmake_minimum_required(VERSION 2.8.12)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/ )

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()


# ==============
project(MonteCarlo CXX)
# ==============

option(ENABLE_MPI "Compile with MPI" OFF)
option(ENABLE_OMP "Compile with OpenMP" OFF)


if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_FOUND)
	set(CMAKE_CXX_COMPILER mpicxx)
	add_definitions(-DHAVE_MPI)
  endif(MPI_FOUND)
endif(ENABLE_MPI)

if(ENABLE_OMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		add_definitions(-DHAVE_OMP)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(OPENMP_FOUND)
endif(ENABLE_OMP)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(Src)


