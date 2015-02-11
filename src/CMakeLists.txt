include_directories(${CMAKE_SOURCE_DIR}/Src)


add_executable(MonteCarlo
Message.cpp
Mesh.cpp
Point.cpp
Database.cpp
Main_Monte.cpp)


if(HAVE_MPI)
target_link_libraries(MonteCarlo ${MPI_CXX_LIBRARIES})
endif(HAVE_MPI)
