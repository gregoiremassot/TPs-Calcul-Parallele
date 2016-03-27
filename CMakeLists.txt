cmake_minimum_required (VERSION 2.6)
project (code_DF)

# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
set (EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")

include_directories ("${PROJECT_SOURCE_DIR}/geometries")
add_subdirectory (geometries)
set (EXTRA_LIBS ${EXTRA_LIBS} geometries)

include_directories ("${PROJECT_SOURCE_DIR}/output")
add_subdirectory (output)
set (EXTRA_LIBS ${EXTRA_LIBS} output)

include_directories ("${PROJECT_SOURCE_DIR}/solver")
add_subdirectory(solver)
set (EXTRA_LIBS ${EXTRA_LIBS} solver)

# add the executable
add_executable(code_DF main.cc)
target_link_libraries (code_DF  ${EXTRA_LIBS})
