cmake_minimum_required(VERSION 3.22)

# set the project name
project(edmd C)

# Compiler flags
SET(GCC_COMPILE_FLAGS "-Wall -Ofast -funroll-loops -march=native")
# SET(GCC_COMPILE_FLAGS "-Wall -g")
string(APPEND CMAKE_C_FLAGSE "${GCC_COMPILE_FLAGS}")
# Library directory
add_subdirectory(PseudoRNG)

# add the executable
add_executable(edmd mdSingle.c)

# For PRNG
target_link_libraries(edmd PUBLIC PseudoRNG)
target_include_directories(edmd PUBLIC 
"${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/PseudoRNG"
)
# For math functions
target_link_libraries(edmd PUBLIC m)