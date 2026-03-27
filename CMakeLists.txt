# minimum version based on Ubuntu 22.04 on 2022-09-16
cmake_minimum_required(VERSION 3.24)

# set the project name
project(dbase VERSION 0.0 LANGUAGES Fortran)

# executables will be built and placed in bin/ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# modules will compile and be stored in build/modfiles directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modfiles)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# source code will be in src/ folder
add_subdirectory(src)
