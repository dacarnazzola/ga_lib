# explicitly name all files in modules/ directory
set(MODULE_FILES 
    kinds.f90
    constants.f90
    system.f90
    timing.f90
    random.f90
    statistics.f90
    sorting.f90
    vector_math.f90
    array_utils.f90
    math_utils.f90
    matrix_math.f90
   )

# compile module source code files to lib_modules
add_library(lib_modules STATIC ${MODULE_FILES})
