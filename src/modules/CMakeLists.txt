# explicitly name all files in modules/ directory
set(MODULE_FILES 
    kinds.f90
    constants.f90
    system.f90
    random.f90
    util.f90
    selection.f90
   )

# compile module source code files to lib_modules
add_library(lib_modules OBJECT ${MODULE_FILES})
