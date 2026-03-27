# compile test executables
add_executable(test_kinds test_kinds.f90)
target_link_libraries(test_kinds lib_modules)

add_executable(test_debug test_debug.f90)
target_link_libraries(test_debug lib_modules)

add_executable(test_random test_random.f90)
target_link_libraries(test_random lib_modules)

add_executable(test_statistics test_statistics.f90)
target_link_libraries(test_statistics lib_modules)

add_executable(test_sorting test_sorting.f90)
target_link_libraries(test_sorting lib_modules)

add_executable(test_vector_math test_vector_math.f90)
target_link_libraries(test_vector_math lib_modules)
