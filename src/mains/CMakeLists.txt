# compile main executables
add_executable(main main.f90)
target_link_libraries(main lib_modules)

add_executable(ex_particle_filter ex_particle_filter.f90)
target_link_libraries(ex_particle_filter lib_modules)

add_executable(ex_sr_ukf ex_sr_ukf.f90)
target_link_libraries(ex_sr_ukf lib_modules)

add_executable(ex_eoir_uav ex_eoir_uav.f90)
target_link_libraries(ex_eoir_uav lib_modules)
