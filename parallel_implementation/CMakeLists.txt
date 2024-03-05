cmake_minimum_required(VERSION 3.19)
project(parallel_lbm CXX)
find_package(HPX REQUIRED)

# Careful: Never include both experimental and regular files! This causes multiple definitions.

set(HEADER_FILES include/access.hpp
                 include/boundaries.hpp 
                 include/collision.hpp
                 include/defines.hpp
                 include/macroscopic.hpp
                 #include/shift_sequential_new.hpp
                 include/simulation.hpp
                 #include/swap_sequential.hpp
                 include/two_lattice_parallel.hpp
                 #include/two_step_sequential.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES src/access.cpp
                 src/boundaries.cpp 
                 src/collision.cpp
                 src/defines.cpp
                 src/macroscopic.cpp
                 #src/shift_sequential_new.cpp
                 src/simulation.cpp
                 #src/swap_sequential.cpp
                 src/two_lattice_parallel.cpp
                 #src/two_step_sequential.cpp
                 )
                 
add_executable(parallel_lbm main_two_lattice_parallel.cpp ${SOURCE_FILES})
target_link_libraries(parallel_lbm HPX::hpx HPX::wrap_main HPX::iostreams_component)