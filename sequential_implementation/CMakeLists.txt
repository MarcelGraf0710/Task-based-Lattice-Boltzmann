cmake_minimum_required(VERSION 3.19)
project(sequential_lbm CXX)

# Careful: Never include both experimental and regular files! This causes multiple definitions.

set(HEADER_FILES include/access.hpp
                 include/boundaries.hpp 
                 include/collision.hpp
                 include/defines.hpp
                 include/macroscopic.hpp
                 #include/shift_sequential.hpp
                 include/simulation.hpp
                 include/swap_sequential_experimental.hpp
                 include/two_lattice_sequential.hpp
                 include/two_step_sequential.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES src/access.cpp
                 src/boundaries.cpp 
                 src/collision.cpp
                 src/defines.cpp
                 src/macroscopic.cpp
                 #src/shift_sequential.cpp
                 src/simulation.cpp
                 src/swap_sequential_experimental.cpp
                 src/two_lattice_sequential.cpp
                 src/two_step_sequential.cpp
                 )
                 
add_executable(sequential_lbm main_two_step.cpp ${SOURCE_FILES})
