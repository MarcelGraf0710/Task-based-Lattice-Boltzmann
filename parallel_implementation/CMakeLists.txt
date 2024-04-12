cmake_minimum_required(VERSION 3.19)
project(parallel_lbm CXX)
find_package(HPX REQUIRED)

# Careful: Never include both experimental and regular files! This causes multiple definitions.

set(HEADER_FILES include/access.hpp
                 include/boundaries.hpp 
                 include/collision.hpp
                 include/defines.hpp
                 include/macroscopic.hpp
                 include/shift_sequential.hpp
                 include/simulation.hpp
                 #include/swap_sequential.hpp
                 include/two_lattice_sequential.hpp
                 include/two_lattice_parallel.hpp
                 include/parallel_framework.hpp
                 include/parallel_two_lattice_framework.hpp
                 include/parallel_two_step_framework.hpp
                 #include/parallel_swap_framework.hpp
                 include/parallel_shift_framework.hpp
                 include/two_step_sequential.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES src/access.cpp
                 src/boundaries.cpp 
                 src/collision.cpp
                 src/defines.cpp
                 src/macroscopic.cpp
                 src/shift_sequential.cpp
                 src/simulation.cpp
                 #src/swap_sequential.cpp
                 src/two_lattice_parallel.cpp
                 src/two_step_sequential.cpp
                 src/two_lattice_parallel.cpp
                 src/two_lattice_sequential.cpp
                 src/parallel_framework.cpp
                 src/parallel_two_lattice_framework.cpp
                 src/parallel_two_step_framework.cpp
                 # src/parallel_swap_framework.cpp
                 src/parallel_shift_framework.cpp
                 )

### Two-lattice algorithms

# add_executable(sequential_two_lattice main_two_lattice.cpp ${SOURCE_FILES})
# target_link_libraries(sequential_two_lattice HPX::hpx HPX::wrap_main HPX::iostreams_component)

# add_executable(parallel_two_lattice main_two_lattice_framework.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_two_lattice HPX::hpx HPX::wrap_main)

### Two-step algorithms

add_executable(sequential_two_step main_two_step.cpp ${SOURCE_FILES})
target_link_libraries(sequential_two_step HPX::hpx HPX::wrap_main HPX::iostreams_component)

# add_executable(parallel_two_step main_two_step_framework.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_two_step HPX::hpx HPX::wrap_main HPX::iostreams_component)

### Swap algorithms

# add_executable(test_swap main_swap_framework.cpp ${SOURCE_FILES})
# target_link_libraries(test_swap HPX::hpx HPX::wrap_main HPX::iostreams_component)

### Shift algorithms

# add_executable(sequential_shift main_shift.cpp ${SOURCE_FILES})
# target_link_libraries(sequential_shift HPX::hpx HPX::wrap_main HPX::iostreams_component)

# add_executable(parallel_shift main_shift_framework.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_shift HPX::hpx HPX::wrap_main HPX::iostreams_component)
