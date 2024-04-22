cmake_minimum_required(VERSION 3.19)
project(parallel_lbm CXX)
find_package(HPX REQUIRED)

# IMPORTANT:
# 
# Do not compile sequential and parallel algorithms together!
# Notice that they typically require different specifications anyways.
# If the sequential algorithm is dealt with first and HPX is not linked against, necessary packages will NOT be found!

set(HEADER_FILES ### General
                 include/access.hpp
                 include/boundaries.hpp 
                 include/collision.hpp
                 include/defines.hpp
                 include/file_interaction.hpp
                 include/macroscopic.hpp
                 include/utils.hpp
                 ### Sequential implementations
                 include/simulation.hpp
                 include/sequential_swap.hpp
                 include/sequential_shift.hpp
                 include/sequential_two_lattice.hpp
                 include/sequential_two_step.hpp
                 ### Parallel implementations
                 include/parallel_framework.hpp
                #  include/parallel_two_lattice.hpp
                 include/parallel_two_lattice_framework.hpp
                 include/parallel_two_step_framework.hpp
                 include/parallel_swap_framework.hpp
                 include/parallel_shift_framework.hpp
                 )

set(SOURCE_FILES ###General
                 src/access.cpp
                 src/boundaries.cpp 
                 src/collision.cpp
                 src/defines.cpp
                 src/file_interaction.cpp
                 src/macroscopic.cpp
                 ### Sequential implementations
                 src/simulation.cpp
                 src/sequential_shift.cpp
                 src/sequential_swap.cpp
                 src/sequential_two_lattice.cpp
                 src/sequential_two_step.cpp
                 ### Parallel implementations
                 src/parallel_framework.cpp
                #  src/parallel_two_lattice.cpp
                 src/parallel_two_lattice_framework.cpp
                 src/parallel_two_step_framework.cpp
                 src/parallel_swap_framework.cpp
                 src/parallel_shift_framework.cpp
                 )

### Two-lattice algorithms

# add_executable(sequential_two_lattice sequential_two_lattice_main.cpp ${SOURCE_FILES})

# add_executable(parallel_two_lattice parallel_two_lattice_main.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_two_lattice HPX::hpx HPX::wrap_main)

add_executable(parallel_two_lattice_framework parallel_two_lattice_framework_main.cpp ${SOURCE_FILES})
target_link_libraries(parallel_two_lattice_framework HPX::hpx HPX::wrap_main)

### Two-step algorithms

# add_executable(sequential_two_step sequential_two_step_main.cpp ${SOURCE_FILES})

# add_executable(parallel_two_step parallel_two_step_main.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_two_step HPX::hpx HPX::wrap_main)

### Swap algorithms

# add_executable(sequential_swap sequential_swap_main.cpp ${SOURCE_FILES})

add_executable(parallel_swap parallel_swap_main.cpp ${SOURCE_FILES})
target_link_libraries(parallel_swap HPX::hpx HPX::wrap_main)

### Shift algorithms

# add_executable(sequential_shift sequential_shift_main.cpp ${SOURCE_FILES})

# add_executable(parallel_shift parallel_shift_main.cpp ${SOURCE_FILES})
# target_link_libraries(parallel_shift HPX::hpx HPX::wrap_main)

### Benchmarks etc
# add_executable(test main_global.cpp)
#target_link_libraries(test HPX::hpx HPX::wrap_main)

add_executable(benchmark main_benchmark.cpp ${SOURCE_FILES})
target_link_libraries(benchmark HPX::hpx HPX::wrap_main)