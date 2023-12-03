cmake_minimum_required(VERSION 3.19)
project(my_hpx_project CXX)
set(CMAKE_PREFIX_PATH "~/Documents/spack/opt/spack/linux-ubuntu22.04-zen2/gcc-11.4.0") # change this to the installation path of HPX
find_package(HPX REQUIRED)
add_executable(my_hpx_program helloWorld.cpp)
target_link_libraries(my_hpx_program HPX::hpx HPX::wrap_main HPX::iostreams_component)