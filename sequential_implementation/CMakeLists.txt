cmake_minimum_required(VERSION 3.19)
project(sequential_lbm CXX)
#set(CMAKE_PREFIX_PATH "~/Documents/spack/opt/spack/linux-ubuntu22.04-zen2/gcc-11.4.0") # change this to the installation path of HPX
#find_package(HPX REQUIRED)
file(GLOB source_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.cpp")
file(GLOB header_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
add_executable(sequential_lbm  ${source_files} ${header_files})
#target_link_libraries(my_hpx_program HPX::hpx HPX::wrap_main HPX::iostreams_component)