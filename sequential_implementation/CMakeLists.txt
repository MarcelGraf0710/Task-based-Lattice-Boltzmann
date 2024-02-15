cmake_minimum_required(VERSION 3.19)
project(sequential_lbm CXX)
#target_include_directories(sequential_lbm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
#set(HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#file(GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
set(HEADER_FILES include/access.hpp
                 include/boundaries.hpp 
                 include/collision.hpp
                 include/defines.hpp
                 include/macroscopic.hpp
                 include/shift_sequential.hpp
                 include/simulation.hpp
                 include/swap_sequential.hpp
                 include/two_lattice_sequential.hpp
                 include/two_step_sequential.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES src/access.cpp
                 src/boundaries.cpp 
                 src/collision.cpp
                 src/defines.cpp
                 src/macroscopic.cpp
                 src/simulation.cpp
                 src/two_lattice_sequential.cpp
                 src/two_step_sequential.cpp
                 
                 )
                 
#add_executable(sequential_lbm main.cpp ${SOURCE_FILES} ${HEADER_FILES} )
add_executable(sequential_lbm main_two_step.cpp ${SOURCE_FILES})
