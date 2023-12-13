cmake_minimum_required(VERSION 3.19)
project(sequential_lbm CXX)
#target_include_directories(sequential_lbm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
#set(HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#file(GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
set(HEADER_FILES #include/access.hpp
                 include/boundaries.hpp 
                 include/defines.hpp
                 #macroscopic.hpp
                 include/simulation_in_cool.hpp
                 #update.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES #src/access.cpp
                 #boundaries.cpp 
                 src/defines.cpp
                 #macroscopic.cpp
                 #simulation.cpp
                 src/simulation_in_cool.cpp
                 #update.cpp
                 #src/utils.cpp
                 )
add_executable(sequential_lbm main.cpp ${HEADER_FILES} ${SOURCE_FILES})
