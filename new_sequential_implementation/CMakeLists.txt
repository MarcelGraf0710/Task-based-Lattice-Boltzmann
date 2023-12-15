cmake_minimum_required(VERSION 3.19)
project(sequential_lbm CXX)
#target_include_directories(sequential_lbm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
#set(HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#file(GLOB HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "*.hpp")
#set(SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
set(HEADER_FILES include/access.hpp
                 include/boundaries.hpp 
                 include/defines.hpp
                 include/macroscopic.hpp
                 include/update.hpp
                 include/simulation.hpp
                 include/utils.hpp
                 )

set(SOURCE_FILES src/access.cpp
                 src/boundaries.cpp 
                 src/defines.cpp
                 src/macroscopic.cpp
                 src/simulation.cpp
                 src/update.cpp
                 )
add_executable(sequential_lbm main.cpp ${HEADER_FILES} ${SOURCE_FILES})
