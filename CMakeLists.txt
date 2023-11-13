cmake_minimum_required(VERSION 3.8)
project(explorerPZ LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Finding the pz package if not in neopz buildtree
if(NOT CMAKE_IS_PZ_BUILDTREE)
    find_package(NeoPZ REQUIRED)
endif()

add_executable(sparsemat main.cpp)

target_link_libraries(sparsemat NeoPZ::pz)

