cmake_minimum_required(VERSION 3.21)
project(Neural_Network)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS  "-Wall -Werror -pedantic --coverage")
add_subdirectory(Google_tests)
add_executable(Neural_Network main.cpp)