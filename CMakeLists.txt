cmake_minimum_required(VERSION 3.22)
project(mis C)

set(CMAKE_C_STANDARD 99)


add_executable(mis bspmattoics.c)

target_link_libraries(mis m)

