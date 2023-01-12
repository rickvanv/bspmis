#TODO: make path variables for MulticoreBSP and BSPedupack folders
cmake_minimum_required(VERSION 3.22)
project(bspmis C)

set(CMAKE_C_STANDARD 99)
set(SOURCE_FILES
        /home/rick/Documents/BSPedupack2.0/bspedupack.c
        /home/rick/Documents/BSPedupack2.0/bspsparse_input.c
        )
#set(HEADER_FILES /home/rick/Documents/MulticoreBSP-for-C/bsp.h /home/rick/Documents/BSPedupack2.0/bspedupack.h /home/rick/Documents/BSPedupack2.0/bspedupack.h /home/rick/Documents/BSPedupack2.0/bspsparse_input.h)

add_library(LIBS ${SOURCE_FILES})
target_include_directories(LIBS PUBLIC /home/rick/Documents/BSPedupack2.0)
target_include_directories(LIBS PUBLIC /home/rick/Documents/MulticoreBSP-for-C)
target_include_directories(LIBS PUBLIC /home/rick/Documents/MulticoreBSP-for-C/lib)

add_executable(bspmis bspmis.c)
target_link_libraries(bspmis LIBS)
target_link_libraries(bspmis /home/rick/Documents/MulticoreBSP-for-C/lib/libmcbsp2.0.4.a)
