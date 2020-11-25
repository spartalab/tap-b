cmake_minimum_required(VERSION 3.13)
project(tap_b C CXX)
set(CMAKE_CXX_STANDARD 14)
include_directories(include)
if (WIN32)
    set(CMAKE_CXX_FLAGS -Wall -Wextra -Wwrite-strings -Wno-parentheses -Winline -Wpedantic -Warray-bounds)
else()
    find_package(Threads REQUIRED)
    set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
    set(THREADS_PREFER_PTHREAD_FLAG TRUE)
    set(CMAKE_CXX_FLAGS -pthread -Wall -Wextra -Wwrite-strings -Wno-parentheses -Winline -Wpedantic -Warray-bounds)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(tap_nctcog_par
        src/bush.c
        src/cpthread.c
        src/convexcombination.c
        src/datastructures.c
        src/fileio.c
        src/main.c
        src/networks.c
        src/parallel_bush.c
        src/tap.c
        src/thpool.c
        src/utils.c
        flows.txt)
add_executable(tap_nctcog
        src/bush.c
        src/convexcombination.c
        src/datastructures.c
        src/fileio.c
        src/main.c
        src/networks.c
        src/parallel_bush.c
        src/tap.c
        src/utils.c
        flows.txt)
add_executable(tap_par
        src/bush.c
        src/cpthread.c
        src/convexcombination.c
        src/datastructures.c
        src/fileio.c
        src/main.c
        src/networks.c
        src/parallel_bush.c
        src/tap.c
        src/thpool.c
        src/utils.c
        flows.txt)
add_executable(tap
        src/bush.c
        src/convexcombination.c
        src/datastructures.c
        src/fileio.c
        src/main.c
        src/networks.c
        src/parallel_bush.c
        src/tap.c
        src/utils.c
        flows.txt)
target_compile_definitions(tap_nctcog_par PRIVATE -DPARALLELISM -DNCTCOG_ENABLED)
target_compile_definitions(tap_nctcog PRIVATE -DNCTCOG_ENABLED)
target_compile_definitions(tap_par PRIVATE -DPARALLELISM)

if(!WIN32)
target_link_libraries(tap_nctcog_par ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap_nctcog ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap_par ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap ${CMAKE_THREAD_LIBS_INIT} m)
endif()