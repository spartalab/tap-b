cmake_minimum_required(VERSION 3.16)
project(tap_b C CXX)
find_package(Threads REQUIRED)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
set(CMAKE_CXX_FLAGS -pthread -Wall -Wextra -Wwrite-strings -Wno-parentheses -Winline -Wpedantic -Warray-bounds)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(include)
add_definitions(-DPARALLELISM -DNCTCOG_ENABLED)
add_executable(tap_nctcog_par
        src/bush.c
        src/bush.c_minmem
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
        src/bush.c_minmem
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
add_executable(tap_par
        src/bush.c
        src/bush.c_minmem
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
        src/bush.c_minmem
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
target_compile_definitions(tap_nctcog_par PUBLIC -DPARALLELISM -DNCTCOG_ENABLED)
target_compile_definitions(tap_nctcog PUBLIC -DNCTCOG_ENABLED)
target_compile_definitions(tap_par PUBLIC -DPARALLELISM)

target_link_libraries(tap_nctcog_par ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap_nctcog ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap_par ${CMAKE_THREAD_LIBS_INIT} m)
target_link_libraries(tap ${CMAKE_THREAD_LIBS_INIT} m)
