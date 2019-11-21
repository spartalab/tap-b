cmake_minimum_required(VERSION 3.13)
project(FinalProject)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(FinalProject
    include/bush.h
    include/datastructures.h
    include/fileio.h
    include/main.h
    include/parallel_bush.h
    include/networks.h
    include/tap.h
    include/utils.h
    include/thpool.h
    net/Anaheim_net.txt
    net/Anaheim_trips.txt
    net/Austin_sdb_net.txt
    net/Austin_sdb_trips.txt
    net/Barcelona_net.txt
    net/Barcelona_trips.txt
    net/Braess_net.txt
    net/Braess_trips.txt
    net/ChicagoRegional_net.txt
    net/ChicagoRegional_trips.txt
    net/ChicagoSketch_net.txt
    net/ChicagoSketch_trips.txt
    net/SiouxFalls_net.txt
    net/SiouxFalls_trips.txt
    src/bush.c
    src/datastructures.c
    src/fileio.c
    src/main.c
    src/networks.c
    src/parallel_bush.c
    src/tap.c
    src/utils.c
    src/thpool.c
    Makefile
    README.md )
