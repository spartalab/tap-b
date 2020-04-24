# TAP-B

Fast, efficient C implementation of Algorithm-B for the traffic assignment problem (Robert Dial). Implementation by Dr. Steve Boyles. Parallelization done by Rishabh Thakkar and Karthik Velayutham. This project is part of UT Austin's SPARTA lab. 

As of 12/26/19, this iteration now includes "batching" for very large OD-matrices.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. There are two main options: normal and parallel.

## Configurations

Specific configurations such as convergence gap and maximum number of iterations are things that must be modified in "main.c". 

### Installation/Usage

First make sure to clone the repo as such:

```
git clone https://github.com/spartalab/tap-b.git
```

Then make sure to switch over to the "parallel-pool" branch.

To build an executable of TAP-B, simply run the following:

```
make
```

To build an executable of TAP-B that is parallelized run the following:

```
make parallel
```

To run the normal executable run the following:

```
./bin/tap  net/<network>_net.txt  net/<network>_trips.txt
```

To run the parallelized binary run the following:

```
./bin/tap  net/<network>_net.txt  net/<network>_trips.txt  <num_of_threads>
```

NOTE: if you do not state the number of threads, you will be assigned to a default number of threads based on the number of cores on your machine. If you don't like that, make sure to input the number of threads so that you machine does not get overworked!



### NCTCOG Specific Usage

To make the NCTCOG Specific Executable use for the serial and parallel versions respectively: 
```
make nctcog
```
```
make nctcog-par
```
NCTCOG's network file format can be interpreted as given with all of the details about the conic delay function.

NCTCOG's trip file argument can be a path to a csv file with trip details in the following format:
`<origin_node>,<destination_node>,<class1_demand>,<class2_demand>,...<class10_demand>`

NCTCOG's trip file argument can also be `STREAM` which implies that the OD matrix will be piped in as binary values
in the same order of information as described above in the CSV format.

If the trip file argument is an empty string (`""`) then the OD matrix will be read from the 
stored binary matrices.

To run the normal executable run the following:

```
./bin/tap  <path_to_NCTCOG_net>  <path_to_NCTCOG_trips|STREAM|"">
```

To run the parallelized binary run the following:

```
./bin/tap  <path_to_NCTCOG_net>  <path_to_NCTCOG_trips|STREAM|"">  <num_of_threads>
```

NOTE: if you do not state the number of threads, you will be assigned to a default number of threads based on the number of cores on your machine. If you don't like that, make sure to input the number of threads so that you machine does not get overworked!

## Authors

* **Steve Boyles** - *Primary implementation* 
* **Rishabh Thakkar** - *Parallelized implementation* 
* **Karthik Velayutham** - *Parallelized implementation* 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Keshav Pingali's CS377P: Programming for Performance class with lectures and study about parallel programming
