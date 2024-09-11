# TAP-B

Fast, efficient C implementation of algorithms for the traffic assignment problem.
Default solution method is Algorithm B, developed by Robert Dial.
Link-based solution methods (MSA, Frank-Wolfe, conjugate and biconjugate Frank-Wolfe) are also provided.
Main implementations are by Steve Boyles.
Parallelization capabilities were added by Rishabh Thakkar and Karthik Velayutham.
This project is part of UT Austin's SPARTA (Simulation, Pricing, Adaptive Routing, and Traffic Assignment) lab.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
The main option is to compile the code either for parallelization (default) or serial implementation.
Primary implementation and testing has been done in a Linux ecosystem.
Windows and Mac testing is ongoing; if you encounter difficulties in compiling, try to make the serial build.
If you want to use the parallel version in Windows, the easiest method is to use WSL (the Windows Subsystem for Linux).

## Configurations

Input parameters are specified through a text user interface (TUI).
See `net/SiouxFalls_parameters.txt` for a simple example of such a file; at a minimum it must specify the network and trip table files, and at least one convergence criterion (small gap, max time, or max iterations).
To see all possible parameter settings, see `net/params_example.txt`

### Installation/Usage

First make sure to clone the repo as such:

```
git clone https://github.com/spartalab/tap-b.git
```

To build an executable of TAP-B, simply run the following:

```
make
```

This creates a "standard" build with parallelization, and without compiler optimization.
If you have difficulties with the parallel version, run

```
make serial
```

to turn off parallelization options.

If you want additional performance, run

```
make release
```

to turn on additional compiler optimizations.

If you want to perform additional development, run

```
make debug
```

which turns off optimization but adds debugging symbols.

