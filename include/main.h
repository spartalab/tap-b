/*
 * main.h -- Header file for main program.
 */

#ifndef MAIN_H
#define MAIN_H

#include "fileio.h"
#include "bush.h"
#include "networks.h"
#include "tap.h"
#include "datastructures.h"
#include "utils.h"
#include "convexcombination.h"

/* Comment out this line to disable debug logs */
#define DEBUG_MODE

void inferNCTCOGNetwork(network_type *network);
void setBatches(network_type *network, int batchSize, bool warmStart);

void main_TNTP(int argc, char* argv[]);
void main_NCTCOG(int argc, char* argv[]);
void main_NCTCOGFW(int argc, char* argv[]);
void main_NCTCOGFW_conic_calculator(int argc, char* argv[]);
void main_NCTCOGFW_conic_derivative(int argc, char* argv[]);
void main_FWtest(int argc, char* argv[]);

#endif
