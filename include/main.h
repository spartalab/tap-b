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

void setBatches(network_type *network, int batchSize, bool warmStart);

#ifdef PARALLELISM
int getNumCores();
#endif

#endif
