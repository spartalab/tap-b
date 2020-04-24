/*
 * parallel_bush.h -- This is the header file for the bush and merge data 
 * structures that are used specifically for the parallelization workloads.
 * 
 * All the method implementation details are the exact same as those in bush.h,
 * the only difference is that this file defines functions that implement separate
 * data structures for each thread, which are then merged for the final result 
 * used to solve equilibrium for the network.
 * 
 */


#ifndef PARALLEL_BUSH_H
#define PARALLEL_BUSH_H

#include <limits.h>
#include <math.h>
#include "tap.h"
#include "networks.h"
#include "datastructures.h"
#include "utils.h"
#include "bush.h"
pthread_mutex_t shift_lock;

/* Parallized main Algorithm B helper functions */
void updateBushB_par(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters);
bool updateFlowsB_par(int origin, network_type *network, bushes_type *bushes,
                      algorithmBParameters_type *parameters);


/* Parallelized counter part of bush manipulation */
void scanBushes_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, scan_type LPrule);
void reconstructMerges_par(int origin, network_type *network, bushes_type *bushes);
bool rescanAndCheck_par(int origin, network_type *network, bushes_type *bushes,
                        algorithmBParameters_type *parameters);
void updateFlowPass_par(int origin, network_type *network, bushes_type *bushes,
                        algorithmBParameters_type *parameters);
void calculateBushFlows_par(int origin,network_type *network,bushes_type *bushes);
void pushBackFlowSimple_par(int j, int origin, network_type *network,
                            bushes_type *bushes);
void pushBackFlowMerge_par(merge_type *merge, network_type *network,
                           bushes_type *bushes, int t_id);
void rectifyMerge_par(int j, merge_type *merge, bushes_type *bushes, int t_id);
void newtonFlowShift_par(int j, merge_type *merge, int origin,
                         network_type *network, bushes_type *bushes,
                         algorithmBParameters_type *parameters);

/* Parallelized counterpart of utility functions */
void exactCostUpdate_par(int ij, double shift, network_type *network, int class);
void classUpdate_par(int hi, int class, double shift,  network_type *network);
void checkFlows_par(network_type *network, bushes_type *bushes, int t_id);


#endif
