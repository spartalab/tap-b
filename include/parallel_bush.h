#ifndef PARALLEL_BUSH_H
#define PARALLEL_BUSH_H

#include <limits.h>
#include <math.h>
#include "tap.h"
#include "networks.h"
#include "datastructures.h"
#include "utils.h"
#include "bush.h"

/* Parallized main Algorithm B helper functions */
void updateBushB_par(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters);
void updateFlowsB_par(int origin, network_type *network, bushes_type *bushes,
                      algorithmBParameters_type *parameters);


/* Parallelized counter part of bush manipulation */
void scanBushes_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, bool longestUsed);
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
void exactCostUpdate_par(int ij, double shift, network_type *network);
void checkFlows_par(network_type *network, bushes_type *bushes, int t_id);


#endif