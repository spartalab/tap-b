/*
 * tap.h -- Header file for general-purpose traffic assignment functions --
 * calculating the Beckmann function, gap measures, TSTT, and so on.  Also
 * contains special routines for quickly calculating special cases of BPR
 * functions.
 */
#ifndef TAP_H
#define TAP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "networks.h"
#include "utils.h"
#include "datastructures.h"

typedef enum {
    BPR,
    BPR_DER,
    BPR_INT,
} cost_type;

typedef enum {
    RELATIVE_GAP_1, /* Definition of relative gap from class */
    RELATIVE_GAP_2, /* An alternative using lower Beckmann bounds */
    AEC,            /* Average excess cost */
    AEC_OB,         /* Special implementation of AEC for origin-based algos */
    MEC,            /* Maximum excess cost */
    PASS            /* Dummy -- skip gap calculation */
} gap_type;


//////////////////////////////////
// General network calculations //
//////////////////////////////////

double BeckmannFunction(network_type *network);
double calculateGap(network_type *network, gap_type gapFunction);
void finalizeNetwork(network_type *network);
double linkCost(arc_type *arc, cost_type costFunction);
double SPTT(network_type *network);
double TSTT(network_type *network);
void updateAllCosts(network_type *network);
void updateAllCostDers(network_type *network);

double generalBPRcost(struct arc_type *arc);
double generalBPRder(struct arc_type *arc);
double generalBPRint(struct arc_type *arc, bool includeFixedCost);
double linearBPRcost(struct arc_type *arc);
double linearBPRder(struct arc_type *arc);
double linearBPRint(struct arc_type *arc, bool includeFixedCost);
double quarticBPRcost(struct arc_type *arc);
double quarticBPRder(struct arc_type *arc);
double quarticBPRint(struct arc_type *arc, bool includeFixedCost);
double conicCost(struct arc_type *arc);
double conicDer(struct arc_type *arc);
double conicInt(struct arc_type *arc, bool includeFixedCost);


int arcNumber(network_type *network, arc_type *arc);
double averageExcessCost(network_type *network);
double relativeGap1(network_type *network);
double relativeGap2(network_type *network);

double classCost(network_type *network, int class, double timeFactor,
                 double tollFactor, double distanceFactor);
double classRevenue(network_type *network, int class);
double classDistance(network_type *network, int class);
double classTravelTime(network_type *network, int class);
double classGeneralizedCost(network_type *network, int class);

void makeStronglyConnectedNetwork(network_type *network);

#endif
