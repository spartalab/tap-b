#ifndef CONVEXCOMBINATIONS_H
#define CONVEXCOMBINATIONS_H

#include <limits.h>
#include <stdio.h>
#include <time.h>
#include "fileio.h"
#include "networks.h"
#include "tap.h"
#include "utils.h"

typedef enum {
    MSA,
    FRANK_WOLFE,
    CONJUGATE_FRANK_WOLFE,
    BICONJUGATE_FRANK_WOLFE
} CCalgorithm_type;

typedef enum {
    HEAP_DIJKSTRA,
    PAPE,
    PAPE_WS
} shortestPath_type;

/* Uses a different set of gap functions... with convex combination algorithms
 * it's easy to precompute SPTT, so no need to do so separately in another
 * routine.  */
typedef struct CCparameters {
    //searchDirection_type searchDirection;
    //lineSearch_type      lineSearch;
    shortestPath_type    SP_algo;
    //gap_type             gapFunction;
    double convergenceGap;
    double demandMultiplier;
    double maxTime;
    int    maxIterations;
    int    maxLineSearchIterations;
    bool   calculateBeckmann;
    double CFWmaxweight;
    bool   warmStart;
    char   flowsFile[STRING_SIZE];
    void   (*searchDirection)(network_type*, double**, double**, double**,
                              double, double, double*, struct CCparameters*);
    double (*lineSearch)(network_type*, double**, int, struct CCparameters*);
    double (*gapFunction)(network_type*, double, double);
#if PARALLELISM
    int numThreads;
#endif
} CCparameters_type;

CCparameters_type initializeCCparameters(CCalgorithm_type algo);
void convexCombinations(network_type *network, CCparameters_type *parameters);

void initializeFlows(network_type *network, double **direction,
                     CCparameters_type *parameters);
void shiftFlows(network_type *network, double **direction, double stepSize);

double MSALineSearch(network_type *network, double **targetFlows, int iteration,
                     CCparameters_type *parameters);
double bisection(network_type *network, double **targetFlows, int iteration,
                 CCparameters_type *parameters);
double NewtonSearch(network_type *network, double **targetFlows, int iteration,
                    CCparameters_type *parameters);
double NewtonIteration(network_type *network, double **direction, double lmin,
                       double lmax);

void AONdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters);
void CFWdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters);
void BFWdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters);

double allOrNothing(network_type *network, double **flows, int originZone,
                    int class, CCparameters_type *parameters);
double allOrNothing_par(network_type *network, double *flows, int originZone,
                        int class, CCparameters_type *parameters);
void topoOrderTree(network_type *network, int *SPOrder, arc_type **SPTree); 
void updateLinks(network_type *network, bool rectifyLinks);

double CCrelativeGap(network_type *network, double tstt, double sptt);
double CCaverageExcessCost(network_type *network, double tstt, double sptt);
double evaluateLinkCost(arc_type *arc, double flow);

double calculateConjugacy(network_type *network, double **direction1,
                          double **direction2);
void checkConjugacy(int minVerbosity, network_type *network, double tolerance,
                    double **direction1, double **direction2);
#endif
