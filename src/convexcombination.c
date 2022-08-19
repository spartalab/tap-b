/* This file contains code for link-based algorithms for traffic assignment
 * problem: MSA, Frank-Wolfe, conjugate Frank-Wolfe, and biconjugate
 * Frank-Wolfe.
 *
 * These four routines all have a common structure; all that is different is
 * the target-seeking and line-search behavior.
 */

#include "convexcombination.h"
#include <time.h>

#if PARALLELISM
#include "thpool.h"
#include <pthread.h>
#endif

#if PARALLELISM
//Struct for thread arguments
struct thread_args {
    int id;
    int clss;
    double sptt;
    double **targetFlows;
//    double *targetFlows;
    network_type *network;
    CCparameters_type *parameters;
};

void allOrNothingPool(void* pVoid) {
    struct thread_args *args = (struct thread_args *) pVoid;
    int r = args->id;
    double **targetFlows = args->targetFlows;
//    double *targetFlows = args->targetFlows;
    network_type *network = args->network;
    CCparameters_type *parameters = args->parameters;
    int c = args->clss;
    args->sptt = allOrNothing(network, targetFlows, r, c, parameters);
//    args->sptt = allOrNothing_par(network, targetFlows, r, c, parameters);
}
threadpool thpool;
#endif
/* Initializes convex combination algorithm parameters.  Argument is a
 * CCalgorithm_type enum that sets default search direction and line search
 * functions.  Other parameters are set to reasonable defaults.
 *
 * You must set at least one of the convergence criteria (maxTime,
 * maxIterations, convergenceGap) before calling convexCombinations or it
 * will run forever. 
 */
CCparameters_type initializeCCparameters(CCalgorithm_type algo) {
    CCparameters_type parameters;

    /* Algorithmic choices: search direction, line search, shortest path */
    switch (algo) {
    case MSA:
        parameters.searchDirection = &AONdirection;
        parameters.lineSearch = &MSALineSearch;
        break;
    case FRANK_WOLFE:
        parameters.searchDirection = &AONdirection;
        parameters.lineSearch = &NewtonSearch;
        break;
    case CONJUGATE_FRANK_WOLFE:
        parameters.searchDirection = &CFWdirection;
        parameters.lineSearch = &NewtonSearch;
        break;
    case BICONJUGATE_FRANK_WOLFE:
        parameters.searchDirection = &BFWdirection;
        parameters.lineSearch = &NewtonSearch;
        break;
    default:
        warning(LOW_NOTIFICATIONS, "Unknown algorithm type for CC "
                                   "parameters; setting to Frank-Wolfe");
        parameters.searchDirection = &AONdirection;
        parameters.lineSearch = &bisection;                      
    }
    parameters.SP_algo = HEAP_DIJKSTRA;
    
    /* Convergence criteria */
    parameters.gapFunction = &CCrelativeGap;
    parameters.maxTime = INFINITY;
    parameters.maxIterations = INT_MAX;
    parameters.convergenceGap = 0;

    /* Miscellanea */
    parameters.demandMultiplier = 1;
    if (parameters.lineSearch == &NewtonSearch) {
        parameters.maxLineSearchIterations = 1;
    } else {
        parameters.maxLineSearchIterations = 10;
    }
    parameters.calculateBeckmann = TRUE;
    parameters.warmStart = FALSE;
    parameters.CFWmaxweight = 1e-8;
    snprintf(parameters.flowsFile, sizeof(parameters.flowsFile), "flows.txt");

    return parameters;
}

void convexCombinations(network_type *network, CCparameters_type *parameters) {
    bool converged = FALSE;
    int iteration = 0;
    int ij, c;
    double tstt, sptt, gap;
    double elapsedTime;
    struct timespec tick, tock;

#if PARALLELISM
    thpool = thpool_init(parameters->numThreads);
#endif
    /* Initialize step sizes so first one is a pure AON direction regardless of
     * search direction choice (two if using biconjugate) */
    double stepSize = 0, oldStepSize = 1, oldOldStepSize = 1;
    /* Set up matrices (#arcs x [#classes+1]) for search directions.
     * FW/MSA only use the first of these; CFW the first two; BFW all three.
     * The last column (#classes + 1) is used to store the aggregated search
     * direction across all classes; this saves effort in computing shifts */
    declareMatrix(double, direction, network->numArcs, network->numClasses + 1);
    declareMatrix(double, oldDirection, network->numArcs,
                  network->numClasses + 1);
    declareMatrix(double, oldOldDirection, network->numArcs,
                  network->numClasses + 1);
    double **temp; /* used for swapping */

    elapsedTime = 0;

    /* If warmStart is true, assume that the flows in network are feasible */
    if (parameters->warmStart == FALSE) {
        initializeFlows(network, direction, parameters);
    }

    while (converged == FALSE) {
        /* Find search direction with whatever algorithm and parameters are
         * relevant */
        clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

        parameters->searchDirection(network, direction, oldDirection,
                                    oldOldDirection, oldStepSize,
                                    oldOldStepSize, &sptt, parameters);
        
        /* This gives us the information needed to report gap and check
         * convergence */
        tstt = TSTT(network);
        gap = parameters->gapFunction(network, tstt, sptt);
        clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
        elapsedTime += (double)((1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec)) * 1.0/1000000000; /* Exclude gap calculations from run time */

        displayMessage(LOW_NOTIFICATIONS, "Iteration %d: gap %.15f, obj %.15g,"
                      " time %.3f\n",iteration, gap, BeckmannFunction(network),
                      elapsedTime);
        if (elapsedTime > parameters->maxTime) converged = TRUE;
        if (iteration >= parameters->maxIterations) converged = TRUE;
        if (gap < parameters->convergenceGap) converged = TRUE;

        /* Now find step size and shift flows (the line search function
           also does the shifting) */
        stepSize = parameters->lineSearch(network, direction, iteration,
                                          parameters);

        
        /* Now set up for next iteration by scaling old directions,
         * shifting down old step sizes, and permuting array flows (no need to
         * copy, just swap pointers to alias them) */
        for (ij = 0; ij < network->numArcs; ij++) {
            /* Inner loop goes all the way up to numClasses to get the
               aggregate column right too */
            for (c = 0; c <= network->numClasses; c++) {
                direction[ij][c] *= (1 - stepSize);
                oldDirection[ij][c] *= (1 - stepSize);
            }
        }
        
        oldOldStepSize = oldStepSize; 
        oldStepSize = stepSize;

        temp = oldOldDirection;
        oldOldDirection = oldDirection;
        oldDirection = direction;
        direction = temp; /* This will be overwritten next time around */
        
        iteration++;
    }
    
    deleteMatrix(direction, network->numArcs);
    deleteMatrix(oldDirection, network->numArcs);
    deleteMatrix(oldOldDirection, network->numArcs);
}

/* Generate an initial solution if none is provided; do an all-or-nothing
 * assignment based on free-flow costs. */
void initializeFlows(network_type *network, double **direction,
                     CCparameters_type *parameters) {
    int ij, c;
    double stepSize, sptt;

    /* Initialize class flows to zero, then get total flows and costs */
    for (ij = 0; ij < network->numArcs; ij++) {
        for (c = 0; c < network->numClasses; c++) {
            network->arcs[ij].classFlow[c] = 0;
        }
    }
    updateLinks(network, TRUE);

    /* Calculate an all-or-nothing assignment based on free-flow times and
     * shift all flow in that direction.  Passing a number of parameters as
     * zero/null because a lot of features in AONdirection are not needed for
     * this. */
    AONdirection(network, direction, NULL, NULL, 0, 0, &sptt, parameters);
    stepSize = 1;
    shiftFlows(network, direction, stepSize);
}

/* Shift flows in a given direction with the given step size.  Then update
 * costs and everything else */
void shiftFlows(network_type *network, double **direction, double stepSize) {
    int ij, c;

    for (ij = 0; ij < network->numArcs; ij++) {
        for (c = 0; c < network->numClasses; c++) {
            network->arcs[ij].classFlow[c] += stepSize * direction[ij][c];
        }
        network->arcs[ij].flow += stepSize * direction[ij][network->numClasses];
    }

    updateLinks(network, FALSE); /* No need to rectify, the same convex combo
                                    applies to class and total flows */
}

/* Ensure link flows are sum of class flows; update costs and derivatives */
void updateLinks(network_type *network, bool rectifyFlows) {
    int ij, c;

    if (rectifyFlows == TRUE) {
        for (ij = 0; ij < network->numArcs; ij++) {
            network->arcs[ij].flow = 0;
            for (c = 0; c < network->numClasses; c++) {
                network->arcs[ij].flow += network->arcs[ij].classFlow[c];
            }
        }
    }

    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].cost =
                network->arcs[ij].calculateCost(&network->arcs[ij]);
        network->arcs[ij].der =
                network->arcs[ij].calculateDer(&network->arcs[ij]);
    }
    
}

/* targetFlows and parameters are not used in this function */
double MSALineSearch(network_type *network, double **direction, int iteration,
                     CCparameters_type *parameters) {
    /* Add 2 to iteration so the step size sequence is
     *    1/2, 1/3, 1/4, ....
     * (notice that the first iteration is 'iteration 0' */
    return 1.0 / (iteration + 2);
    /* Suppress compiler warnings for unused variables */
    displayMessage(FULL_DEBUG, "%p %p %p", network, direction, parameters);
}

/* the iteration argument is not used in this function */
double bisection(network_type *network, double **direction, int iteration,
                 CCparameters_type *parameters) {
    int ij, c, k; /* ij indexes arcs; c classes; k is an iteration counter */
    /* lambda is step size; lmin and lmax are lower and upper bounds; der
     * is the derivative test at the midpoint */
    double lambda = 1, lmin = 0, lmax = 1, der; 
    double tempFlow;

    for (k = 0; k < parameters->maxLineSearchIterations; k++) {
        lambda = (lmax + lmin) / 2;
        der = 0;
        for (c = 0; c < network->numClasses; c++) {
            changeFixedCosts(network, c);
            for (ij = 0; ij < network->numArcs; ij++) {
                /* Get total flow based on complete shift */
                tempFlow = network->arcs[ij].flow 
                            + lambda * direction[ij][network->numClasses];
                der += evaluateLinkCost(&network->arcs[ij], tempFlow)
                        * direction[ij][c];
            }
        }
        if (der > 0) lmax = lambda; else lmin = lambda;
    }

    shiftFlows(network, direction, lambda);
    return lambda;
    
    /* Suppress compiler warnings for unused variables */
    displayMessage(FULL_DEBUG, "%d", iteration);
}

/* Does a Newton line search over (possibly multiple) iterations; calls
 * NewtonIteration successively, updating link flows each time */
double NewtonSearch(network_type *network, double **direction, int iteration,
                    CCparameters_type *parameters) {
    int k;
    double stepSize, totalStepSize = 0;
    double lmin = 0, lmax = 1; /* Initial bounds on Newton */

    /* Do a single step regardless... */
    stepSize = NewtonIteration(network, direction, lmin, lmax);
    totalStepSize = stepSize;

    /* ...and if we need to do more, shift flows and keep going (note k
     * starts at 1 to reflect the iteration already done */
    for (k = 1; k < parameters->maxLineSearchIterations; k++) {
        shiftFlows(network, direction, stepSize);
        /* Adjust valid Newton bounds based on last step size */
        lmin -= stepSize;
        lmax -= stepSize;
        stepSize = NewtonIteration(network, direction, lmin, lmax);
        totalStepSize += stepSize;
    }

    shiftFlows(network, direction, stepSize); /* one final shift */
    return totalStepSize; /* return a final step size for a last shift in main
                        convex combinations routine */

    /* Suppress compiler warnings for unused variables */
    displayMessage(FULL_DEBUG, "%d", iteration);
}

/* Does a *single* Newton iteration starting from the current flows */
double NewtonIteration(network_type *network, double **direction, double lmin,
                       double lmax) {
    int ij, c;
    double lambda = 0, numer = 0, denom = 0;
    for (c = 0; c < network->numClasses; c++) {
        changeFixedCosts(network, c);
        for (ij = 0; ij < network->numArcs; ij++) {
            numer += network->arcs[ij].calculateCost(&network->arcs[ij]) 
                    * direction[ij][c];
        }
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        denom += network->arcs[ij].calculateDer(&network->arcs[ij])
                * direction[ij][network->numClasses]
                * direction[ij][network->numClasses];
    }                    
    
    if (denom != 0) { /* Typical case */
        lambda = - numer / denom;
        lambda = min(lambda, lmax);
        lambda = max(lambda, lmin);
    } else { /* Denominator is zero; indicates constant link performance
                functions, so shift everything to shortest paths */
         lambda = 1;
    }
    return lambda;
}


/* This is the "critical loop", calculating an all-or-nothing direction from
 * the current flows.  It calls allOrNothing, which in turn calls a shortest
 * path routine.
 *
 * This is the main opportunity for parallelizing... the allOrNothing calls
 * (the loop over each origin and class) can all be done independently of each
 * other.  They all write to the targetFlows matrix (incrementing its values
 * with the flows from that origin) but do not read from it.  Hopefully this
 * makes it easy to parallelize, we can consider other structures if this will
 * be a problem.
 */
void AONdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters) {
    int ij, r, c; /* Indices: ij for arcs, r for origin, c for class */
    double **targetFlows = direction; /* alias to save memory */

    /* Initialize the flow shifts, including "last column" of totals */
    for (ij = 0; ij < network->numArcs; ij++) {
        for (c = 0; c < network->numClasses; c++) {
            direction[ij][c] = 0;
        }
        direction[ij][network->numClasses] = 0;
    }

    *sptt = 0; /* calculate incrementally */

#if PARALLELISM
    struct thread_args args[network->numZones];
    for (r = 0; r < network->numZones; ++r) {
        args[r].id = r;
        args[r].network = network;
//        args[r].targetFlows = calloc(network->numArcs, sizeof(double));
        args[r].targetFlows = targetFlows;
        args[r].clss = -1;
        args[r].parameters = parameters;
        args[r].sptt = -1;
    }
#endif

    for (c = 0; c < network->numClasses; c++) {
        changeFixedCosts(network, c);
#if PARALLELISM
        for (r = 0; r < network->numZones; ++r) {
            args[r].clss = c;
            thpool_add_work(thpool, (void (*)(void *)) allOrNothingPool, (void*)&args[r]);
        }
        thpool_wait(thpool);
        for (r = 0; r < network->numZones; ++r) {
            if (args[r].sptt < 0) {
                fatalError("SPTT is negative for origin %d is %f", r, args[r].sptt);
            }
//            for (int i = 0; i < network->numArcs; ++i) {
//                targetFlows[i][c] += args[r].targetFlows[i];
//                args[r].targetFlows[i] = 0.0;
//            }
            *sptt += args[r].sptt;
        }
#else
        for (r = 0; r < network->numZones; r++) {
            *sptt += allOrNothing(network, targetFlows, r, c, parameters);
        }
#endif
    }
//#if PARALLELISM
//    for (r = 0; r < network->numZones; ++r) {
//        deleteVector(args[r].targetFlows);
//    }
//#endif

    for (ij = 0; ij < network->numArcs; ij++) {
        for (c = 0; c < network->numClasses; c++) {
            direction[ij][c] = targetFlows[ij][c]
                                - network->arcs[ij].classFlow[c];
            direction[ij][network->numClasses] += direction[ij][c];
        }
    }

    return;
    /* Suppress compiler warnings */
    displayMessage(FULL_DEBUG, "%p %p %f %f", oldDirection, oldOldDirection,
                   oldStepSize, oldOldStepSize);
}

void CFWdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters) {
    int ij, c;
    double **AON = direction; /* alias to save memory */
    double numer, denom, alpha;

    /* First find all-or-nothing direction, ignoring irrelevant params  */
    AONdirection(network, AON, NULL, NULL, 0, 0, sptt, parameters);
    
    if (oldStepSize == 1) return; /* Do pure AON step if last step was full */
    /* ^-- this code is OK since AONdirection and direction point to the same
     * values */

    /* Now calculate parameters for conjugacy */
    numer = 0; denom = 0;
    for (c = 0; c < network->numClasses; c++) {
        changeFixedCosts(network, c);
        for (ij = 0; ij < network->numArcs; ij++) {
            numer += network->arcs[ij].der
                        * AON[ij][c] * oldDirection[ij][c];
            denom += network->arcs[ij].der 
                        * oldDirection[ij][c]
                        * (AON[ij][c] - oldDirection[ij][c]);
        }
    }
    if (denom != 0) {
        alpha = min(numer / denom, 1 - parameters->CFWmaxweight);
    } else {
        alpha = 0;
    }

    /* Initialize the "last column" of total flow shifts */
    for (ij = 0; ij < network->numArcs; ij++) {
        direction[ij][network->numClasses] = 0;
    }
    for (c = 0; c < network->numClasses; c++) {
        for (ij = 0; ij < network->numArcs; ij++) {
            direction[ij][c] = alpha * oldDirection[ij][c]
                            + (1 - alpha) * AON[ij][c];
            direction[ij][network->numClasses] += direction[ij][c];
        }
    }

    //displayMessage(FULL_NOTIFICATIONS, "Conjugacy level: %f\n",     
    //               calculateConjugacy(network, direction, oldDirection));

    return;
    /* Suppress compiler warnings */
    displayMessage(FULL_DEBUG, "%p %f", oldOldDirection, oldOldStepSize);
}

void BFWdirection(network_type *network, double **direction, 
                           double **oldDirection, double **oldOldDirection,
                           double oldStepSize, double oldOldStepSize,
                           double *sptt,  CCparameters_type *parameters) {
    int ij, c;
    double **AON= direction; /* alias to save memory */
    double numer, denom, mu, nu, beta0, beta1, beta2; 

    /* First find all-or-nothing direction, ignoring irrelevant params  */
    AONdirection(network, AON, NULL, NULL, 0, 0, sptt, parameters);

    /* Do pure AON step if either of last two steps were full */
    if (oldStepSize == 1 || oldOldStepSize == 1) return;
        
    /* Now calculate parameters for conjugacy */
    numer = 0; denom = 0;
    for (c = 0; c < network->numClasses; c++) {
        changeFixedCosts(network, c);
        for (ij = 0; ij < network->numArcs; ij++) {
            numer += network->arcs[ij].der
                        * AON[ij][c] * oldOldDirection[ij][c];
            denom += network->arcs[ij].der 
                        * oldOldDirection[ij][c]
                        * (oldOldDirection[ij][c] - oldDirection[ij][c]);
        }
    }
    if (denom == 0) return; /* Do pure AON step in this exceptional case */
    mu = - numer / denom * (1 - oldStepSize);
    
    numer = 0; denom = 0;
    for (c = 0; c < network->numClasses; c++) {
        changeFixedCosts(network, c);
        for (ij = 0; ij < network->numArcs; ij++) {
            numer += network->arcs[ij].der
                        * AON[ij][c] * oldDirection[ij][c];
            denom += network->arcs[ij].der 
                        * oldDirection[ij][c] * oldDirection[ij][c];
        }
    }
    if (denom == 0) return; /* Do pure AON step in this exceptional case */
    nu = mu * oldStepSize / (1 - oldStepSize) - numer / denom;
    
    beta0 = 1 / (1 + mu + nu);
    beta1 = nu * beta0;
    beta2 = mu * beta0;
    
    /* Initialize the "last column" of total flow shifts */
    for (ij = 0; ij < network->numArcs; ij++) {
        direction[ij][network->numClasses] = 0;
    }

    for (c = 0; c < network->numClasses; c++) {
        for (ij = 0; ij < network->numArcs; ij++) {
            direction[ij][c] = beta0 * AON[ij][c]
                            + (beta1 + beta2) * oldDirection[ij][c]
                            + beta2 * beta2 / (1 - oldStepSize)
                             * (oldOldDirection[ij][c] - oldDirection[ij][c]);
            direction[ij][network->numClasses] += direction[ij][c];
        }
    }

    //displayMessage(FULL_NOTIFICATIONS, "Conjugacy level 1: %f\n",     
    //               calculateConjugacy(network, direction, oldDirection));

    //displayMessage(FULL_NOTIFICATIONS, "Conjugacy level 2: %f\n",     
    //               calculateConjugacy(network, direction, oldOldDirection));

    
}

/* Finds an all-or-nothing assignment from a given origin.  Involves three main
 * steps: shortest-path finding; identifying a topological order on the tree;
 * and performing a descending pass to calculate flows.  
 *
 * In the process, also calculates the SPTT for this origin (it comes almost
 * for free from the shortest path algo) and returns this for use in gap
 * calculations 
 */
double allOrNothing(network_type *network, double **flows, int originZone,
                    int class, CCparameters_type *parameters) {
    int curnode, backnode, backarc, i;
    double *remainingVehicles;
    double originSPTT = 0;
    declareVector(arc_type *, SPTree, network->numArcs);
    declareVector(int, SPOrder, network->numNodes);
    declareVector(double, SPLabels, network->numNodes);
    int origin = nodeclass2origin(network, originZone, class);

    /* Find all-to-one shortest paths from origin */
    switch (parameters->SP_algo) {
    case HEAP_DIJKSTRA:
        heapDijkstraNew(originZone, SPLabels, SPTree, network);    
        break;
    case PAPE:
        BellmanFordNew(originZone, SPLabels, SPTree, NULL, network, DEQUE);    
        break;
    case PAPE_WS:
        fatalError("Warm-started PAPE not available in this implementation of "
                   "convex combinations (storing trees takes too much space).");
        break;
    default:
        fatalError("Unknown shortest path algorithm %d\n", parameters->SP_algo);    
    }
  
    /* Calculate shortest path time (for gap) */
    for (i = 0; i < network->numZones; i++) {
        if (SPTree[i] != NULL) { /* Ordinary case, node is reachable */
            originSPTT += SPLabels[i] * network->demand[origin][i];
        } else { /* No path found... only an issue if there is demand */
            if (network->demand[origin][i] > 0 && i != originZone) {
                fatalError("No path found from %d to %d but demand exists!",
                            originZone, i);
            }
        } 
    }
    
    topoOrderTree(network, SPOrder, SPTree);

    /* Now load vehicles onto this tree in reverse topological order -- only one
       sweep of the tree is needed to find all flow from this origin.
       remainingVehicles gives the number of vehicles at each node which have
       not yet been fully assigned to their path.
       
       In the TNTP file format, origins are numbered first.  The second loop
       thus picks up where the first one left off.
      */
    remainingVehicles = SPLabels; /* Re-use array to save memory, we don't need
                                     the shortest path labels anymore */

    for (i = 0; i < network->numZones; i++)
        remainingVehicles[i] = network->demand[origin][i];
    for (; i < network->numNodes; i++)
        remainingVehicles[i] = 0;

    /* Here is the main loop, in reverse topological order */
    for (i = network->numNodes - 1; i > 0; i--) {
        curnode = SPOrder[i];
        if (curnode == originZone) break;
        if (SPTree[curnode] != NULL) { /* Usual case, can push vehicles back */
            backnode = SPTree[curnode]->tail;
            backarc = ptr2arc(network, SPTree[curnode]);
#ifdef PARALLELISM
            pthread_mutex_lock(&network->arc_muts[backarc]);
#endif
            flows[backarc][class] += remainingVehicles[curnode];
#ifdef PARALLELISM
            pthread_mutex_unlock(&network->arc_muts[backarc]);
#endif
            remainingVehicles[backnode] += remainingVehicles[curnode];
        } else { /* No path found... only an issue if there is demand */
            if (remainingVehicles[curnode] > 0) {
                fatalError("allOrNothing: no path from %d to %d despite "
                           " demand %f existing there!", origin, curnode,
                           remainingVehicles[curnode]);
            }
        } 
        remainingVehicles[curnode] = 0;
    }
    
    deleteVector(SPOrder);
    deleteVector(SPTree);
    deleteVector(SPLabels);
    return originSPTT;
}
/* Finds an all-or-nothing assignment from a given origin.  Involves three main
 * steps: shortest-path finding; identifying a topological order on the tree;
 * and performing a descending pass to calculate flows.
 *
 * In the process, also calculates the SPTT for this origin (it comes almost
 * for free from the shortest path algo) and returns this for use in gap
 * calculations
 */
double allOrNothing_par(network_type *network, double *flows, int originZone,
                    int class, CCparameters_type *parameters) {
    int curnode, backnode, backarc, i;
    double *remainingVehicles;
    double originSPTT = 0;
    declareVector(arc_type *, SPTree, network->numArcs);
    declareVector(int, SPOrder, network->numNodes);
    declareVector(double, SPLabels, network->numNodes);
    int origin = nodeclass2origin(network, originZone, class);

    /* Find all-to-one shortest paths from origin */
    switch (parameters->SP_algo) {
    case HEAP_DIJKSTRA:
        heapDijkstraNew(originZone, SPLabels, SPTree, network);
        break;
    case PAPE:
        BellmanFordNew(originZone, SPLabels, SPTree, NULL, network, DEQUE);
        break;
    case PAPE_WS:
        fatalError("Warm-started PAPE not available in this implementation of "
                   "convex combinations (storing trees takes too much space).");
        break;
    default:
        fatalError("Unknown shortest path algorithm %d\n", parameters->SP_algo);
    }

    /* Calculate shortest path time (for gap) */
    for (i = 0; i < network->numZones; i++) {
        if (SPTree[i] != NULL) { /* Ordinary case, node is reachable */
            originSPTT += SPLabels[i] * network->demand[origin][i];
        } else { /* No path found... only an issue if there is demand */
            if (network->demand[origin][i] > 0 && i != originZone) {
                fatalError("No path found from %d to %d but demand exists!",
                            originZone, i);
            }
        }
    }

    topoOrderTree(network, SPOrder, SPTree);

    /* Now load vehicles onto this tree in reverse topological order -- only one
       sweep of the tree is needed to find all flow from this origin.
       remainingVehicles gives the number of vehicles at each node which have
       not yet been fully assigned to their path.

       In the TNTP file format, origins are numbered first.  The second loop
       thus picks up where the first one left off.
      */
    remainingVehicles = SPLabels; /* Re-use array to save memory, we don't need
                                     the shortest path labels anymore */

    for (i = 0; i < network->numZones; i++)
        remainingVehicles[i] = network->demand[origin][i];
    for (; i < network->numNodes; i++)
        remainingVehicles[i] = 0;

    /* Here is the main loop, in reverse topological order */
    for (i = network->numNodes - 1; i > 0; i--) {
        curnode = SPOrder[i];
        if (curnode == originZone) break;
        if (SPTree[curnode] != NULL) { /* Usual case, can push vehicles back */
            backnode = SPTree[curnode]->tail;
            backarc = ptr2arc(network, SPTree[curnode]);
            flows[backarc] += remainingVehicles[curnode];
            remainingVehicles[backnode] += remainingVehicles[curnode];
        } else { /* No path found... only an issue if there is demand */
            if (remainingVehicles[curnode] > 0) {
                fatalError("allOrNothing: no path from %d to %d despite "
                           " demand %f existing there!", origin, curnode,
                           remainingVehicles[curnode]);
            }
        }
        remainingVehicles[curnode] = 0;
    }

    deleteVector(SPOrder);
    deleteVector(SPTree);
    deleteVector(SPLabels);
    return originSPTT;
}


/* A specialized version of topologicalOrder (in networks.c) that does not rely
   on a full network struct, but only a list of predecessors defining a tree.
   */
void topoOrderTree(network_type *network, int *SPOrder, arc_type **SPTree) {
    int i, curnode;
    int curOrder;
    int pathPos;
    declareVector(bool, labeled, network->numNodes);
    declareVector(int, pathNodes, network->numNodes);
    
    for (i = 0; i < network->numNodes; i++) {
        labeled[i] = FALSE;
    }
    
    curOrder = 0;
    for (i = 0; i < network->numNodes; i++) {
        /* Trace a path back from i to a labeled node (or the origin) */
        pathPos = 0;
        curnode = i;
        while (labeled[curnode] == FALSE) {
            pathNodes[pathPos++] = curnode;
            labeled[curnode] = TRUE;
            if (SPTree[curnode] != NULL) curnode = SPTree[curnode]->tail;
        }
        /* Now add those nodes to the topological order */
        while (--pathPos >= 0) {
            SPOrder[curOrder++] = pathNodes[pathPos];
        }
    }
    
    if (curOrder != network->numNodes)
        fatalError("refreshTopologicalOrder did not label all nodes!");
    
    deleteVector(labeled);
    deleteVector(pathNodes);
}

double CCrelativeGap(network_type *network, double tstt, double sptt) {
    return (tstt - sptt) / sptt;

    /* Suppress compiler warnings */
    displayMessage(FULL_DEBUG, "%p", network);
}

double CCaverageExcessCost(network_type *network, double tstt, double sptt) {
    return (tstt - sptt) / network->totalODFlow;
}

/* Evaluate a link performance function WITHOUT changing the underlying flow */
double evaluateLinkCost(arc_type *arc, double flow) {
    double oldFlow = arc->flow;
    double cost;
    arc->flow = flow;
    cost = arc->calculateCost(arc);
    arc->flow = oldFlow;
    return cost;
}

/* This function is useful in error-checking... for the given directions,
 * calculates direction1 * H * direction2, where H is the diagonal matrix of
 * link travel time derivatives.  If the two directions are conjugate w.r.t. H,
 * this should be zero.
 *
 * Note that the same derivative applies across all user classes, thanks to the
 * assumptions of constant tolls and the use of a tollFactor (1/VOT) rather
 * than VOT itself, in the definition of link cost.
 *
 * */
double calculateConjugacy(network_type *network, double **direction1,
                          double **direction2) {
    int ij, c;
    double s = 0;
    for (ij = 0; ij < network->numArcs; ij++) {
        for (c = 0; c < network->numClasses; c++) {
            s += direction1[ij][c] * network->arcs[ij].der * direction2[ij][c];
        }
    }

    return s;
}

void checkConjugacy(int minVerbosity, network_type *network, double tolerance,
                    double **direction1, double **direction2) {
    double conj;

    if (verbosity < minVerbosity) return;

    conj = fabs(calculateConjugacy(network, direction1, direction2));
    if (conj > tolerance) {
        displayMessage(minVerbosity, "Conjugacy check failed with value %f\n",
                       conj);
    }
}

