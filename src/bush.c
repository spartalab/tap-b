#include "bush.h"
#include "thpool.h"
#include <pthread.h> /*used in other parts of the assignment */
#define NUM_THREADS 2
#define PAR 1
//Struct for thread arguments
struct thread_args {
    int id;
    bushes_type *bushes;
    network_type *network;
    algorithmBParameters_type *parameters;
};


void updateBushPool(void* pVoid) {
    struct thread_args *args = (struct thread_args *) pVoid;
    int id = args->id;
    bushes_type *bushes = args->bushes;
    network_type *network = args->network;
    algorithmBParameters_type *parameters = args->parameters;
    updateBushB_par(id, network, bushes, parameters, id);
}

void updateFlowsPool(void* pVoid) {
    struct thread_args *args = (struct thread_args *) pVoid;
    int id = args->id;
    bushes_type *bushes = args->bushes;
    network_type *network = args->network;
    algorithmBParameters_type *parameters = args->parameters;
    updateFlowsB_par(id, network, bushes, parameters, id);
}

void AlgorithmB(network_type *network, algorithmBParameters_type *parameters) {
    /* Strong connectivity check */
    makeStronglyConnectedNetwork(network);

    /* Allocate memory for bushes */
    int origin, i, iteration = 0;
    bushes_type *bushes = createBushes(network);

    double elapsedTime = 0, gap = INFINITY;
    struct thread_args args[network->numZones];
    for (int j = 0; j < network->numZones; ++j) {
        args[j].id = j;
        args[j].network = network;
        args[j].parameters = parameters;
        args[j].bushes = bushes;
    }
    threadpool thpool = thpool_init(NUM_THREADS);
    /* Initialize */
    clock_t stopTime = clock(); /* used for timing */
    initializeBushesB(network, bushes, parameters);
    /* Iterate */
    do {
        iteration++;
        for (int j = 0; j < network->numZones; ++j) {
            thpool_add_work(thpool, (void*) updateBushPool, (void*)&args[j]);
        }
        thpool_wait(thpool);
        for (int j = 0; j < network->numZones; ++j) {
            thpool_add_work(thpool, (void*) updateFlowsPool, (void*)&args[j]);
        }
        thpool_wait(thpool);

        /* Shift flows -- Inner iterations*/
        for (i = 0; i < parameters->innerIterations; i++) {
            for (int j = 0; j < network->numZones; ++j) {
                thpool_add_work(thpool, (void*) updateFlowsPool, (void*)&args[j]);

            }
            thpool_wait(thpool);
        }

        /* Check gap and report progress */
        elapsedTime += ((double)(clock() - stopTime)) / CLOCKS_PER_SEC; /* Exclude gap calculations from run time */
        gap = calculateGap(network, parameters->gapFunction);
        displayMessage(LOW_NOTIFICATIONS, "Iteration %ld: gap %.15f, Beckmann %.13g, time %.3f s.\n", iteration, gap, BeckmannFunction(network), elapsedTime);
        stopTime = clock();
    } while (elapsedTime < parameters->maxTime && iteration < parameters->maxIterations && gap > parameters->convergenceGap);

    /* Clean up */
    deleteBushes(network, bushes);

}


/*
 * initializeAlgorithmBParameters -- sets default values for algorithm
 * parameter choices.
 */
algorithmBParameters_type initializeAlgorithmBParameters() {
    algorithmBParameters_type parameters;

    parameters.gapFunction = RELATIVE_GAP_1;
    parameters.convergenceGap = 0;
    parameters.maxTime = INFINITY;
    parameters.maxIterations = LONG_MAX;

    parameters.innerIterations = 20;
    parameters.shiftReps = 1;
    parameters.rescanAfterShift = FALSE;

    parameters.thresholdGap = 0;
    parameters.minCostDifference = 0;
    parameters.minLinkFlowShift = 0;

    parameters.minLinkFlow = 1e-14;
    parameters.minDerivative = 1e-6;
    parameters.newtonStep = 1;
    parameters.numNewtonShifts = 1;

    parameters.SPQueueDiscipline = DEQUE;

    parameters.createInitialBush = &initialBushShortestPath;
    parameters.topologicalOrder = &genericTopologicalOrder;
    parameters.linkShiftB = &exactCostUpdate;

    return parameters;
}

/*
 * initialBushShortestPath -- Creates initial bushes based on the shortest path
 * tree at free-flow times.  Arguments are the origin for the bush, pointers to
 * the network and bush data structures (see comments on network.h and bush.h),
 * and the parameters struct.
 */
void initialBushShortestPath(int origin, network_type *network,
                             bushes_type *bushes,
                             algorithmBParameters_type *parameters) {
    /* Find shortest path and set bush preds */
    arcIndexBellmanFord(origin, bushes->SPcost, bushes->pred[origin], network,
                        parameters->SPQueueDiscipline);

    /* Set bush topological order */
    bushes->numMerges[origin] = 0;
    bushes->merges[origin] = newVector(0, merge_type *);
    parameters->topologicalOrder(origin, network, bushes, parameters);
}

/*
 * initialBushBFS -- Create initial bushes based on the breadth-first search
 * tree.  Likely to be much worse than shortest paths (see
 * initialBushShortestPath), but allows initialization in time linear in
 * network size.
 *
 * Not yet implemented.
 */
void initialBushBFS(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    fatalError("initialBushBFS is not yet implemented, call something else.");

    /* Suppress warnings */
   displayMessage(FULL_DEBUG, "%d%p%p%p", origin, network, bushes, parameters);
}

/*
 * genericTopologicalOrder -- Find a topological order using the standard
 * algorithm (finding and marking nodes with no marked predecessors).
 * Arguments are the origin corresponding to the bush, and the network, bush,
 * and parameters data structures.
 */
void genericTopologicalOrder(int origin, network_type *network,
                             bushes_type *bushes,
                             algorithmBParameters_type *parameters) {
    arcListElt *curArc;
    int i, j, m,  next, highestMerge = 0;
    
    declareVector(int, indegree, network->numNodes);
    for (i = 0; i < network->numNodes; i++) {
        indegree[i] = 1; /* By default non-origin nodes are assumed to have 1
                            incoming link; merges and origin handled below */
         bushes->bushOrder[origin][i] = NO_PATH_EXISTS;
    }
    for (i = 0; i < network->numNodes; i++) {
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            indegree[i] = bushes->merges[origin][m]->numApproaches;
        }
    }
    indegree[origin] = 0;

    queue_type LIST = createQueue(network->numNodes, network->numNodes);
    next = 0;
    for (i = 0; i < network->numNodes; i++)
        if (indegree[i] == 0)
            enQueue(&LIST, i);
    while (LIST.curelts > 0) {
        i = deQueue(&LIST);
        bushes->bushOrder[origin][next] = i;
        if (isMergeNode(origin, i, bushes) == TRUE) highestMerge = next;
        next++;
        for (curArc = network->nodes[i].forwardStar.head; curArc != NULL;
                curArc = curArc->next) {
            if (isInBush(origin, ptr2arc(network, curArc->arc), network,
                         bushes) == TRUE) {
                j = curArc->arc->head;
                indegree[j]--;
                if (indegree[j] == 0) enQueue(&LIST, j);
            }
        }
    }
    if (next < network->numNodes) {
        int j = 0;
        j += 1;
        fatalError("Graph given to bushTopologicalOrder contains a cycle.");
    }
    bushes->lastMerge[origin] = highestMerge;

    deleteQueue(&LIST);
    deleteVector(indegree);
    /* Suppress warnings -- parameters is not used in this implementation */
    if (0) displayMessage(FULL_DEBUG, "%p", parameters);
}

/*
 * mergeFirstTopologicalOrder -- a specialized topological order algorithm
 * which aims to place merge nodes lower in the ordering.  This results in a
 * lower lastMerge value for the bush, potentially speeding up later
 * computations.
 *
 * Not yet implemented.
 */
void mergeFirstTopologicalOrder(int origin, network_type *network,
                                bushes_type *bushes,
                                algorithmBParameters_type *parameters) {
    fatalError("mergeFirstTopologicalOrder is not yet implemented.");
    /* Suppress warnings */
    displayMessage(FULL_DEBUG,"%d%p%p%p", origin, network, bushes, parameters);
}


/*
 * The next three functions provide an interface for how merge nodes are
 * identified.  Non-merge nodes have a non-negative pred value in the bushes
 * struct, indicating the ID of the predecessor arc in the bush.  Merge nodes
 * have a negative pred value, relating to the index in the array of merge
 * nodes in the bush struct.  Merges 0, 1, 2, ... are indexed by -1, -2, -3 as
 * predecessors (otherwise 0 is ambiguous).
 */

/* 
 * isMergeNode -- determine whether node i is a merge node in the bush
 * corresponding to origin.
 */
bool isMergeNode(int origin, int i, bushes_type *bushes) {
    return (bushes->pred[origin][i] < 0 && i != origin) ? TRUE : FALSE;
}

/*
 * pred2merge -- For a merge node, convert the pred value to the corresponding
 * index in the merge array.
 */
int pred2merge(int ij) {
    return -(ij + 1);
}

/*
 * merge2pred -- Inverts pred2merge; converts an index in the merge array to
 * the corresponding pred value.
 */
int merge2pred(int m) {
    return -(m + 1);
}

/*
 * isInBush -- determines whether a given link (ij) is in a particular bush
 * (the one corresponding to origin).  With the current bush data structures,
 * this requires examining the predecessor of link (i,j)'s head.
 */
bool isInBush(int origin, int ij, network_type *network, bushes_type *bushes) {
    int backarc, j = network->arcs[ij].head, m;
    merge_type *merge;
    if (isMergeNode(origin, j, bushes) == FALSE) {
        return (bushes->pred[origin][j] == ij) ? TRUE : FALSE;
    } else {
        m = pred2merge(bushes->pred[origin][j]);
        merge = bushes->merges[origin][m];
        for (backarc = 0; backarc < merge->numApproaches; backarc ++) {
            if (merge->approach[backarc] == ij) return TRUE;
        }
    }
    return FALSE;
}


/*
 * createBushes -- Initialize the bushes data structure for a particular
 * network, allocating memory as needed.
 */
bushes_type *createBushes(network_type *network) {
    declareScalar(bushes_type, bushes);
    int i;

    bushes->LPcost = newVector(network->numNodes, double);
    bushes->SPcost = newVector(network->numNodes, double);
    bushes->flow = newVector(network->numArcs, double);
    bushes->nodeFlow = newVector(network->numNodes, double);
    bushes->pred = newMatrix(network->numZones, network->numNodes, int);
    bushes->bushOrder = newMatrix(network->numZones, network->numNodes, int);
    bushes->lastMerge = newVector(network->numZones, int);
    bushes->numMerges = newVector(network->numZones, int);
    bushes->merges = newVector(network->numZones, merge_type**);

#if PAR
    bushes->LPcost_par = newMatrix(network->numNodes, network->numNodes,double);
    bushes->SPcost_par = newMatrix(network->numNodes, network->numNodes,double);
    bushes->flow_par = newMatrix(network->numNodes, network->numArcs,double);
    bushes->nodeFlow_par = newMatrix(network->numNodes, network->numNodes,double);
#endif

    for (i = 0; i < network->numZones; i++) {
       bushes->numMerges[i] = 0;
    }
    return bushes;
}

/*
 * deleteBushes -- Cleans up when done, deallocating memory associated with
 * bushes.
 */
void deleteBushes(network_type *network, bushes_type *bushes) {
    int i, m;

    deleteVector(bushes->LPcost);
    deleteVector(bushes->SPcost);
    deleteVector(bushes->flow);
    deleteVector(bushes->nodeFlow);
    deleteMatrix(bushes->pred, network->numZones);
    deleteMatrix(bushes->bushOrder, network->numZones);
    deleteVector(bushes->lastMerge);
    for (i = 0; i < network->numZones; i++) {
        for (m = 0; m < bushes->numMerges[i]; m++) {
            deleteMerge(bushes->merges[i][m]);
        }
        deleteVector(bushes->merges[i]);
    }
    deleteVector(bushes->numMerges);
    deleteVector(bushes->merges);
//#if PAR
//    deleteMatrix(bushes->LPcost_par, network->numNodes);
//    deleteMatrix(bushes->SPcost_par, network->numNodes);
//    deleteMatrix(bushes->flow_par, network->numArcs);
//    deleteMatrix(bushes->nodeFlow_par, network->numNodes);
//#endif
    deleteScalar(bushes);
}


/*
 * initializeBushesB -- set up initial bushes at the beginning of the algorithm
 * (as distinct from createBushesB, which merely allocates memory).  This
 * involves the following steps:
    1. Zero out flow on all links, set travel times to free-flow
    2. Cycle through bushes, call getInitialBush.
    3. Calculate bush link flows with a downward pass
    4. Increase global link flows based on the bush
    5. Optionally recalculate travel times, depending on parameters -- doing so
       produces a more accurate initial solution, at the cost of some more BPR
       function evaluations.  This is controlled by setting calculateCost to
       point to a different function.
    6. At end, update derivatives.
*/
void initializeBushesB(network_type *network, bushes_type *bushes,
                       algorithmBParameters_type *parameters) {
    int i, ij;

    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].flow = 0;
        network->arcs[ij].cost = network->arcs[ij].freeFlowTime;
    }

    for (i = 0; i < network->numZones; i++) {
        /* createInitialBush also sets preds, bushOrder */
        parameters->createInitialBush(i, network, bushes, parameters); 
        calculateBushFlows(i, network, bushes);
        for (ij = 0; ij < network->numArcs; ij++) {
            network->arcs[ij].flow += bushes->flow[ij];
            network->arcs[ij].cost =
                network->arcs[ij].calculateCost(&network->arcs[ij]);
        }
    }

    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].der =
            network->arcs[ij].calculateDer(&network->arcs[ij]);
    }
}

/*
 * scanBushes -- simultaneously find shortest and longest paths for a given
 * bush.  If the longestUsed argument is TRUE, the longest path search will
 * restrict attention to longest *used* paths.  If FALSE, longest bush paths
 * will be calculated regardless of whether there is flow on them.
*/
void scanBushes(int origin, network_type *network, bushes_type *bushes,
                algorithmBParameters_type *parameters, bool longestUsed) {
    int h, i, hi, m, curnode, curarc;
    double tempcost;
    merge_type *merge;

    for (i = 0; i < network->numNodes; i++) {
        bushes->LPcost[i] = -INFINITY;
        bushes->SPcost[i] = INFINITY;
    }

    /* Ensure costs are up to date */
    if (parameters->linkShiftB != &exactCostUpdate) updateAllCosts(network);

    bushes->SPcost[origin] = 0;
    bushes->LPcost[origin] = 0;
    for (curnode = 1; curnode < network->numNodes; curnode++) {
        i = bushes->bushOrder[origin][curnode];
        /* Iterate over incoming links */
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            merge = bushes->merges[origin][m];
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->SPcost[h] + network->arcs[hi].cost;
                if (tempcost < bushes->SPcost[i]) {
                    bushes->SPcost[i] = tempcost;
                    merge->SPlink = curarc;
                }
                tempcost = bushes->LPcost[h] + network->arcs[hi].cost;      
                if (tempcost > bushes->LPcost[i]
                     && (longestUsed == FALSE ||
                         merge->approachFlow[curarc] > parameters->minLinkFlow)
                   ) {
                    bushes->LPcost[i] = tempcost;
                    merge->LPlink = curarc;
                }
            }
        } else { /* Only one incoming bush link, not much to do */
            hi = bushes->pred[origin][i];
            h = network->arcs[hi].tail;
            bushes->LPcost[i] = bushes->LPcost[h] + network->arcs[hi].cost;
            bushes->SPcost[i] = bushes->SPcost[h] + network->arcs[hi].cost;
        }
   }

}

/*
 * scanBushes_par -- simultaneously find shortest and longest paths for a given
 * bush.  If the longestUsed argument is TRUE, the longest path search will
 * restrict attention to longest *used* paths.  If FALSE, longest bush paths
 * will be calculated regardless of whether there is flow on them.
*/
void scanBushes_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, bool longestUsed, int t_id) {
//    displayMessage(FULL_NOTIFICATIONS, "Top of scan %d\n", t_id);
    int h, i, hi, m, curnode, curarc;
    double tempcost;
    merge_type *merge;

    for (i = 0; i < network->numNodes; i++) {
        bushes->LPcost_par[t_id][i] = -INFINITY;
        bushes->SPcost_par[t_id][i] = INFINITY;
    }

    /* Ensure costs are up to date */
    if (parameters->linkShiftB != &exactCostUpdate) updateAllCosts(network);

    bushes->SPcost_par[t_id][origin] = 0;
    bushes->LPcost_par[t_id][origin] = 0;
    for (curnode = 1; curnode < network->numNodes; curnode++) {
        i = bushes->bushOrder[origin][curnode];
        /* Iterate over incoming links */
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            if(&bushes->merges[origin][m] < (merge_type **) 0x1000000) {
                int a = 0;
                displayMessage(FULL_NOTIFICATIONS, "Bad address %x %d\n", &bushes->merges[origin][m], t_id);
            }
            merge = bushes->merges[origin][m];
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->SPcost_par[t_id][h] + network->arcs[hi].cost;
                if (tempcost < bushes->SPcost_par[t_id][i]) {
                    bushes->SPcost_par[t_id][i] = tempcost;
                    merge->SPlink = curarc;
                }
                tempcost = bushes->LPcost_par[t_id][h] + network->arcs[hi].cost;
                if (tempcost > bushes->LPcost_par[t_id][i]
                    && (longestUsed == FALSE ||
                        merge->approachFlow[curarc] > parameters->minLinkFlow)
                        ) {
                    bushes->LPcost_par[t_id][i] = tempcost;
                    merge->LPlink = curarc;
                }
            }
        } else { /* Only one incoming bush link, not much to do */
            hi = bushes->pred[origin][i];
            h = network->arcs[hi].tail;
            bushes->LPcost_par[t_id][i] = bushes->LPcost_par[t_id][h] + network->arcs[hi].cost;
            bushes->SPcost_par[t_id][i] = bushes->SPcost_par[t_id][h] + network->arcs[hi].cost;
        }
    }
//    displayMessage(FULL_NOTIFICATIONS, "Hello end of scan %d\n", t_id);


}

/*
 * updateBushB -- updates the topology of a bush: deleting unused links, and
 * adding links which can serve as shortcuts.  Follows Nie (2009) for how this
 * is done, first checking a stronger criterion for adding bush links, then as
 * a fallback checking a weaker one.
 *
 * New links are marked using the NEW_LINK constant, for use in
 * reconstructMerges (the function which rebuilds the merge data structures
 * after bush updating.)
 */
void updateBushB(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters) {
    int ij, i, j, newArcs = 0;
   
    /* First update labels... ignoring longest unused paths since those will be
     * removed in the next step. */
    scanBushes(origin, network, bushes, parameters, TRUE);
    calculateBushFlows(origin, network, bushes);
   
    /* Make a first pass... */
    for (ij = 0; ij < network->numArcs; ij++) {
        /* Mark links with near-zero contribution for removal by setting to
         * exactly zero */
        if (bushes->flow[ij] < parameters->minLinkFlow) bushes->flow[ij] = 0;
        if (bushes->flow[ij] > 0) continue; /* Link is already in the bush, no
                                               need to add it again */
        /* See if the link provides a shortcut using the strict criterion */
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        if (bushes->LPcost[i] == -INFINITY && bushes->LPcost[j] > -INFINITY)
            continue; /* No path to extend */
        if (bushes->SPcost[i] + network->arcs[ij].cost < bushes->SPcost[j]
            && bushes->LPcost[i] < bushes->LPcost[j]
            && (network->arcs[ij].tail == origin
                || network->arcs[ij].tail >= network->firstThroughNode)) 
        {
            bushes->flow[ij] = NEW_LINK;
            newArcs++;
        /* Never delete shortest path tree... should be OK with floating point
         * comparison since this is how SPcost is calculated */
        } else if (bushes->SPcost[i]+network->arcs[ij].cost==bushes->SPcost[j] 
                   && bushes->flow[ij] == 0) {
            bushes->flow[ij] = NEW_LINK;
        }
    }
   
    /* If strict criterion fails, try a looser one */
    if (newArcs == 0) {
        for (ij = 0; ij < network->numArcs; ij++) {
            i = network->arcs[ij].tail;
            j = network->arcs[ij].head;
            if (bushes->LPcost[i]==-INFINITY && bushes->LPcost[j]>-INFINITY)
                continue; /* No path to extend */
            if (bushes->flow[ij] == 0 && bushes->LPcost[i] < bushes->LPcost[j]
                && (network->arcs[ij].tail == origin
                    || network->arcs[ij].tail >= network->firstThroughNode))
            {
                bushes->flow[ij] = NEW_LINK;
            }
        }
    }      

   /* Finally update bush data structures: delete/add merges, find a new
    * topological order, rectify approach proportions */
    reconstructMerges(origin, network, bushes);
    parameters->topologicalOrder(origin, network, bushes, parameters);
}

/*
 * updateBushB_par -- updates the topology of a bush: deleting unused links, and
 * adding links which can serve as shortcuts.  Follows Nie (2009) for how this
 * is done, first checking a stronger criterion for adding bush links, then as
 * a fallback checking a weaker one.
 *
 * New links are marked using the NEW_LINK constant, for use in
 * reconstructMerges (the function which rebuilds the merge data structures
 * after bush updating.)
 */
void updateBushB_par(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters, int t_id) {
//    displayMessage(FULL_NOTIFICATIONS, "Top of update bushb %d %d\n", origin ,t_id);

    int ij, i, j, newArcs = 0;

    /* First update labels... ignoring longest unused paths since those will be
     * removed in the next step. */
    scanBushes_par(origin, network, bushes, parameters, TRUE, t_id);
    calculateBushFlows_par(origin, network, bushes, t_id);
    /* Make a first pass... */
    for (ij = 0; ij < network->numArcs; ij++) {
        /* Mark links with near-zero contribution for removal by setting to
         * exactly zero */
        if (bushes->flow_par[t_id][ij] < parameters->minLinkFlow) bushes->flow_par[t_id][ij] = 0;
        if (bushes->flow_par[t_id][ij] > 0) continue; /* Link is already in the bush, no
                                               need to add it again */
        /* See if the link provides a shortcut using the strict criterion */
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        if (bushes->LPcost_par[t_id][i] == -INFINITY && bushes->LPcost_par[t_id][j] > -INFINITY)
            continue; /* No path to extend */
        if (bushes->SPcost_par[t_id][i] + network->arcs[ij].cost < bushes->SPcost_par[t_id][j]
            && bushes->LPcost_par[t_id][i] < bushes->LPcost_par[t_id][j]
            && (network->arcs[ij].tail == origin
                || network->arcs[ij].tail >= network->firstThroughNode))
        {
            bushes->flow_par[t_id][ij] = NEW_LINK;
            newArcs++;
            /* Never delete shortest path tree... should be OK with floating point
             * comparison since this is how SPcost is calculated */
        } else if (bushes->SPcost_par[t_id][i]+network->arcs[ij].cost==bushes->SPcost_par[t_id][j]
                   && bushes->flow_par[t_id][ij] == 0) {
            bushes->flow_par[t_id][ij] = NEW_LINK;
        }
    }

    /* If strict criterion fails, try a looser one */
    if (newArcs == 0) {
        for (ij = 0; ij < network->numArcs; ij++) {
            i = network->arcs[ij].tail;
            j = network->arcs[ij].head;
            if (bushes->LPcost_par[t_id][i]==-INFINITY && bushes->LPcost_par[t_id][j]>-INFINITY)
                continue; /* No path to extend */
            if (bushes->flow_par[t_id][ij] == 0 && bushes->LPcost_par[t_id][i] < bushes->LPcost_par[t_id][j]
                && (network->arcs[ij].tail == origin
                    || network->arcs[ij].tail >= network->firstThroughNode))
            {
                bushes->flow_par[t_id][ij] = NEW_LINK;
            }
        }
    }

    /* Finally update bush data structures: delete/add merges, find a new
     * topological order, rectify approach proportions */
    reconstructMerges_par(origin, network, bushes, t_id);
    parameters->topologicalOrder(origin, network, bushes, parameters);

//    displayMessage(FULL_NOTIFICATIONS, "End of update bushb %d\n", t_id);

}
/*
 * reconstructMerges -- update the merge data structures after a bush is
 * updated by adding and removing links.  In particular, for each node we do
 * the following...
 * 
 * 1. Normalize approach proportions for incoming links -- if there are zero
 *    such links for a non-origin, throw an error.
 * 2. Set pred if there is just one incoming link, otherwise create a merge.
 * 
 * Stores merges in a linked list at first, so they can be created in one pass.
 * Then transfer into an array for fast indexing.
 */
void reconstructMerges(int origin, network_type *network, bushes_type *bushes){
    int i, hi, lastApproach, m, arc, numApproaches;
    arcListElt *curArc;
    merge_type *merge;
    mergeDLL *mergeList = createMergeDLL();
    mergeDLLelt *curMerge;

   
    /* Create necessary merges */
    for (i = 0; i < network->numNodes; i++) {
        if (i == origin) continue;
        numApproaches = 0;
        for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                curArc = curArc->next) {
            hi = ptr2arc(network, curArc->arc);
            if (bushes->flow[hi] > 0 || bushes->flow[hi] == NEW_LINK) {
                numApproaches++;
                lastApproach = hi;
            }
        }
        if (numApproaches == 0)
            fatalError("Cannot have non-origin node %d in bush %d without"
                       "incoming contributing links", i, origin);
        if (numApproaches == 1) { /* No merge */
            bushes->pred[origin][i] = lastApproach;
        } else { /* Must create a merge */
            merge = createMerge(numApproaches);
            arc = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                 curArc = curArc->next) {
                hi = ptr2arc(network, curArc->arc);
                if (bushes->flow[hi] > 0 || bushes->flow[hi] == NEW_LINK) {
                    if (bushes->flow[hi] == NEW_LINK) bushes->flow[hi] = 0;
                    merge->approach[arc] = hi;
                    merge->approachFlow[arc] = bushes->flow[hi];
                    arc++;               
                }
            }
            insertMergeDLL(mergeList, merge, i, mergeList->tail);
        }
    }

    /* Now transfer to array, deleting old merges and replacing with new ones*/
    for (m = 0; m < bushes->numMerges[origin]; m++) {
        deleteMerge(bushes->merges[origin][m]);
    }
    deleteVector(bushes->merges[origin]);
    bushes->numMerges[origin] = mergeList->size;
    bushes->merges[origin] = newVector(mergeList->size, merge_type *);
    m = 0;
    for (curMerge = mergeList->head; curMerge != NULL;
         curMerge = curMerge->next) {
        bushes->merges[origin][m] = curMerge->merge;
        bushes->pred[origin][curMerge->node] = merge2pred(m);
        m++;
    }

    deleteMergeDLL(mergeList);
}

/*
 * reconstructMerges_par -- update the merge data structures after a bush is
 * updated by adding and removing links.  In particular, for each node we do
 * the following...
 *
 * 1. Normalize approach proportions for incoming links -- if there are zero
 *    such links for a non-origin, throw an error.
 * 2. Set pred if there is just one incoming link, otherwise create a merge.
 *
 * Stores merges in a linked list at first, so they can be created in one pass.
 * Then transfer into an array for fast indexing.
 */
void reconstructMerges_par(int origin, network_type *network, bushes_type *bushes, int t_id){
//    displayMessage(FULL_NOTIFICATIONS, "Top of update reconstructmerg %d\n", t_id);

    int i, hi, lastApproach, m, arc, numApproaches;
    arcListElt *curArc;
    merge_type *merge;
    mergeDLL *mergeList = createMergeDLL();
    mergeDLLelt *curMerge;


    /* Create necessary merges */
    for (i = 0; i < network->numNodes; i++) {
        if (i == origin) continue;
        numApproaches = 0;
        for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
             curArc = curArc->next) {
            hi = ptr2arc(network, curArc->arc);
            if (bushes->flow_par[t_id][hi] > 0 || bushes->flow_par[t_id][hi] == NEW_LINK) {
                numApproaches++;
                lastApproach = hi;
            }
        }
        if (numApproaches == 0)
            fatalError("Cannot have non-origin node %d in bush %d without"
                       "incoming contributing links", i, origin);
        if (numApproaches == 1) { /* No merge */
            bushes->pred[origin][i] = lastApproach;
        } else { /* Must create a merge */
            merge = createMerge(numApproaches);
            arc = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                 curArc = curArc->next) {
                hi = ptr2arc(network, curArc->arc);
                if (bushes->flow_par[t_id][hi] > 0 || bushes->flow_par[t_id][hi] == NEW_LINK) {
                    if (bushes->flow_par[t_id][hi] == NEW_LINK) bushes->flow_par[t_id][hi] = 0;
                    merge->approach[arc] = hi;
                    merge->approachFlow[arc] = bushes->flow_par[t_id][hi];
                    arc++;
                }
            }
            insertMergeDLL(mergeList, merge, i, mergeList->tail);
        }
    }

    /* Now transfer to array, deleting old merges and replacing with new ones*/
    for (m = 0; m < bushes->numMerges[origin]; m++) {
        deleteMerge(bushes->merges[origin][m]);
    }
    deleteVector(bushes->merges[origin]);
    bushes->numMerges[origin] = mergeList->size;
    bushes->merges[origin] = newVector(mergeList->size, merge_type *);
    m = 0;
    for (curMerge = mergeList->head; curMerge != NULL;
         curMerge = curMerge->next) {
        bushes->merges[origin][m] = curMerge->merge;
        bushes->pred[origin][curMerge->node] = merge2pred(m);
        m++;
    }

    deleteMergeDLL(mergeList);

//    displayMessage(FULL_NOTIFICATIONS, "End of update reconstructmerg %d\n", t_id);

}

/* 
 * findDivergenceNodes -- For a particular bush, find a divergence node
 * corresponding to every merge node.  For AlgorithmB, this means tracing the
 * longest and shortest path trees backward until they intersect.  Setting this
 * node once after bushes are created saves time during flow shifts.
 */
void findDivergenceNodes(int origin, network_type *network,
                         bushes_type *bushes) {
    int i, m, n, SPnode, LPnode, diverge;
    bool found_diverge;
    merge_type *merge;
    declareVector(int, mark, network->numNodes);

    for (i = 0; i < network->numNodes; i++) mark[i] = NO_PATH_EXISTS;

    for (m = 0; m < bushes->numMerges[origin]; m++) {
        i = network->arcs[bushes->merges[origin][m]->approach[0]].head;
        found_diverge = FALSE;
        SPnode = i;
        LPnode = i;
        /* Sometimes there is no longest used path to a node (zero demand),
         * nothing to do. */
        if (bushes->merges[origin][m]->LPlink == NO_PATH_EXISTS) {
            bushes->merges[origin][m]->divergenceNode = NO_PATH_EXISTS;
            continue;
        }
        while (found_diverge == FALSE) {
            if (SPnode != origin) {
                if (isMergeNode(origin, SPnode, bushes) == FALSE) {
                    SPnode = network->arcs[bushes->pred[origin][SPnode]].tail;
                } else {
                    n = pred2merge(bushes->pred[origin][SPnode]);
                    merge = bushes->merges[origin][n];
                    SPnode=network->arcs[merge->approach[merge->SPlink]].tail;
                }
                if (mark[SPnode] == i) {
                    diverge = SPnode;
                    found_diverge = TRUE;
                }
                mark[SPnode] = i;
            }
            if (LPnode != origin) {
                if (isMergeNode(origin, LPnode, bushes) == FALSE) {
                    LPnode = network->arcs[bushes->pred[origin][LPnode]].tail;
                } else {
                    n = pred2merge(bushes->pred[origin][LPnode]);
                    merge = bushes->merges[origin][n];               
                    LPnode=network->arcs[merge->approach[merge->LPlink]].tail;
                }
                if (mark[LPnode] == i) {
                    diverge = LPnode;
                    found_diverge = TRUE;
                }
                mark[LPnode] = i;
            }
        }
        bushes->merges[origin][m]->divergenceNode = diverge;
    }

    deleteVector(mark);
}


/*
 * updateFlowsB -- shift flows on a bush to move closer to equilibrium.
 * This function is fairly short, as it mainly just manages a few other
 * functions that do the actual work.
 * - calculateBushFlows gives current flows (remember, these are not stored
 *   persistently).
 * - rescanAndCheck both updates longest/shortest path labels, but also checks
 *   whether these are different enough to actually do anything -- depending
 *   on parameters we can skip bushes which are closer to equilibrium, to
 *   spend more time on bushes further away.
 * - updateFlowPass then does the actual flow shifting, inside of a loop in
 *   case we want to do multiple shifts per origin.
 */
void updateFlowsB(int origin, network_type *network, bushes_type *bushes,
                  algorithmBParameters_type *parameters) {
    int i;

    /* Recompute bush flows for this origin */
    calculateBushFlows(origin, network, bushes);

    /* Update longest/shortest paths, check whether there is work to do */
    if (rescanAndCheck(origin, network, bushes, parameters) == FALSE) return;

    /* Now do (possibly multiple) shifts per origin */
    for (i = 0; i < parameters->shiftReps; i++) {
        updateFlowPass(origin, network, bushes, parameters);
        /* Uncomment next line for extra validation checking */
        /* checkFlows(network, bushes); */
        if (parameters->rescanAfterShift == TRUE
                && i + 1 < parameters->shiftReps) {
            if (rescanAndCheck(origin, network, bushes, parameters) == FALSE)
                return;
        }
    }

}


/*
 * updateFlowsB -- shift flows on a bush to move closer to equilibrium.
 * This function is fairly short, as it mainly just manages a few other
 * functions that do the actual work.
 * - calculateBushFlows gives current flows (remember, these are not stored
 *   persistently).
 * - rescanAndCheck both updates longest/shortest path labels, but also checks
 *   whether these are different enough to actually do anything -- depending
 *   on parameters we can skip bushes which are closer to equilibrium, to
 *   spend more time on bushes further away.
 * - updateFlowPass then does the actual flow shifting, inside of a loop in
 *   case we want to do multiple shifts per origin.
 */
void updateFlowsB_par(int origin, network_type *network, bushes_type *bushes,
                  algorithmBParameters_type *parameters, int t_id) {
    int i;

    /* Recompute bush flows for this origin */
    calculateBushFlows_par(origin, network, bushes, t_id);

    /* Update longest/shortest paths, check whether there is work to do */
    if (rescanAndCheck_par(origin, network, bushes, parameters, t_id) == FALSE) return;

    /* Now do (possibly multiple) shifts per origin */
    for (i = 0; i < parameters->shiftReps; i++) {
        updateFlowPass_par(origin, network, bushes, parameters, t_id);
        /* Uncomment next line for extra validation checking */
//        checkFlows_par(network, bushes, t_id);
        if (parameters->rescanAfterShift == TRUE
            && i + 1 < parameters->shiftReps) {
            if (rescanAndCheck_par(origin, network, bushes, parameters, t_id) == FALSE)
                return;
        }
    }

}

/*
 * rescanAndCheck -- calls functions updating the longest/shortest path labels
 * on the bush, along with divergence nodes; and then checks whether the bush
 * is close enough to equilibrium to skip for now (in which case the function
 * returns FALSE).  If we need to shift flows, it returns TRUE.
 */
bool rescanAndCheck(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    int i;
    double maxgap = 0;

    scanBushes(origin, network, bushes, parameters, TRUE);
    findDivergenceNodes(origin, network, bushes);
    for (i = 0; i < network->numNodes; i++) {
        maxgap = max(maxgap, fabs(bushes->LPcost[i] -bushes->SPcost[i]));
    }
    if (maxgap < parameters->thresholdGap) return FALSE;

    return TRUE;
}

/*
 * rescanAndCheck -- calls functions updating the longest/shortest path labels
 * on the bush, along with divergence nodes; and then checks whether the bush
 * is close enough to equilibrium to skip for now (in which case the function
 * returns FALSE).  If we need to shift flows, it returns TRUE.
 */
bool rescanAndCheck_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, int t_id) {
//    displayMessage(FULL_NOTIFICATIONS, "Top of rescan and check\n");

    int i;
    double maxgap = 0;

    scanBushes_par(origin, network, bushes, parameters, TRUE, t_id);
    findDivergenceNodes(origin, network, bushes);
    for (i = 0; i < network->numNodes; i++) {
        maxgap = max(maxgap, fabs(bushes->LPcost_par[t_id][i] -bushes->SPcost_par[t_id][i]));
    }
    if (maxgap < parameters->thresholdGap) return FALSE;

    return TRUE;
}

/*
 * updateFlowPass -- does the actual work of shifting flows from longest to
 * shortest paths, through a pass in descending topological order. 
 */
void updateFlowPass(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    int i, j, k, m, node;
    merge_type *merge;

    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++)
        bushes->nodeFlow[i] = network->OD[origin][i].demand;

    for (; i < network->numNodes; i++) bushes->nodeFlow[i] = 0;

    /* Descending pass for flow shifts */
    for (node = network->numNodes - 1; node > 0; node--) {
        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple(j, origin, network, bushes);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            if (merge->LPlink == merge->SPlink
                || fabs(bushes->LPcost[j] - bushes->SPcost[j]
                        < parameters->minCostDifference)
                || bushes->nodeFlow[j] < parameters->minLinkFlowShift
                || merge->LPlink == NO_PATH_EXISTS)
            { /* Also nothing to do */
                pushBackFlowMerge(merge, network, bushes);
            } else {
                for (k = 0; k < parameters->numNewtonShifts; k++) {
                    newtonFlowShift(j,merge,origin,network,bushes,parameters);
                }
            }
        }
    }
}

/*
 * updateFlowPass_par -- does the actual work of shifting flows from longest to
 * shortest paths, through a pass in descending topological order.
 */
void updateFlowPass_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, int t_id) {
    int i, j, k, m, node;
    merge_type *merge;

    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++)
        bushes->nodeFlow_par[t_id][i] = network->OD[origin][i].demand;

    for (; i < network->numNodes; i++) bushes->nodeFlow_par[t_id][i] = 0;

    /* Descending pass for flow shifts */
    for (node = network->numNodes - 1; node > 0; node--) {
        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple_par(j, origin, network, bushes, t_id);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            if (merge->LPlink == merge->SPlink
                || fabs(bushes->LPcost_par[t_id][j] - bushes->SPcost_par[t_id][j]
                        < parameters->minCostDifference)
                || bushes->nodeFlow_par[t_id][j] < parameters->minLinkFlowShift
                || merge->LPlink == NO_PATH_EXISTS)
            { /* Also nothing to do */
                pushBackFlowMerge_par(merge, network, bushes, t_id);
            } else {
                newtonFlowShift_par(j,merge,origin,network,bushes,parameters, t_id);
            }
        }
    }
}

/* 
 * calculateBushFlows -- a key feature of this implementation is that bush
 * flows are not persistently stored, but generated as needed.  This saves
 * memory at the expense of a little extra computation.  calculateBushFlows
 * does this by making a pass in descending topological order, loading flows
 * from the OD matrix, and splitting flows as needed.
 */
void calculateBushFlows(int origin,network_type *network,bushes_type *bushes) {
    int i, j, ij, m, node;
    merge_type *merge;

    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++) {
        bushes->nodeFlow[i] = network->OD[origin][i].demand;
    }
    for (; i < network->numNodes; i++) {
        bushes->nodeFlow[i] = 0;
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        bushes->flow[ij] = 0;
    }

    /* Descending pass for flow calculations  */
    for (node = network->numNodes - 1; node > 0; node--) {
        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple(j, origin, network, bushes);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            pushBackFlowMerge(merge, network, bushes);
        }
    }

}

/*
 * calculateBushFlows_par -- a key feature of this implementation is that bush
 * flows are not persistently stored, but generated as needed.  This saves
 * memory at the expense of a little extra computation.  calculateBushFlows
 * does this by making a pass in descending topological order, loading flows
 * from the OD matrix, and splitting flows as needed.
 */
void calculateBushFlows_par(int origin,network_type *network,bushes_type *bushes, int t_id) {

    int i, j, ij, m, node;
    merge_type *merge;


    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++) {
        bushes->nodeFlow_par[t_id][i] = network->OD[origin][i].demand;
    }
    for (; i < network->numNodes; i++) {
        bushes->nodeFlow_par[t_id][i] = 0;
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        bushes->flow_par[t_id][ij] = 0;
    }

    /* Descending pass for flow calculations  */
    for (node = network->numNodes - 1; node > 0; node--) {

        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple_par(j, origin, network, bushes, t_id);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            pushBackFlowMerge_par(merge, network, bushes, t_id);
        }
    }

}

/*
 * pushBackFlowSimple -- the easy way to "split" flow.  For a non-merge node,
 * just push flow onto the predecessor.
 */
void pushBackFlowSimple(int j, int origin, network_type *network,
                        bushes_type *bushes) {
    int i, ij;

    ij = bushes->pred[origin][j];
    i = network->arcs[ij].tail;
    bushes->flow[ij] = bushes->nodeFlow[j];
    bushes->nodeFlow[i] += bushes->nodeFlow[j];
}

/*
 * pushBackFlowSimple_par -- the easy way to "split" flow.  For a non-merge node,
 * just push flow onto the predecessor.
 */
void pushBackFlowSimple_par(int j, int origin, network_type *network,
                        bushes_type *bushes, int t_id) {

//    displayMessage(FULL_NOTIFICATIONS, "top pushback simple\n");

    int i, ij;

    ij = bushes->pred[origin][j];
    i = network->arcs[ij].tail;
    bushes->flow_par[t_id][ij] = bushes->nodeFlow_par[t_id][j];
    bushes->nodeFlow_par[t_id][i] += bushes->nodeFlow_par[t_id][j];


}

/*
 * pushBackFlowMerge -- the harder way to split flow, when there are multiple
 * approaches to a merge node.
 */
void pushBackFlowMerge(merge_type *merge, network_type *network,
                       bushes_type *bushes) {
    int i, ij, arc;
    double flow;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        ij = merge->approach[arc];
        i = network->arcs[ij].tail;
        flow = merge->approachFlow[arc];
        bushes->flow[ij] = flow;
        bushes->nodeFlow[i] += flow;
    }
}

/*
 * pushBackFlowMerge_par -- the harder way to split flow, when there are multiple
 * approaches to a merge node.
 */
void pushBackFlowMerge_par(merge_type *merge, network_type *network,
                       bushes_type *bushes, int t_id) {

//    displayMessage(FULL_NOTIFICATIONS, "top pushback merge\n");

    int i, ij, arc;
    double flow;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        ij = merge->approach[arc];
        i = network->arcs[ij].tail;
        flow = merge->approachFlow[arc];
        bushes->flow_par[t_id][ij] = flow;
        bushes->nodeFlow_par[t_id][i] += flow;
    }
}

/*
 * rectifyMerge -- to guard against numerical errors at a particular merge
 * node, recalculates the approach proportions.  In the degenerate case of zero
 * node flow, push everything onto the shortest path. 
 */
void rectifyMerge(int j, merge_type *merge, bushes_type *bushes) {
    int arc;
    double totalFlow = 0;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        totalFlow += merge->approachFlow[arc];
    }

    if (totalFlow > 0) {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] *= bushes->nodeFlow[j] / totalFlow;
        }
    } else {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] = 0;
        }
        merge->approachFlow[merge->SPlink] = bushes->nodeFlow[j];
    }

}

/*
 * rectifyMerge_par -- to guard against numerical errors at a particular merge
 * node, recalculates the approach proportions.  In the degenerate case of zero
 * node flow, push everything onto the shortest path.
 */
void rectifyMerge_par(int j, merge_type *merge, bushes_type *bushes, int t_id) {

    int arc;
    double totalFlow = 0;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        totalFlow += merge->approachFlow[arc];
    }

    if (totalFlow > 0) {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] *= bushes->nodeFlow_par[t_id][j] / totalFlow;
        }
    } else {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] = 0;
        }
        merge->approachFlow[merge->SPlink] = bushes->nodeFlow_par[t_id][j];
    }

}

/*
 * newtonFlowShift -- For a given merge node, use Newton's method to shift flow
 * from the longer segment to the shorter one.
 */
void newtonFlowShift(int j, merge_type *merge, int origin,
                     network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters) {
    double flow1, flow2, cost1, cost2, der1, der2, shift;
    int i, hi, m;
    merge_type *segmentMerge;

    /* Calculate information for longest path segment */
    i = j;
    cost1 = 0;
    der1 = 0;
    flow1 = merge->approachFlow[merge->LPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            flow1=min(flow1,segmentMerge->approachFlow[segmentMerge->LPlink]);
        }
        cost1 += network->arcs[hi].cost;
        der1 += network->arcs[hi].der;
        i = network->arcs[hi].tail;
    }

   /* Do same for shortest path segment */
    i = j;
    cost2 = 0;
    der2 = 0;
    flow2 = merge->approachFlow[merge->SPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            flow2=min(flow2,segmentMerge->approachFlow[segmentMerge->SPlink]);
        }
        cost2 += network->arcs[hi].cost;
        der2 += network->arcs[hi].der;
        i = network->arcs[hi].tail;      
    }

    /* Determine Newton shift, truncating if need be to ensure feasibility */
    if (der1 + der2 == 0) der1 = parameters->minDerivative;
    shift = parameters->newtonStep * (cost1 - cost2) / (der1 + der2);
    shift = min(shift, flow1);
    shift = max(shift, -flow2);

    /* Now apply shift */
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            segmentMerge->approachFlow[segmentMerge->LPlink] -= shift;
        }
        bushes->flow[hi] -= shift;
        parameters->linkShiftB(hi, -shift, network);
        i = network->arcs[hi].tail;      
    }
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            segmentMerge->approachFlow[segmentMerge->SPlink] += shift;         
        }
        bushes->flow[hi] += shift;
        parameters->linkShiftB(hi, shift, network);
        i = network->arcs[hi].tail;      
    }
}

/*
 * newtonFlowShift -- For a given merge node, use Newton's method to shift flow
 * from the longer segment to the shorter one.
 */
void newtonFlowShift_par(int j, merge_type *merge, int origin,
                     network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters, int t_id) {
    double flow1, flow2, cost1, cost2, der1, der2, shift;
    int i, hi, m;
    merge_type *segmentMerge;

    /* Calculate information for longest path segment */
    i = j;
    cost1 = 0;
    der1 = 0;
    flow1 = merge->approachFlow[merge->LPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            flow1=min(flow1,segmentMerge->approachFlow[segmentMerge->LPlink]);
        }
        cost1 += network->arcs[hi].cost;
        der1 += network->arcs[hi].der;
        i = network->arcs[hi].tail;
    }

    /* Do same for shortest path segment */
    i = j;
    cost2 = 0;
    der2 = 0;
    flow2 = merge->approachFlow[merge->SPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            flow2=min(flow2,segmentMerge->approachFlow[segmentMerge->SPlink]);
        }
        cost2 += network->arcs[hi].cost;
        der2 += network->arcs[hi].der;
        i = network->arcs[hi].tail;
    }

    /* Determine Newton shift, truncating if need be to ensure feasibility */
    if (der1 + der2 == 0) der1 = parameters->minDerivative;
    shift = parameters->newtonStep * (cost1 - cost2) / (der1 + der2);
    shift = min(shift, flow1);
    shift = max(shift, -flow2);

    /* Now apply shift */
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            segmentMerge->approachFlow[segmentMerge->LPlink] -= shift;
        }
        bushes->flow_par[t_id][hi] -= shift;
//        bushes->nodeFlowshift_par[origin][hi] -= shift;
//        displayMessage(FULL_NOTIFICATIONS, "LP Bush flow shift for origin %d and arc %d is %f\n", origin, bushes->nodeFlowshift_par[origin][hi]);
//        parameters->linkShiftB(hi, -shift, network);
        exactCostUpdate_par(hi, -shift, network);
        i = network->arcs[hi].tail;
    }
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            segmentMerge->approachFlow[segmentMerge->SPlink] += shift;
        }
        bushes->flow_par[t_id][hi] += shift;
//        bushes->nodeFlowshift_par[origin][hi] += shift;
//        displayMessage(FULL_NOTIFICATIONS, "SP Bush flow shift for origin %d and arc %d is %f\n", origin, hi ,bushes->nodeFlowshift_par[origin][hi]);
        exactCostUpdate_par(hi, shift, network);

//        parameters->linkShiftB(hi, shift, network);
        i = network->arcs[hi].tail;
    }

}

/*
 * checkFlows -- an error-checking routine which can be invoked to see if flow
 * conservation is properly maintained on the bushes.  In the release
 * implementation, all calls to this function are commented out.
 */
#define FLOW_TOLERANCE 1e-5
void checkFlows(network_type *network, bushes_type *bushes) {
    int r, i, ij, kl;
    arcListElt *curArc;
    double balance;
   
    declareVector(double, flowCheck, network->numArcs);
    for (ij = 0; ij < network->numArcs; ij++) {
        flowCheck[ij] = 0;
    }
   
    for (r = 0; r < network->numZones; r++) {
        calculateBushFlows(r, network, bushes);
        /* First check bush flow consistency */
        for (ij = 0; ij < network->numArcs; ij++) {
            flowCheck[ij] += bushes->flow[ij];
            if (bushes->flow[ij] < 0) {
                for (kl = 0; kl < network->numArcs; kl++) {            
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                   network->arcs[kl].tail + 1,
                                   network->arcs[kl].head + 1,
                                   bushes->flow[kl]);
                }            
                fatalError("Flow validation failed: origin %d, link (%ld,%ld)"
                           "has negative flow.", r, network->arcs[ij].tail + 1,
                           network->arcs[ij].head + 1);
            }
        }
        for (i = 0; i < network->numNodes; i++) {
            if (i == r) continue;
            balance = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                    curArc = curArc->next) {
                balance += bushes->flow[ptr2arc(network, curArc->arc)];
            }
            for (curArc = network->nodes[i].forwardStar.head; curArc != NULL;
                    curArc = curArc->next) {
                balance -= bushes->flow[ptr2arc(network, curArc->arc)];
            }         
            if (i < network->numZones) balance -= network->OD[r][i].demand;
            if (fabs(balance) > FLOW_TOLERANCE) {
                for (kl = 0; kl < network->numArcs; kl++) {            
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                    network->arcs[kl].tail + 1,
                                    network->arcs[kl].head + 1,
                                    bushes->flow[kl]);
                }            
                fatalError("Flow validation failed: origin %d, node %ld"
                           "violates conservation.", r, i + 1);
         
            }
        }
      
    }
   
    /* Then check overall link flow consistency */
    for (ij = 0; ij < network->numArcs; ij++) {
        if (fabs(flowCheck[ij] - network->arcs[ij].flow) > FLOW_TOLERANCE) {
            fatalError("Flow validation failed: link (%ld,%ld) has aggregate"
                       "flow %f, but %f when disaggregated.",
                       network->arcs[ij].tail + 1,
                       network->arcs[ij].head + 1,
                       network->arcs[ij].flow, flowCheck[ij]);
        }
    }
   
    deleteVector(flowCheck);
}

/*
 * checkFlows_par -- an error-checking routine which can be invoked to see if flow
 * conservation is properly maintained on the bushes.  In the release
 * implementation, all calls to this function are commented out.
 */
void checkFlows_par(network_type *network, bushes_type *bushes, int t_id) {
    int r, i, ij, kl;
    arcListElt *curArc;
    double balance;

    declareVector(double, flowCheck, network->numArcs);
    for (ij = 0; ij < network->numArcs; ij++) {
        flowCheck[ij] = 0;
    }

    for (r = 0; r < network->numZones; r++) {
        calculateBushFlows_par(r, network, bushes, t_id);
        /* First check bush flow consistency */
        for (ij = 0; ij < network->numArcs; ij++) {
            flowCheck[ij] += bushes->flow_par[t_id][ij];
            if (bushes->flow_par[t_id][ij] < 0) {
                for (kl = 0; kl < network->numArcs; kl++) {
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                   network->arcs[kl].tail + 1,
                                   network->arcs[kl].head + 1,
                                   bushes->flow[kl]);
                }
                fatalError("Flow validation failed: origin %d, link (%ld,%ld)"
                           "has negative flow.", r, network->arcs[ij].tail + 1,
                           network->arcs[ij].head + 1);
            }
        }
        for (i = 0; i < network->numNodes; i++) {
            if (i == r) continue;
            balance = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                 curArc = curArc->next) {
                balance += bushes->flow_par[t_id][ptr2arc(network, curArc->arc)];
            }
            for (curArc = network->nodes[i].forwardStar.head; curArc != NULL;
                 curArc = curArc->next) {
                balance -= bushes->flow_par[t_id][ptr2arc(network, curArc->arc)];
            }
            if (i < network->numZones) balance -= network->OD[r][i].demand;
            if (fabs(balance) > FLOW_TOLERANCE) {
                for (kl = 0; kl < network->numArcs; kl++) {
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                   network->arcs[kl].tail + 1,
                                   network->arcs[kl].head + 1,
                                   bushes->flow_par[t_id][kl]);
                }
                fatalError("Flow validation failed: origin %d, node %ld"
                           "violates conservation.", r, i + 1);

            }
        }

    }

    /* Then check overall link flow consistency */
    for (ij = 0; ij < network->numArcs; ij++) {
        if (fabs(flowCheck[ij] - network->arcs[ij].flow) > FLOW_TOLERANCE) {
            fatalError("Flow validation failed: link (%ld,%ld) has aggregate"
                       "flow %f, but %f when disaggregated.",
                       network->arcs[ij].tail + 1,
                       network->arcs[ij].head + 1,
                       network->arcs[ij].flow, flowCheck[ij]);
        }
    }

    deleteVector(flowCheck);
}

/*
 * exactCostUpdate -- The most precise way to update the cost on a link after
 * changing its flow, by explicitly recomputing the BPR function and its
 * derivative.
 */
void exactCostUpdate(int ij, double shift, network_type *network) {
//    pthread_mutex_lock(&flow_mut);
    network->arcs[ij].flow += shift;
    network->arcs[ij].cost=network->arcs[ij].calculateCost(&network->arcs[ij]);
    network->arcs[ij].der = network->arcs[ij].calculateDer(&network->arcs[ij]);
//    pthread_mutex_unlock(&flow_mut);
}

void exactCostUpdate_par(int ij, double shift, network_type *network) {
    pthread_mutex_lock(&network->arc_muts[ij]);
    network->arcs[ij].flow += shift;
    network->arcs[ij].cost=network->arcs[ij].calculateCost(&network->arcs[ij]);
    network->arcs[ij].der = network->arcs[ij].calculateDer(&network->arcs[ij]);
    pthread_mutex_unlock(&network->arc_muts[ij]);
}

/*
 * linearCostUpdate -- A faster way to approximately update costs, using a
 * linear approximation to the BPR function (and keeping the derivative
 * unchanged.)
 */
void linearCostUpdate(int ij, double shift, network_type *network) {
    network->arcs[ij].flow += shift;
    network->arcs[ij].cost += shift * network->arcs[ij].der;
}

/*
 * noCostUpdate -- The laziest way to "update" costs -- do nothing except shift
 * the flow.  Be careful to update the cost explicitly somewhere else.
 */
void noCostUpdate(int ij, double shift, network_type *network) {
    network->arcs[ij].flow += shift;
}

/************************
 ** MERGE LINKED LISTS **
 ************************/

merge_type *createMerge(int numApproaches) {
    declareScalar(merge_type, merge);
    merge->numApproaches = numApproaches;
    merge->approach = newVector(numApproaches, int);
    merge->approachFlow = newVector(numApproaches, double);
    merge->LPlink = NO_PATH_EXISTS;
    merge->SPlink = NO_PATH_EXISTS;
    merge->divergenceNode = NO_PATH_EXISTS;
    return merge;
}

void deleteMerge(merge_type *merge) {
    deleteVector(merge->approach);
    deleteVector(merge->approachFlow);
    deleteScalar(merge);
}

void displayMerge(int minVerbosity, merge_type *merge, network_type *network) {
    int a;
   
    if (verbosity < minVerbosity) return;
   
    displayMessage(minVerbosity, "Merge with %d approaches: \n",
                   merge->numApproaches);
    for (a = 0; a < merge->numApproaches; a++) {
        displayMessage(minVerbosity, "(%ld,%ld) %f %s %s\n",
                        network->arcs[merge->approach[a]].tail,
                        network->arcs[merge->approach[a]].head,
                        merge->approachFlow[a],
                        merge->SPlink == a ? "SP" : "  ",
                        merge->LPlink == a ? "LP" : "  ");
        }
    displayMessage(minVerbosity,"Divergence node: %d\n",merge->divergenceNode);
}

mergeDLL *createMergeDLL() {
    declareScalar(mergeDLL, newdll);
    newdll->head = NULL;
    newdll->tail = NULL;
    newdll->size = 0;
    return newdll;
}

mergeDLLelt *insertMergeDLL(mergeDLL *list, merge_type *merge, int i,
                            mergeDLLelt *after) {
    declareScalar(mergeDLLelt, newNode);
    newNode->merge = merge;
    newNode->node = i;
    if (after != NULL) {
        newNode->prev = after;
        newNode->next = after->next;
    if (list->tail != after)
        newNode->next->prev = newNode;
    else
        list->tail = newNode;
        after->next = newNode;
    } else {
        newNode->prev = NULL;
        newNode->next = list->head;
        if (list->tail != after)
            newNode->next->prev = newNode;
        else
            list->tail = newNode;
        list->head = newNode;
    }
    list->size++;
    return newNode;
}

void deleteMergeDLL(mergeDLL *list) {
    while (list->head != NULL)
        deleteMergeDLLelt(list, list->tail);
    killScalar(list);
}

void deleteMergeDLLelt(mergeDLL *list, mergeDLLelt *elt) {
    if (list->tail != elt) {
        if (list->head != elt)
            elt->prev->next = elt->next;
        else
            list->head = elt->next;
        elt->next->prev = elt->prev;
    } else {
        list->tail = elt->prev;
        if (list->head != elt)
            elt->prev->next = elt->next;
        else
            list->head = elt->next;
    }
    list->size--;
    killScalar(elt);
}

void displayMergeDLL(int minVerbosity, mergeDLL *list) {
    mergeDLLelt *curnode = list->head;
    displayMessage(minVerbosity,"Start of the list: %p\n", (void *)list->head);
    while (curnode != NULL) {
        displayMessage(minVerbosity, "%p %d %p %p\n", (void *)curnode,
                curnode->node, (void *)curnode->prev, (void *)curnode->next);
        curnode = (*curnode).next;
    }
    displayMessage(minVerbosity, "End of the list: %p\n", (void *)list->tail);
}
