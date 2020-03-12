/*
 * bush.c -- Contains all of the components of Algorithm B, including bush
 * creation, updating, and flow shifting.
 *
 * See comments on bush.h for descriptions of the data structures used for
 * bushes.
 */
#include "bush.h"
#include <time.h>

#if PARALLELISM
    #include "thpool.h"
    #include <pthread.h> 
    #include "parallel_bush.h"
#endif

#if PARALLELISM
//Struct for thread arguments
struct thread_args {
    int id;
    bool update_flows_ret;
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
    updateBushB_par(id, network, bushes, parameters);
}

void updateFlowsPool(void* pVoid) {
    struct thread_args *args = (struct thread_args *) pVoid;
    int id = args->id;
    bushes_type *bushes = args->bushes;
    network_type *network = args->network;
    algorithmBParameters_type *parameters = args->parameters;
    args->update_flows_ret |= updateFlowsB_par(id, network, bushes, parameters);
}
threadpool thpool;
#endif

/*
 * AlgorithmB -- master function controlling overall flow of the algorithm.
 * Arguments are pointers to a network, and to a struct of algorithm
 * parameters. 
 */
void AlgorithmB(network_type *network, algorithmBParameters_type *parameters) {
#if PARALLELISM
pthread_mutex_init(&shift_lock, NULL);
thpool = thpool_init(parameters->numThreads);
#endif
    /* Strong connectivity check */
    makeStronglyConnectedNetwork(network);

    /* Allocate memory for bushes */
    int batch, iteration = 0, lastClass = IS_MISSING;
    displayMessage(FULL_NOTIFICATIONS, "Creating Initial bushes\n");
    bushes_type *bushes = createBushes(network);
    struct timespec tick, tock;
    char beckmannString[STRING_SIZE];

    double elapsedTime = 0, gap = INFINITY, batchGap = INFINITY;

    /* Initialize */
    clock_t stopTime = clock(); /* used for timing */
    initializeAlgorithmB(network, &bushes, parameters);
    displayMessage(LOW_NOTIFICATIONS, "Initialization done in %.3f s.\n",
        ((double)(clock() - stopTime)) / CLOCKS_PER_SEC);
    if(parameters->calculateBeckmann == TRUE)
        displayMessage(DEBUG, "Initial Beckmann: %f\n", BeckmannFunction(network));

    
    /* Iterate */
    parameters->innerIterations = 1;
    while (elapsedTime < parameters->maxTime
             && iteration < parameters->maxIterations
             && (iteration == 0 || gap > parameters->convergenceGap)) {
        iteration++;
        gap = 0; /* Will accumulate total gap across batches for averaging */
        clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
        /* Iterate over batches of origins */
        for (batch = 0; batch < network->numBatches; batch++) {
            /* Do main work for this batch */
            displayMessage(FULL_NOTIFICATIONS, "Loading Batch...\n");
            loadBatch(batch, network, &bushes, parameters);
            displayMessage(FULL_NOTIFICATIONS, "Loaded Batch...\nUpdating Batch Bush...\n");
            updateBatchBushes(network, bushes, &lastClass, parameters);
            displayMessage(FULL_NOTIFICATIONS, "Updated Batch Bush...\nUpdating Batch Flows...\n");
            updateBatchFlows(network, bushes, &lastClass, parameters);
            displayMessage(FULL_NOTIFICATIONS, "Updated Batch Flows...\nStoring batch ...\n");
            storeBatch(batch, network, bushes, parameters);
            displayMessage(FULL_NOTIFICATIONS, "Stored Batch\n");

            /* Check gap and report progress. */
            clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
            elapsedTime += (double)((1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec)) * 1.0/1000000000; /* Exclude gap calculations from run time */
            stopTime = clock();
            displayMessage(FULL_NOTIFICATIONS, "Calculating batch relative gap...\n");

            batchGap = bushRelativeGap(network, bushes, parameters);
            displayMessage(FULL_NOTIFICATIONS, "Calculated batch relative gap...\n");
            gap += batchGap;
            if (parameters->includeGapTime == FALSE) stopTime = clock(); 
            if (parameters->calculateBeckmann == TRUE) {
                sprintf(beckmannString, "obj %.15g, ",
                        BeckmannFunction(network));
            } else {
                beckmannString[0] = '\0';
            }

            if (network->numBatches > 1) {
                displayMessage(LOW_NOTIFICATIONS, "*Batch %ld: batchgap %.15f,"
                       " %stime %.3f s.\n", batch, batchGap,
                           beckmannString, elapsedTime);  
            }
                    
        }
        /* Report information from the entire iteration (all batches) */
        gap /= network->numBatches;
        displayMessage(LOW_NOTIFICATIONS, "Iteration %d:%s gap %.15f, "
                       "%stime %.3f, %d shifts\n",
                        iteration, 
                        network->numBatches > 1 ? " est." : "",
                        gap,
                        beckmannString,
                        elapsedTime,
                        parameters->numFlowShifts);
        parameters->innerIterations = 20; /* Tuneable parameter */
        parameters->thresholdAEC = 0.25 * gap; /* Tuneable parameter */

    } 

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
    parameters.maxIterations = INT_MAX;

    parameters.innerIterations = 20;
    parameters.shiftReps = 1;
    parameters.rescanAfterShift = FALSE;

    parameters.thresholdGap = 0;
    parameters.thresholdAEC = 0;
    parameters.minCostDifference = 0;
    parameters.minLinkFlowShift = 0;

    parameters.minLinkFlow = 1e-14;
    parameters.minDerivative = 1e-6;
    parameters.newtonStep = 1;
    parameters.numNewtonShifts = 1;

    parameters.numFlowShifts = 0;

    parameters.demandMultiplier = 1;

    parameters.storeMatrices = FALSE;
    parameters.storeBushes = FALSE;
    parameters.reuseFirstBush = FALSE;

    snprintf(parameters.batchStem, sizeof(parameters.batchStem), "batch");
    snprintf(parameters.matrixStem, sizeof(parameters.matrixStem), "matrix");
    snprintf(parameters.flowsFile, sizeof(parameters.flowsFile), "flows.txt");

    parameters.warmStart = FALSE;

    parameters.SPQueueDiscipline = DEQUE;
    
    parameters.includeGapTime = TRUE;

    parameters.createInitialBush = &initialBushShortestPath;
    parameters.topologicalOrder = &genericTopologicalOrder;
    parameters.linkShiftB = &exactCostUpdate;

    return parameters;
}

void initializeAlgorithmB(network_type *network, bushes_type **bushes,
                          algorithmBParameters_type *parameters) {

    int batch, c, ij, origin;
    char batchFileName[STRING_SIZE];    
    displayMessage(FULL_NOTIFICATIONS, "Initializing Algorithm B\n");
    /* 1. Initialize flows to zero and costs to free-flow */
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].flow = 0;
        network->arcs[ij].cost = network->arcs[ij].freeFlowTime;
        for (c = 0; c < network->numClasses; c++) {
            network->arcs[ij].classFlow[c] = 0;
        }
    }
    displayMessage(FULL_NOTIFICATIONS, "Initialed flows\n");
    /* 2. Now create bushes, either reading from a file (if warmStart parameter
          is TRUE); or creating entirely from scratch (default); or by
          creating one set of bushes from scratch and reusing for other classes
          (can only use if reuseFirstBush parameter is TRUE *and* batch size
          is the same as the number of origins.  The following code handles
          all of these cases*/
    if (network->numBatches != network->numClasses &&
            parameters->reuseFirstBush == TRUE) {
        fatalError("Cannot reuse bushes unless batches match classes.");
    }
              
    for (batch = 0; batch < network->numBatches; batch++) {
        /* Set up new batch */
        network->curBatch = batch;
        sprintf(batchFileName, "%s%d.bin", parameters->batchStem,
                network->curBatch);

        if (network->numBatches > 1 || parameters->storeMatrices == TRUE) {
                displayMessage(FULL_NOTIFICATIONS, "Reading matrix %d\n", network->curBatch);
            readBinaryMatrix(network, parameters);
        }
        displayMessage(FULL_NOTIFICATIONS, "Read matrix %d\n", network->curBatch);

        /* Form bush structure: either read from file or recreate */
        if (parameters->warmStart == TRUE) { /* Read file and rectify */
            displayMessage(FULL_NOTIFICATIONS, "Reading batch %d\n", batch);
            readBushes(network, bushes, batchFileName);
            displayMessage(FULL_NOTIFICATIONS, "Read batch %d\n", batch);
        } else { /* No warm start, have to re-initialize */
            /* Do we have to create a bush from scratch, or can we reuse
             * the first? */
            displayMessage(FULL_NOTIFICATIONS, "Creating batch %d\n", batch);
            if (batch == 0 || parameters->reuseFirstBush == FALSE) {
                initializeBushesB(network, *bushes, parameters);
            }
            displayMessage(FULL_NOTIFICATIONS, "Created batch %d\n", batch);
        }
        /* Now load demand onto these bushes to calculate flows.  Need to do
         * this (regardless of fresh initialization) for two reasons:
         *  1. In case demand changes (warm starting)
         *  2. To get the initial arc flows correct
         */
        for (origin = 0; origin < network->batchSize; origin++) {
            rectifyBushFlows(origin, network, *bushes);
            c = origin2class(network, origin);
            for (ij = 0; ij < network->numArcs; ij++) {
                network->arcs[ij].flow += (*bushes)->flow[ij];
                network->arcs[ij].classFlow[c] += (*bushes)->flow[ij];
                network->arcs[ij].cost =
                    network->arcs[ij].calculateCost(&network->arcs[ij]);
            }
        }
        sprintf(batchFileName, "%s%d.bin", parameters->batchStem, batch);
        if (network->numBatches > 1 || parameters->storeBushes == TRUE) {
            writeBushes(network, *bushes, batchFileName);
        }
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].der =
            network->arcs[ij].calculateDer(&network->arcs[ij]);
    }    
}

void loadBatch(int batch, network_type *network, bushes_type **bushes,
               algorithmBParameters_type *parameters) {
    char batchFileName[STRING_SIZE];
    
    network->curBatch = batch;
    sprintf(batchFileName, "%s%d.bin", parameters->batchStem, batch);
    if (network->numBatches > 1) {
        /* Even if storeBushes == TRUE, if there is only one batch
         * there is no point in reading it again */
        readBushes(network, bushes, batchFileName);
    }
    if (network->numBatches > 1 || parameters->storeMatrices == TRUE) {
        readBinaryMatrix(network, parameters);
    }               
}               

void storeBatch(int batch, network_type *network, bushes_type *bushes,
                algorithmBParameters_type *parameters) {
    char batchFileName[STRING_SIZE];
    
    sprintf(batchFileName, "%s%d.bin", parameters->batchStem, batch);                
    if (network->numBatches > 1 || parameters->storeBushes == TRUE) {
        writeBushes(network, bushes, batchFileName);
    }
}

void updateBatchBushes(network_type *network, bushes_type *bushes,
                       int *lastClass, algorithmBParameters_type *parameters) {
                          
    
#if PARALLELISM
    int c; 
    struct thread_args args[network->batchSize];
    for (int j = 0; j < network->batchSize; ++j) {
        args[j].id = j;
        args[j].network = network;
        args[j].parameters = parameters;
        args[j].bushes = bushes;
        args[j].update_flows_ret = FALSE;
    }

    for (int j = 0; j < network->batchSize; ++j) {
        if (outOfOrigins(network, j) == TRUE) break;
        bushes->updateBush[j] = TRUE;
        c = origin2class(network, j);
        if (c != *lastClass) {
            changeFixedCosts(network, c);
        }
        thpool_add_work(thpool, (void (*)(void *)) updateBushPool, (void*)&args[j]);
        *lastClass = c;
    }
    thpool_wait(thpool);

    for (int j = 0; j < network->batchSize; ++j) {
        if (outOfOrigins(network, j) == TRUE) break;
        bushes->updateBush[j] = TRUE;
        c = origin2class(network, j);
        if (c != *lastClass) {
            changeFixedCosts(network, c);
        }
        thpool_add_work(thpool, (void (*)(void *)) updateFlowsPool, (void*)&args[j]);
        *lastClass = c;
    }
    thpool_wait(thpool);
#else
    int origin, c;
    for (origin = 0; origin < network->batchSize; origin++) {
        if (outOfOrigins(network, origin) == TRUE) break;
        bushes->updateBush[origin] = TRUE;
        c = origin2class(network, origin);
        if (c != *lastClass) {
            changeFixedCosts(network, c);
        }
        updateBushB(origin, network, bushes, parameters);
        updateFlowsB(origin, network, bushes, parameters);
        *lastClass = c;
    }
#endif
}

void updateBatchFlows(network_type *network, bushes_type *bushes,
                      int *lastClass, algorithmBParameters_type *parameters) {
    int i, c;                       
    bool doneAny;
    for (i = 0; i < parameters->innerIterations; i++) {
        doneAny = FALSE;
            
 #if PARALLELISM
         struct thread_args args[network->batchSize];
         for (int j = 0; j < network->batchSize; ++j) {
             args[j].id = j;
             args[j].network = network;
             args[j].parameters = parameters;
             args[j].bushes = bushes;
             args[j].update_flows_ret = FALSE;
         }
    
         for (int j = 0; j < network->batchSize; ++j) {
             if (outOfOrigins(network, j) == TRUE) break;
             if (bushes->updateBush[j] == FALSE) continue;
             c = origin2class(network, j);
             if (c != *lastClass) {
                 changeFixedCosts(network, c);
             }
             thpool_add_work(thpool, (void (*)(void *)) updateFlowsPool, (void*)&args[j]);
         }
         thpool_wait(thpool);

         for (int j = 0; j < network->batchSize; ++j) {
             doneAny |= args[j].update_flows_ret;
         }
#else
        int origin;
        for (origin = 0; origin < network->batchSize; origin++) {
            if (outOfOrigins(network, origin) == TRUE) break;
            if (bushes->updateBush[origin] == FALSE) continue;
            c = origin2class(network, origin);
            if (c != *lastClass) {
                changeFixedCosts(network, c);
            }
            doneAny |= updateFlowsB(origin,network,bushes,parameters);
            *lastClass = c;
        }
#endif
        if (doneAny == FALSE) break;
    }
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
    int ij, originNode = origin2node(network, origin);
    arcIndexBellmanFord(originNode, bushes->SPcost, bushes->pred[origin],
                        network, parameters->SPQueueDiscipline);
//    displayMessage(FULL_NOTIFICATIONS, "Finished shortest path bush preds %d\n", origin);
    for (ij = 0; ij < network->numNodes; ij++) {
        if (ij != originNode && bushes->pred[origin][ij] == IS_MISSING) {
            if (verbosity >= DEBUG) {
                for (ij = 0; ij < network->numNodes; ij++) {
                    displayMessage(DEBUG, "%d %d %f\n", ij + 1,
                                    bushes->pred[origin][ij],
                                    bushes->SPcost[ij]);
                }
                for (ij = 0; ij < network->numArcs; ij++) {
                    displayMessage(DEBUG, "(%d,%d) %f %f\n", 
                                   network->arcs[ij].tail+1,
                                   network->arcs[ij].head+1,
                                   network->arcs[ij].flow,
                                   network->arcs[ij].cost);
                }
            }
            fatalError("Cannot find initial connected bush.");
        }
    }

    /* Set bush topological order */
    bushes->numMerges[origin] = 0;
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
//    displayMessage(FULL_NOTIFICATIONS, "Starting topo order with indegree vector\n");
    declareVector(int, indegree, network->numNodes);
    for (i = 0; i < network->numNodes; i++) {
        indegree[i] = 1; /* By default non-origin nodes are assumed to have 1
                            incoming link; merges and origin handled below */
         bushes->bushOrder[origin][i] = NO_PATH_EXISTS;
    }
    /*displayMessage(DEBUG, "Now working with origin %d, really %d\n",
                  origin, origin2node(network, origin));*/
    for (i = 0; i < network->numNodes; i++) {
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            //displayMessage(LOW_NOTIFICATIONS, "Node %d is merge node %d\n", i, m);
            indegree[i] = bushes->merges[origin][m]->numApproaches;
        }
    }
    indegree[origin2node(network, origin)] = 0;

//    displayMessage(FULL_NOTIFICATIONS, "Making topo order q\n");
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
        displayMessage(LOW_NOTIFICATIONS, "next: %d, network->numNodes: %d\n", 
                           next, network->numNodes);
        fatalError("Graph given to bushTopologicalOrder contains a cycle.");
    }
    bushes->lastMerge[origin] = highestMerge;

//    displayMessage(FULL_NOTIFICATIONS, "Deleting topo order q\n");
    deleteQueue(&LIST);
    deleteVector(indegree);
//    displayMessage(FULL_NOTIFICATIONS, "Deleted topo order q\n");

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
    return (bushes->pred[origin][i] < 0 &&
            i != origin2node(bushes->network, origin)) ? TRUE : FALSE;
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

/***********************
 * CUSTOM GAP ROUTINES *
 ***********************/

/* 
 * bushSPTT -- specialized SPTT finding using bush structures as a warm start
 */
double bushSPTT(network_type *network, bushes_type *bushes,
              algorithmBParameters_type *parameters) {
    int r, j, c, lastClass = IS_MISSING, originNode;
    double sptt = 0;
    for (r = 0; r < network->batchSize; r++) {
        if (outOfOrigins(network, r) == TRUE) break;
        originNode = origin2node(network, r);
        c = origin2class(network, r);
        if (c != lastClass) {
            changeFixedCosts(network, c);
        }
        scanBushes(r, network, bushes, parameters, NO_LONGEST_PATH);
        BellmanFord_NoLabel(originNode, bushes->SPcost, network, DEQUE,
                            bushes->SPcost, bushes->bushOrder[r]);
        for (j = 0; j < network->numZones; j++) {
            sptt += network->demand[r][j] * bushes->SPcost[j];
        }
        lastClass = c;
    }
    return sptt;
}

double bushTSTT(network_type *network, bushes_type *bushes) {
    int r, ij, c, lastClass = IS_MISSING;
    double tstt = 0;
    for (r = 0; r < network->batchSize; r++) {
        if (outOfOrigins(network, r) == TRUE) break;
        c = origin2class(network, r);
        if (c != lastClass) {
            changeFixedCosts(network, c);
        }
        calculateBushFlows(r, network, bushes);
        for (ij = 0; ij < network->numArcs; ij++) {
            tstt += bushes->flow[ij] * network->arcs[ij].cost;
        }
        lastClass = c;
    }
    return tstt;
}

double bushRelativeGap(network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters) {
    
    double sptt = bushSPTT(network, bushes, parameters);
    double tstt = network->numBatches == 1 ?
                                         TSTT(network) :
                                         bushTSTT(network, bushes);
    displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest "
                          "path TSTT: %f\n", tstt, sptt);
    /* if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and "
                                "denom are %f %f\n", tstt, sptt); */
    if (sptt == 0 && tstt == 0) {
        warning(LOW_NOTIFICATIONS, "No flow or demand on bush\n");
        return 0;
    } else if(sptt == 0) {
        warning(LOW_NOTIFICATIONS, "SPTT is zero\n");
    }
    return (tstt / sptt - 1);
}

double bushAEC(network_type *network, bushes_type *bushes,
             algorithmBParameters_type *parameters) {
    double sptt = bushSPTT(network, bushes, parameters);
    double tstt = TSTT(network);
    if (tstt < sptt) warning(LOW_NOTIFICATIONS, "ANegative gap.  TSTT and SPTT "
                                                "are %f %f\n", tstt, sptt);
    if (network->totalODFlow == 0) {
        warning(LOW_NOTIFICATIONS, "No flow or demand on bush\n");
        return 0;
    }
    return ((tstt - sptt) / network->totalODFlow);
}

double bushMEC(network_type *network, bushes_type *bushes,
             algorithmBParameters_type *parameters) {
    double mec = 0;
    int j, r, c, lastClass = IS_MISSING, originNode;
    for (r = 0; r < network->batchSize; r++) {
        if (outOfOrigins(network, r) == TRUE) break;
        originNode = origin2node(network, r);
        c = origin2class(network, r);
        if (c != lastClass) {
            changeFixedCosts(network, c);
        }
        scanBushes(r, network, bushes, parameters, LONGEST_USED_PATH);
        BellmanFord_NoLabel(originNode, bushes->SPcost, network, DEQUE,
                            bushes->SPcost, bushes->bushOrder[r]);
        for (j = 0; j < network->numZones; j++) {
            mec = max(mec, bushes->LPcost[j] - bushes->SPcost[j]);
        }
    }
    return mec;
}


/*
 * createBushes -- Initialize the bushes data structure for a particular
 * network, allocating memory as needed.
 */
bushes_type *createBushes(network_type *network) {
    declareScalar(bushes_type, bushes);
    int i;

    bushes->network = network;
    bushes->LPcost = newVector(network->numNodes, double);
    bushes->SPcost = newVector(network->numNodes, double);
    bushes->flow = newVector(network->numArcs, double);
    bushes->nodeFlow = newVector(network->numNodes, double);
    bushes->pred = newMatrix(network->batchSize, network->numNodes, int);
    bushes->bushOrder = newMatrix(network->batchSize, network->numNodes, int);
    bushes->lastMerge = newVector(network->batchSize, int);
    bushes->numMerges = newVector(network->batchSize, int);
    bushes->merges = newVector(network->batchSize, merge_type**);
    bushes->updateBush = newVector(network->batchSize, bool);

#if PARALLELISM
    bushes->LPcost_par = newMatrix(network->batchSize, network->numNodes,double);
    bushes->SPcost_par = newMatrix(network->batchSize, network->numNodes,double);
    bushes->flow_par = newMatrix(network->batchSize, network->numArcs,double);
    bushes->nodeFlow_par = newMatrix(network->batchSize, network->numNodes,double);
#endif
    for (i = 0; i < network->batchSize; i++) {
       bushes->numMerges[i] = 0;
       bushes->merges[i] = newVector(1, merge_type*);
       bushes->updateBush[i] = TRUE;
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
    deleteMatrix(bushes->pred, network->batchSize);
    deleteMatrix(bushes->bushOrder, network->batchSize);
    deleteVector(bushes->lastMerge);
    for (i = 0; i < network->batchSize; i++) {
        for (m = 0; m < bushes->numMerges[i]; m++) {
            deleteMerge(bushes->merges[i][m]);
        }
        deleteVector(bushes->merges[i]);
    }
    deleteVector(bushes->numMerges);
    deleteVector(bushes->merges);
    deleteVector(bushes->updateBush);
#if PARALLELISM
    deleteMatrix(bushes->LPcost_par, network->batchSize);
    deleteMatrix(bushes->SPcost_par, network->batchSize);
    deleteMatrix(bushes->flow_par, network->batchSize);
    deleteMatrix(bushes->nodeFlow_par, network->batchSize);
#endif

    deleteScalar(bushes);
}


/*
 * initializeBushesB -- set up initial bushes at the beginning of the algorithm
 * (as distinct from createBushes, which merely allocates memory).  This
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
    int c, origin, lastClass = IS_MISSING;

    for (origin = 0; origin < network->batchSize; origin++) {
        if (outOfOrigins(network, origin) == TRUE) break;
        c = origin2class(network, origin);
//        displayMessage(FULL_NOTIFICATIONS, "Working on origin %d with class %d\n", origin, c);
        if (c != lastClass) {
            changeFixedCosts(network, c);
        }
        lastClass = c;
        /* createInitialBush also sets preds, bushOrder */
        parameters->createInitialBush(origin, network, bushes, parameters);
//        displayMessage(FULL_NOTIFICATIONS, "Created Initial bush %d\n", origin);
        calculateBushFlows(origin, network, bushes);
    }

}

/*
 * scanBushes -- simultaneously find shortest and longest paths for a given
 * bush.  If the longestUsed argument is TRUE, the longest path search will
 * restrict attention to longest *used* paths.  If FALSE, longest bush paths
 * will be calculated regardless of whether there is flow on them.
*/
void scanBushes(int origin, network_type *network, bushes_type *bushes,
                algorithmBParameters_type *parameters, scan_type LPrule) {
    int h, i, hi, m, curnode, curarc;
    double tempcost;
    merge_type *merge;

    for (i = 0; i < network->numNodes; i++) {
        bushes->LPcost[i] = -INFINITY;
        bushes->SPcost[i] = INFINITY;
    }

    /* Ensure costs are up to date */
    if (parameters->linkShiftB != &exactCostUpdate) updateAllCosts(network);

    bushes->SPcost[origin2node(network, origin)] = 0;
    bushes->LPcost[origin2node(network, origin)] = 0;
    for (curnode = 1; curnode < network->numNodes; curnode++) {
        i = bushes->bushOrder[origin][curnode];
        /* Iterate over incoming links */
        //displayMessage(DEBUG, "Scanning node %d:%d (%d)\n", origin, i, curnode);
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            merge = bushes->merges[origin][m];
            /* Find shortest incoming link */
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->SPcost[h] + network->arcs[hi].cost;
                if (tempcost < bushes->SPcost[i]) {
                    bushes->SPcost[i] = tempcost;
                    merge->SPlink = curarc;
                }
            }
            /* Find longest (need to separate out depending on LPrule) */
            if (LPrule == NO_LONGEST_PATH) continue;
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->LPcost[h] + network->arcs[hi].cost;     
                if (tempcost > bushes->LPcost[i]
                  && (LPrule == LONGEST_BUSH_PATH
                    || (LPrule == LONGEST_USED_PATH
                      && merge->approachFlow[curarc] > parameters->minLinkFlow)
                    || (LPrule == LONGEST_USED_OR_SP
                      &&(merge->approachFlow[curarc] > parameters->minLinkFlow
                      || merge->SPlink == curarc)))) {
                    
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
    scanBushes(origin, network, bushes, parameters, LONGEST_USED_OR_SP);
    calculateBushFlows(origin, network, bushes);
  
    /* Make a first pass... */
    for (ij = 0; ij < network->numArcs; ij++) {
        /* Mark links with near-zero contribution for removal by setting to
         * exactly zero */
        if (bushes->flow[ij] < parameters->minLinkFlow) {
            if (isInBush(origin, ij, network, bushes) == TRUE)
                displayMessage(FULL_DEBUG, "Attempting to delete (%d,%d)\n", 
                           network->arcs[ij].tail+1, network->arcs[ij].head+1);
            bushes->flow[ij] = 0;
        }
        if (bushes->flow[ij] > 0) continue; /* Link is already in the bush, no
                                               need to add it again */
        /* See if the link provides a shortcut using the strict criterion */
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        if (bushes->LPcost[i] == -INFINITY && bushes->LPcost[j] > -INFINITY)
            continue; /* No path to extend */
        if (bushes->SPcost[i] + network->arcs[ij].cost < bushes->SPcost[j]
            && bushes->LPcost[i] < bushes->LPcost[j]
            && (network->arcs[ij].tail == origin2node(network, origin)
                || network->arcs[ij].tail >= network->firstThroughNode)) 
        {
            bushes->flow[ij] = NEW_LINK;
            newArcs++;
        /* Never delete shortest path tree... should be OK with floating point
         * comparison since this is how SPcost is calculated */
        } else if (bushes->SPcost[i]+network->arcs[ij].cost==bushes->SPcost[j] 
                   && bushes->flow[ij] == 0
                   && isInBush(origin, ij, network, bushes) == TRUE) {
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
                && (network->arcs[ij].tail == origin2node(network, origin)
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
        if (i == origin2node(network, origin)) continue;
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
            if (SPnode != origin2node(network, origin)) {
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
            if (LPnode != origin2node(network, origin)) {
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
 * - rescanAndCheck both updates longtest/shortest path labels, but also checks
 *   whether these are different enough to actually do anything -- depending
 *   on parameters we can skip bushes which are closer to equilibrium, to 
 *   spend more time on bushes further away.
 * - updateFlowPass then does the actual flow shifting, inside of a loop in
 *   case we want to do multiple shifts per origin.
 */
bool updateFlowsB(int origin, network_type *network, bushes_type *bushes,
                  algorithmBParameters_type *parameters) {
    int i;
 
    /* Recompute bush flows for this origin */
    calculateBushFlows(origin, network, bushes);

    /* Update longest/shortest paths, check whether there is work to do */
    if (rescanAndCheck(origin, network, bushes, parameters) == FALSE) {
        bushes->updateBush[origin] = FALSE;
        displayMessage(DEBUG, "bailing out\n");
        return FALSE;
    }

    /* Now do (possibly multiple) shifts per origin */
    for (i = 0; i < parameters->shiftReps; i++) {
        updateFlowPass(origin, network, bushes, parameters);
        parameters->numFlowShifts++;
        /* Uncomment next line for extra validation checking */
        /* checkFlows(network, bushes); */
        if (parameters->rescanAfterShift == TRUE
                && i + 1 < parameters->shiftReps) {
            if (rescanAndCheck(origin, network, bushes, parameters) == FALSE)
                return TRUE;
        }
    }
    
    return TRUE;
}

/*
 * rescanAndCheck -- calls functions updating the longest/shortest path labels
 * on the bush, aint with divergence nodes; and then checks whether the bush
 * is close enough to equilibrium to skip for now (in which case the function
 * returns FALSE).  If we need to shift flows, it returns TRUE.
 */
bool rescanAndCheck(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    int i, j, ij;
    double maxgap = 0;
    double bushSPTT = 0, bushExcess = 0;

    scanBushes(origin, network, bushes, parameters, LONGEST_USED_PATH);
    findDivergenceNodes(origin, network, bushes);
    for (i = 0; i < network->numZones; i++) {
        bushSPTT += network->demand[origin][i] * bushes->SPcost[i];
        maxgap = max(maxgap, fabs(bushes->LPcost[i] -bushes->SPcost[i]));
    }
    for (; i < network->numNodes; i++) {
        maxgap = max(maxgap, fabs(bushes->LPcost[i] -bushes->SPcost[i]));
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        bushExcess += bushes->flow[ij] * (network->arcs[ij].cost +
                        bushes->SPcost[i] - bushes->SPcost[j]);
    }
    displayMessage(DEBUG, "Scanning %d, gap is %f\n", origin, 
                    bushExcess / bushSPTT);
    displayMessage(DEBUG, "Max gap, threshold, threshold AEC: %f %f %f\n",
                    maxgap, parameters->thresholdGap, parameters->thresholdAEC);
    if (maxgap < parameters->thresholdGap) return FALSE;
    if (bushExcess / bushSPTT < parameters->thresholdAEC) return FALSE;

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
        bushes->nodeFlow[i] = network->demand[origin][i];

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
                || (fabs(bushes->LPcost[j] - bushes->SPcost[j])
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
    for (i = 0; i < network->numZones; i++)
        bushes->nodeFlow[i] = network->demand[origin][i];
    for (; i < network->numNodes; i++)
        bushes->nodeFlow[i] = 0;
    for (ij = 0; ij < network->numArcs; ij++)
        bushes->flow[ij] = 0;

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
 * rectifyBushFlows -- identical to calculateBushFlows but rectifying each
 * merge along the way.  Especially useful when warm-starting or initializing
 * one class with another's OD matrix... it will keep the same proportions
 * but with new demand.
 */
void rectifyBushFlows(int origin,network_type *network,bushes_type *bushes) {
    int i, j, ij, m, node;
    merge_type *merge;

    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++)
        bushes->nodeFlow[i] = network->demand[origin][i];
    for (; i < network->numNodes; i++)
        bushes->nodeFlow[i] = 0;
    for (ij = 0; ij < network->numArcs; ij++)
        bushes->flow[ij] = 0;

    /* Descending pass for flow calculations  */
    for (node = network->numNodes - 1; node > 0; node--) {
        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple(j, origin, network, bushes);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            rectifyMerge(j, merge, bushes);
            pushBackFlowMerge(merge, network, bushes);
        }
    }

}


/*
 * newtonFlowShift -- For a given merge node, use Newton's method to shift flow
 * from the inter segment to the shorter one.
 */
void newtonFlowShift(int j, merge_type *merge, int origin,
                     network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters) {
    double flow1, flow2, cost1, cost2, der1, der2, shift;
    int i, hi, m, c;
    merge_type *segmentMerge;

    c = origin2class(network, origin);

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
        network->arcs[hi].classFlow[c] -= shift;
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
        network->arcs[hi].classFlow[c] += shift;
        parameters->linkShiftB(hi, shift, network);
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
   
    for (r = 0; r < network->batchSize; r++) {
        if (outOfOrigins(network, r) == TRUE) break;
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
            if (i == origin2node(network, r)) continue;
            balance = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                    curArc = curArc->next) {
                balance += bushes->flow[ptr2arc(network, curArc->arc)];
            }
            for (curArc = network->nodes[i].forwardStar.head; curArc != NULL;
                    curArc = curArc->next) {
                balance -= bushes->flow[ptr2arc(network, curArc->arc)];
            }         
            if (i < network->numZones) balance -= network->demand[r][i];
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
            fatalError("Flow validation failed: link (%ld,%ld) has aggregate "
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
    network->arcs[ij].flow += shift;
    network->arcs[ij].cost=network->arcs[ij].calculateCost(&network->arcs[ij]);
    network->arcs[ij].der = network->arcs[ij].calculateDer(&network->arcs[ij]);
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


void writeBushes(network_type *network, bushes_type *bushes, char *filename) {
    int m, origin, check;
    merge_type *merge;
    FILE *batchFile = openFile(filename, "wb");
        
    /* Write batch size as header */
    check = network->batchSize;
    fwrite(&check, sizeof(check), 1, batchFile);
    
    /* Now write each origin in turn */
    for (origin = 0; origin < network->batchSize; origin++) {
        if (outOfOrigins(network, origin) == TRUE) break;

        /* First write the origin ID as a check */
        check = origin + network->batchSize * network->curBatch;
        fwrite(&check, sizeof(check), 1, batchFile);
        
        /* Now write last merge and number of merges*/
        fwrite(&bushes->lastMerge[origin], sizeof(bushes->lastMerge[0]),
               1, batchFile);
        fwrite(&bushes->numMerges[origin], sizeof(bushes->numMerges[0]),
               1, batchFile);

        /* Now write the topological order */
        fwrite(bushes->bushOrder[origin], sizeof(bushes->bushOrder[origin][0]),
               network->numNodes, batchFile);
        
        /* Now write the predecessors */
        fwrite(bushes->pred[origin], sizeof(bushes->pred[origin][0]),
               network->numNodes, batchFile);
               
        /* Finally write the merge information */
        for (m = 0; m < bushes->numMerges[origin]; m++) {
            merge = bushes->merges[origin][m];        
            fwrite(&merge->numApproaches, sizeof(merge->numApproaches),
                   1, batchFile);
            fwrite(merge->approach, sizeof(merge->approach[0]),
                   merge->numApproaches, batchFile);
            fwrite(merge->approachFlow, sizeof(merge->approachFlow[0]),
                   merge->numApproaches, batchFile);                   
            fwrite(&merge->SPlink, sizeof(merge->SPlink), 1, batchFile);
            fwrite(&merge->LPlink, sizeof(merge->LPlink), 1, batchFile);            
            fwrite(&merge->divergenceNode, sizeof(merge->divergenceNode),
                   1, batchFile);  
        }
    }
    
    fclose(batchFile);
}

void readBushes(network_type *network, bushes_type **bushes, char *filename) {
    int m, origin, check, numApproaches;
    merge_type *merge;
    FILE *batchFile = openFile(filename, "rb");
    
    /* First do a check */
    my_fread(&check, sizeof(check), 1, batchFile);
    if (check != network->batchSize) {
        fatalError("File %s has the wrong batch size for this network.",
                   filename);
    }
    
    /* First clear out all existing bush data */
    deleteBushes(network, *bushes);
    *bushes = createBushes(network);
    
    /* Now read each origin in turn */
    for (origin = 0; origin < network->batchSize; origin++) {
        if (outOfOrigins(network, origin) == TRUE) break;
        /* First read the origin ID and check */
        my_fread(&check, sizeof(check), 1, batchFile);
        if (check != origin + network->batchSize * network->curBatch) {
            fatalError("Reading wrong origin from file %s.", filename);
        }        
        
        /* Now read last merge and number of merges*/
        my_fread(&(*bushes)->lastMerge[origin],sizeof((*bushes)->lastMerge[0]),
               1, batchFile);
        my_fread(&(*bushes)->numMerges[origin],sizeof((*bushes)->numMerges[0]),
               1, batchFile);

        /* Now read the topological order */
        my_fread((*bushes)->bushOrder[origin],
               sizeof((*bushes)->bushOrder[origin][0]),
               network->numNodes, batchFile);
        
        /* Now read the predecessors */
        my_fread((*bushes)->pred[origin], sizeof((*bushes)->pred[origin][0]),
               network->numNodes, batchFile);
               
        /* Finally read the merge information */
        deleteVector((*bushes)->merges[origin]);
        (*bushes)->merges[origin] = newVector((*bushes)->numMerges[origin],
                                           merge_type *);
        for (m = 0; m < (*bushes)->numMerges[origin]; m++) {
            my_fread(&numApproaches, sizeof(numApproaches), 1, batchFile);
            merge = createMerge(numApproaches); 
            my_fread(merge->approach, sizeof(merge->approach[0]),
                   merge->numApproaches, batchFile);
            my_fread(merge->approachFlow, sizeof(merge->approachFlow[0]),
                   merge->numApproaches, batchFile);                   
            my_fread(&merge->SPlink, sizeof(merge->SPlink), 1, batchFile);
            my_fread(&merge->LPlink, sizeof(merge->LPlink), 1, batchFile);
            my_fread(&merge->divergenceNode, sizeof(int),
                   1, batchFile);  
            (*bushes)->merges[origin][m] = merge;        
        }
    }
            
    fclose(batchFile);
}

/*
 * Reads a binary OD matrix into the network, for network->curBatch
 */
void readBinaryMatrix(network_type *network,
                      algorithmBParameters_type *parameters) {
    int check, r;
    char filename[STRING_SIZE];
    FILE* matrixFile;

    sprintf(filename, "%s%d.bin", parameters->matrixStem, network->curBatch);
    matrixFile = openFile(filename, "rb");
    my_fread(&check, sizeof(check), 1, matrixFile);
    if (check != network->curBatch) fatalError("Reading wrong binary matrix.");
//    displayMessage(FULL_NOTIFICATIONS, "Reading binary matrix %d\n", check);
    for (r = 0; r < network->batchSize; r++) {
        my_fread(network->demand[r], sizeof(network->demand[r][0]),
                  network->numZones, matrixFile);
    }
//    displayMessage(FULL_NOTIFICATIONS, "finished reading binary matrix %d\n", check);
    fclose(matrixFile);
}
