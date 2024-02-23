#include "fileio.h"

///////////////////////////
// Reading network files //
///////////////////////////

network_type *readParametersFile(algorithmBParameters_type *thisRun,
                                 char *filename) {
    /* Read AlgorithmB parameters from a file.  Also read and load
     * network data from the indicated files, and return the pointer
     * to the network data structure created in this way. */
	int status;
    int c;
	char fullLine[STRING_SIZE * 2]; // Double-length to allow concatenation
    char filePath[STRING_SIZE], dataPath[STRING_SIZE];
    char networkFileName[STRING_SIZE];
    stringList_type *tripsFile = NULL, *newTripsFile = NULL;
	char metadataTag[STRING_SIZE], metadataValue[STRING_SIZE];
	FILE *parametersFile = openFile(filename, "r");
    int numBatches = 1;
    network_type *network = newScalar(network_type);
#ifdef PARALLELISM
    thisRun->numThreads = getNumCores();
#endif

	/* Initialize (set mandatory values to missing, mandatory strings to length
	   zero, others to defaults).  Note that some parameters are initialized
       by initializeAlgorithmBParameters, assuming that this was called for
       thisRun when it was created.  These default values are not repeated
       here. */
    networkFileName[0] = '\0';
    strcpy(filePath, "net/");
    strcpy(dataPath, "net/");
    strcpy(thisRun->flowsFile, "flows.txt");
    network->numClasses = 0;

	/* Process parameter file */
	while (!feof(parametersFile)) {
		do {
			if (fgets(fullLine, STRING_SIZE, parametersFile) == NULL) break;
			status = parseMetadata(fullLine, metadataTag, metadataValue);
		} while (status == BLANK_LINE || status == COMMENT);
		if        (strcmp(metadataTag, "NETWORK FILE") == 0) {
			strcpy(networkFileName, metadataValue);
		} else if (strcmp(metadataTag, "TRIPS FILE") == 0) {
			newTripsFile = newScalar(stringList_type);
            strcpy(newTripsFile->string, metadataValue);
            newTripsFile->prev = tripsFile;
            tripsFile = newTripsFile; 
            network->numClasses++;
		} else if (strcmp(metadataTag, "CONVERGENCE GAP") == 0) {
            thisRun->convergenceGap = atof(metadataValue);
		} else if (strcmp(metadataTag, "MAX ITERATIONS") == 0) {
            thisRun->maxIterations = atoi(metadataValue);
		} else if (strcmp(metadataTag, "MAX RUN TIME") == 0) {
            thisRun->maxTime = atof(metadataValue);
		} else if (strcmp(metadataTag, "FILE PATH") == 0) {
            strcpy(filePath, metadataValue);
		} else if (strcmp(metadataTag, "FLOWS FILE") == 0) {
            strcpy(thisRun->flowsFile, metadataValue);
		} else if (strcmp(metadataTag, "GAP FUNCTION") == 0) {
			if    (strcmp(metadataValue, "RELATIVE GAP") == 0)
				thisRun->gapFunction = RELATIVE_GAP_1;
			else if (strcmp(metadataValue, "RELATIVE GAP 2") == 0)
				thisRun->gapFunction = RELATIVE_GAP_2;
			else if (strcmp(metadataValue, "AVERAGE EXCESS COST") == 0)
				thisRun->gapFunction = AEC;
			else if (strcmp(metadataValue, "MAXIMUM EXCESS COST") == 0)
				thisRun->gapFunction = MEC;
			else
				fatalError("Unknown gap function %s\n", metadataValue);
		} else if (strcmp(metadataTag, "OMIT BECKMANN") == 0) {
            thisRun->calculateBeckmann = FALSE;
		} else if (strcmp(metadataTag, "DEMAND MULTIPLIER") == 0) {
            thisRun->demandMultiplier = atof(metadataValue);
		} else if (strcmp(metadataTag, "STORE MATRICES") == 0) {
            thisRun->storeMatrices = TRUE;
		} else if (strcmp(metadataTag, "DATA PATH") == 0) {
            strcpy(dataPath, metadataValue);
		} else if (strcmp(metadataTag, "NUMBER OF BATCHES") == 0) {
            numBatches = atoi(metadataValue);
		} else if (strcmp(metadataTag, "BATCH STEM") == 0) {
            strncpy(thisRun->batchStem, metadataValue,
                    sizeof(thisRun->batchStem) - 1);
		} else if (strcmp(metadataTag, "MATRIX STEM") == 0) {
            strncpy(thisRun->matrixStem, metadataValue, 
                    sizeof(thisRun->matrixStem) - 1);
		} else if (strcmp(metadataTag, "WARM START") == 0) {
            thisRun->warmStart = TRUE;
		} else if (strcmp(metadataTag, "NUMBER OF THREADS") == 0) {
#ifdef PARALLELISM
            thisRun->numThreads = atoi(metadataValue);
#else
            warning(LOW_NOTIFICATIONS, "This is the serial build of tap-b, "
                                       "ignoring specified number of threads.");
#endif
		} else if (strcmp(metadataTag, "INNER ITERATIONS") == 0) {
            thisRun->innerIterations = atoi(metadataValue);
		} else if (strcmp(metadataTag, "SHIFT REPETITIONS") == 0) {
            thisRun->shiftReps = atoi(metadataValue);
		} else if (strcmp(metadataTag, "RESCAN AFTER SHIFT") == 0) {
            thisRun->rescanAfterShift = TRUE;
		} else if (strcmp(metadataTag, "THRESHOLD GAP") == 0) {
            thisRun->thresholdGap = atof(metadataValue);
		} else if (strcmp(metadataTag, "THRESHOLD AEC") == 0) {
            thisRun->thresholdAEC = atof(metadataValue);
		} else if (strcmp(metadataTag, "MIN COST DIFFERENCE") == 0) {
            thisRun->minCostDifference = atof(metadataValue);
		} else if (strcmp(metadataTag, "MIN LINK FLOW SHIFT") == 0) {
            thisRun->minLinkFlowShift = atof(metadataValue);
		} else if (strcmp(metadataTag, "MIN LINK FLOW") == 0) {
            thisRun->minLinkFlow = atof(metadataValue);
		} else if (strcmp(metadataTag, "MIN DERIVATIVE") == 0) {
            thisRun->minDerivative = atof(metadataValue);
		} else if (strcmp(metadataTag, "NEWTON STEP") == 0) {
            thisRun->newtonStep = atof(metadataValue);
		} else if (strcmp(metadataTag, "NEWTON SHIFTS") == 0) {
            thisRun->numNewtonShifts = atoi(metadataValue);
		} else if (strcmp(metadataTag, "SHORTEST PATH QUEUE DISCIPLINE") == 0) {
			if    (strcmp(metadataValue, "DEQUE") == 0)
				thisRun->SPQueueDiscipline = DEQUE;
			else if (strcmp(metadataValue, "FIFO") == 0)
				thisRun->SPQueueDiscipline = FIFO;
			else if (strcmp(metadataValue, "LIFO") == 0)
				thisRun->SPQueueDiscipline = LIFO;
			else
				fatalError("Unknown queue discipline %s\n", metadataValue);
		} else if (strcmp(metadataTag, "LINK COST UPDATE") == 0) {
			if    (strcmp(metadataValue, "EXACT") == 0)
				thisRun->linkShiftB = &exactCostUpdate;
			else if (strcmp(metadataValue, "LINEAR") == 0)
				thisRun->linkShiftB = &linearCostUpdate;
			else if (strcmp(metadataValue, "NONE") == 0)
				thisRun->linkShiftB = &noCostUpdate;
			else
				fatalError("Unknown link cost update %s\n", metadataValue);
		} else if (strcmp(metadataTag, "STORE BUSHES") == 0) {
            thisRun->storeBushes = TRUE;
		} else if (strcmp(metadataTag, "REUSE FIRST BUSH") == 0) {
            thisRun->reuseFirstBush = TRUE;
		} else if (strcmp(metadataTag, "EXCLUDE GAP TIME") == 0) {
            thisRun->includeGapTime = FALSE;
		} else if (strcmp(metadataTag, "INITIAL BUSH") == 0) {
			if    (strcmp(metadataValue, "SP_TREE") == 0)
				thisRun->createInitialBush = &initialBushShortestPath;
			else
				fatalError("Unknown initial bush method %s\n", metadataValue);
		} else if (strcmp(metadataTag, "TOPOLOGICAL ORDER") == 0) {
			if    (strcmp(metadataValue, "STANDARD") == 0)
				thisRun->topologicalOrder = &genericTopologicalOrder;
			else
				fatalError("Unknown topological sort %s\n", metadataValue);
		} else {
			warning(MEDIUM_NOTIFICATIONS,
			        "Ignoring unknown metadata tag in parameters file - %s\n",
			        metadataTag);
		}
    }
    fclose(parametersFile);

    /* Check mandatory elements are present and validate input */
	if (strlen(networkFileName) == 0)
		fatalError("Missing network file!");
	if (network->numClasses <= 0)
		fatalError("Missing demand file!");
	if (thisRun->maxIterations == INT_MAX
	        && thisRun->maxTime == INFINITY
	        && thisRun->convergenceGap == 0)
		warning(LOW_NOTIFICATIONS,
		        "No termination criteria specified... program will run until "
		        "interrupted manually.\n");
    if (thisRun->demandMultiplier < 0)
        fatalError("Negative demand multiplier.");
    if (numBatches < 1)
        fatalError("Must use at least one batch.");
    if (thisRun->numThreads < 1)
        fatalError("Must use at least one thread.");
    if (thisRun->innerIterations < 0)
        warning(LOW_NOTIFICATIONS, "Negative number of inner iterations.");
    if (thisRun->shiftReps < 0)
        warning(LOW_NOTIFICATIONS, "Negative number of shift repetitions.");

    /* Concatenate paths, read network */
    int saveDigits = ceil(log10(network->numClasses)) + 4; // #classes + ".bin"
    snprintf(fullLine, 2*STRING_SIZE, "%s%s", filePath, networkFileName);
    strncpy(networkFileName, fullLine, STRING_SIZE - 1);
    snprintf(fullLine, 2*STRING_SIZE - saveDigits, "%s%s",
            dataPath, thisRun->batchStem);
    strncpy(thisRun->batchStem, fullLine, STRING_SIZE - 1);
    snprintf(fullLine, 2*STRING_SIZE - saveDigits, "%s%s",
            dataPath, thisRun->matrixStem);
    strncpy(thisRun->matrixStem, fullLine, STRING_SIZE - 1);
    char **tripsFiles = (char **)malloc(network->numClasses * sizeof(char *));
    if (tripsFiles == NULL) fatalError("Can't set up trip file array.\n");
    for (c = network->numClasses - 1; c >= 0; c--) {
        tripsFiles[c] = (char *)malloc(2 * STRING_SIZE * sizeof(char));
        if (tripsFiles[c] == NULL) fatalError("Can't setup trip file array.\n");
        snprintf(tripsFiles[c], 2*STRING_SIZE - 1, "%s%s",
                 filePath, tripsFile->string);
        newTripsFile = tripsFile; /* Save for garbage collection */
        tripsFile = tripsFile->prev;
        deleteScalar(newTripsFile);
    }
    
    /* garbage collection */
    readOBANetwork(network, networkFileName, tripsFiles, network->numClasses,
                   thisRun->demandMultiplier);
    if (thisRun->storeMatrices == TRUE) {
        writeBinaryMatrices(network, thisRun->matrixStem);
    }
    for (c = network->numClasses - 1; c >= 0; c--) {
        free(tripsFiles[c]);
    }
    free(tripsFiles);
    setBatches(network, network->numOrigins / numBatches, thisRun->warmStart,
               thisRun->storeMatrices, thisRun->matrixStem);
    return network;
}

/*
 * Write OD matrices into binary files, one for each batch.
 * For use in later warm-starts or assignments.
 */
void writeBinaryMatrices(network_type *network, char *matrixStem) {
    int batch, r, origin;
    char filename[STRING_SIZE];
    FILE* matrixFile;

    for (batch = 0; batch < network->numBatches; batch++) {
        sprintf(filename, "%s%d.bin", matrixStem, batch);
        matrixFile = openFile(filename, "wb");
        fwrite(&batch, sizeof(batch), 1, matrixFile); /* Header */
        for (r = 0; r < network->batchSize; r++) {
            if (outOfOrigins(network, r) == TRUE) break;
            origin = r + batch * network->batchSize;
            fwrite(network->demand[origin], sizeof(network->demand[origin][0]),
                   network->numZones, matrixFile);
        }
        fclose(matrixFile);
    }
}

void readOBANetwork(network_type *network, char *linkFileName, 
                    char **tripFileName, int numClasses,
                    double defaultDemandMultiplier) {
    int i, j, r = 0, c;
    int check;
    int numParams, status;
    double demand, totalDemandCheck = 0, statedDemand = 0;
    double defaultTollFactor, defaultDistanceFactor;
    double demandMultiplier = IS_MISSING;

    char fullLine[STRING_SIZE], trimmedLine[STRING_SIZE], *token;
    char metadataTag[STRING_SIZE], metadataValue[STRING_SIZE];

    FILE *linkFile = openFile(linkFileName, "r");
    FILE *tripFile;

    network->numZones = IS_MISSING;
    network->numArcs = IS_MISSING;
    network->numNodes = IS_MISSING;
    network->numClasses = numClasses;
    network->firstThroughNode = IS_MISSING;
    defaultTollFactor = IS_MISSING;
    defaultDistanceFactor = IS_MISSING;

    /* Read link file metadata */
    bool endofMetadata = FALSE;
    do {
        if (fgets(fullLine, STRING_SIZE, linkFile) == NULL) 
            fatalError("Link file %s ended (or other I/O error) before "
                        "metadata complete.", linkFileName);
        status = parseMetadata(fullLine, metadataTag, metadataValue);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if         (strcmp(metadataTag, "NUMBER OF ZONES") == 0) {
            network->numZones = atoi(metadataValue);
        } else if (strcmp(metadataTag, "NUMBER OF LINKS") == 0) {
            network->numArcs = atoi(metadataValue);
        } else if (strcmp(metadataTag, "NUMBER OF NODES") == 0) {
            network->numNodes = atoi(metadataValue);
        } else if (strcmp(metadataTag, "FIRST THRU NODE") == 0) {
            network->firstThroughNode = atoi(metadataValue) - 1;
        } else if (strcmp(metadataTag, "DISTANCE FACTOR") == 0) {
            defaultDistanceFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "TOLL FACTOR") == 0) {
            defaultTollFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "END OF METADATA") == 0) {
            endofMetadata = TRUE;
        } else {
            warning(MEDIUM_NOTIFICATIONS, "Ignoring unknown metadata tag %s "
                    "in link file %s", metadataTag, linkFileName);
        }
    } while (endofMetadata == FALSE);

    /* Check input for completeness and correctness */
    if (network->numZones == IS_MISSING) 
        fatalError("Link file %s does not contain number of zones.", 
                linkFileName);
    if (network->numNodes == IS_MISSING) 
        fatalError("Link file %s does not contain number of nodes.", 
                linkFileName);
    if (network->numArcs == IS_MISSING) 
        fatalError("Link file %s does not contain number of links.", 
                linkFileName);
    if (network->firstThroughNode == IS_MISSING) {
        warning(LOW_NOTIFICATIONS, "Link file %s does not contain first "
                "through node, setting to 1 as default.\n", linkFileName);
        network->firstThroughNode = 0;
    }
    if (defaultDistanceFactor == IS_MISSING) {
        defaultDistanceFactor = 0;
    }
    if (defaultTollFactor == IS_MISSING) {
        defaultTollFactor = 0;
    }
    if (network->numZones < 1) 
        fatalError("Link file %s does not contain a positive number of nodes.",
                linkFileName);
    if (network->numArcs < 1) 
        fatalError("Link file %s does not contain a positive number of links.",
                linkFileName);
    if (network->numNodes < 1) 
        fatalError("Link file %s does not contain a positive number of nodes.",
                linkFileName);

    displayMessage(MEDIUM_NOTIFICATIONS, "Global distance/toll factors: "
            "%lf %lf\n", defaultDistanceFactor, defaultTollFactor);

    network->numOrigins = network->numZones * network->numClasses;
    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->demand = newMatrix(network->numOrigins, network->numZones,double);
    network->tollFactor = newVector(network->numClasses, double);
    network->distanceFactor = newVector(network->numClasses, double);
#ifdef PARALLELISM
    network->arc_muts = newVector(network->numArcs, pthread_mutex_t);
#endif

    /* Default batching for network reading; can adjust later */
    network->batchSize = network->numOrigins;
    network->numBatches = 1;
    network->curBatch = 0;
    
    for (i = 0; i < network->numOrigins; i++) {
        for (j = 0; j < network->numZones; j++) {
            network->demand[i][j] = 0;
        }
    }

    /* Read link data */
    for (i = 0; i < network->numArcs; i++) {
        if (fgets(fullLine, STRING_SIZE, linkFile) == NULL)
            fatalError("Link file %s ended (or other I/O error) before link "
                    "data complete.", linkFileName);
        status = parseLine(fullLine, trimmedLine);
        if (status == BLANK_LINE || status == COMMENT) {
            i--;
            continue;
        }
        numParams=sscanf(trimmedLine,"%d %d %lf %lf %lf %lf %lf %lf %lf %d",
            &network->arcs[i].tail,
            &network->arcs[i].head,
            &network->arcs[i].capacity,
            &network->arcs[i].length,
            &network->arcs[i].freeFlowTime,
            &network->arcs[i].alpha,
            &network->arcs[i].beta,
            &network->arcs[i].speedLimit,
            &network->arcs[i].toll,
            &network->arcs[i].linkType);
        if (numParams != 10) 
            fatalError("Link file %s has an error in this line:\n\"%s\"",
                    linkFileName, fullLine);
        if (network->arcs[i].tail < 1 
                || network->arcs[i].tail > network->numNodes) 
            fatalError("Arc tail %d out of range in network file %s.", 
                    i, linkFileName);
        if (network->arcs[i].head < 1 
                || network->arcs[i].head > network->numNodes) 
            fatalError("Arc head %d out of range in network file %s.", 
                    i, linkFileName);
        if (network->arcs[i].length < 0) 
            warning(FULL_NOTIFICATIONS, 
                    "Arc length %d negative in network file %s.\n%s", i,
                    linkFileName, fullLine);
        if (network->arcs[i].freeFlowTime < 0) 
            fatalError("Arc free flow time %d negative in network file "
                    "%s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].alpha < 0) fatalError("Alpha %d negative in "
                "network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].beta < 0) fatalError("Beta %d negative in "
                "network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].speedLimit < 0) warning(FULL_NOTIFICATIONS, 
                "Speed limit %d negative in network file %s.\n%s", i, 
                linkFileName, fullLine);
        if (network->arcs[i].toll < 0) warning(FULL_NOTIFICATIONS, "Toll %d "
                "negative in network file %s.\n%s", i, linkFileName, fullLine);
        if (network->arcs[i].capacity <= 0) fatalError("Capacity %d "
                "nonpositive in network file %s.\n%s", i, linkFileName, 
                fullLine);
        network->arcs[i].tail--;
        network->arcs[i].head--;
        network->arcs[i].flow = 0;
        network->arcs[i].cost = network->arcs[i].freeFlowTime;
        if (network->arcs[i].beta == 1) {
           network->arcs[i].calculateCost = &linearBPRcost;
           network->arcs[i].calculateDer = &linearBPRder;           
           network->arcs[i].calculateInt = &linearBPRint;           
        } else if (network->arcs[i].beta == 4) {
           network->arcs[i].calculateCost = &quarticBPRcost;
           network->arcs[i].calculateDer = &quarticBPRder;           
           network->arcs[i].calculateInt = &quarticBPRint;           
        } else {
           network->arcs[i].calculateCost = &generalBPRcost;
           network->arcs[i].calculateDer = &generalBPRder;           
           network->arcs[i].calculateInt = &generalBPRint;           
        }           
        network->arcs[i].classFlow = newVector(network->numClasses, double);
        network->arcs[i].classCost = newVector(network->numClasses, double);
        network->arcs[i].classToll = newVector(network->numClasses, double);        
        for (c = 0; c < network->numClasses; c++) {
            /* Assume all costs face same toll unless otherwise stated.
               (e.g., to make a toll-insensitive class you can make its toll
               zero) */
            network->arcs[i].classToll[c] = network->arcs[i].toll;
        }
        /* classFlow and classCost are initialized in finalizeNetwork */
    }

    for (c = 0; c < network->numClasses; c++) {
        totalDemandCheck = 0; statedDemand = 0;
        tripFile = openFile(tripFileName[c], "r");
        /* Verify trip table metadata */
        endofMetadata = FALSE;
        network->totalODFlow = 0;
        network->tollFactor[c] = defaultTollFactor;
        network->distanceFactor[c] = defaultDistanceFactor;
        do {
            if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) 
                fatalError("Trip file %s ended (or other I/O error) before "
                        "metadata complete.", tripFileName);
            status = parseMetadata(fullLine, metadataTag, metadataValue);
            if (status == BLANK_LINE || status == COMMENT) continue;
            if         (strcmp(metadataTag, "NUMBER OF ZONES") == 0) {
                check = atoi(metadataValue);
                if (check != network->numZones) fatalError("Number of zones in"
                        "trip and link files do not match.");
            } else if (strcmp(metadataTag, "TOTAL OD FLOW") == 0) {
                statedDemand = atof(metadataValue);
            } else if (strcmp(metadataTag, "DEMAND MULTIPLIER") == 0) {
                demandMultiplier = atof(metadataValue);
            } else if (strcmp(metadataTag, "DISTANCE FACTOR") == 0) {
                network->distanceFactor[c] = atof(metadataValue);
            } else if (strcmp(metadataTag, "TOLL FACTOR") == 0) {
                network->tollFactor[c] = atof(metadataValue);
            } else if (strcmp(metadataTag, "END OF METADATA") == 0) {
                endofMetadata = TRUE;
            } else {
                warning(MEDIUM_NOTIFICATIONS, "Ignoring unknown metadata tag "
                        "%s in trips file %s", metadataTag, tripFileName[c]);
            }
        } while (endofMetadata == FALSE);
        if (demandMultiplier == IS_MISSING)
            demandMultiplier = defaultDemandMultiplier;

        /* Now read trip table */
        while (!feof(tripFile)) {
            if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) break;
            status = parseLine(fullLine, trimmedLine);
            if (status == BLANK_LINE || status == COMMENT) continue;
            if (strstr(trimmedLine, "Origin") != NULL) {
                // i indexes current origin
                sscanf(strstr(trimmedLine, "Origin")+6,"%d", &i);  
                if (i <= 0 || i > network->numZones) fatalError("Origin %d is"
                        "out of range in trips file %s", i, tripFileName[c]);
                i--;
                r = nodeclass2origin(network, i, c);
                continue;
            }
            token = strtok(trimmedLine , ";");
            while (token != NULL && strlen(token) > 1) {
                numParams = sscanf(token, "%d : %lf", &j, &demand);
                if (numParams < 2) break;
                if (j <= 0 || j > network->numZones) fatalError("Destination "
                        "%d is out of range in trips file %s\n%s\n%s", j, 
                        tripFileName[c], fullLine, token);
                j--;
                network->demand[r][j] = demand * demandMultiplier;
                if (demand < 0) fatalError("Negative demand from origin %d to "
                        "destination %d in class %d", i, j, c);
                totalDemandCheck += demand;
                network->totalODFlow += network->demand[r][j];
                token = strtok(NULL, ";");
            }
            blankInputString(trimmedLine, STRING_SIZE);
        }
        displayMessage(MEDIUM_NOTIFICATIONS, "Class %d distance/toll factors: "
                                             "%f %f\n", c+1,
                                             network->distanceFactor[c],
                                             network->tollFactor[c]);
        displayMessage(MEDIUM_NOTIFICATIONS,
                       "Read %f trips for class %d (%f expected)\n",
                       totalDemandCheck, c + 1, statedDemand);
        if (fabs(totalDemandCheck / statedDemand - 1) > 0.01) {
            warning(LOW_NOTIFICATIONS, "Class %d demand differs from stated "
                                       "value by more than 1%\n", c+1);
        }
        fclose(tripFile);
    }

    displayMessage(MEDIUM_NOTIFICATIONS, "%d nodes, %d arcs, and %d zones\n",
            network->numNodes, network->numArcs, network->numZones);

    fclose(linkFile);
    displayMessage(LOW_NOTIFICATIONS, "File read and memory allocated.\n");

    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS, "Forward and reverse star lists "
            "generated.\n");

}

/*
setBatches: Re-partitions the network into origin batches of the given
size.  This should NEVER be called during the middle of a run, or
unpredictable things may happen.

As a result, we write the binary matrices HERE.
*/
void setBatches(network_type *network, int batchSize, bool warmStart,
                bool storeMatrices, char *matrixStem) {
    network->batchSize = batchSize;
    if (network->numOrigins % batchSize != 0) {
        fatalError("Number of origins (%d) must be divisible by the batch "
                   "size(%d)\n", network->numOrigins, batchSize);
    }
    network->numBatches = (network->numOrigins - 1) / batchSize + 1;
    network->curBatch = 0;

    if (warmStart == FALSE
            && (network->numBatches > 1 || storeMatrices == TRUE ))
        writeBinaryMatrices(network, matrixStem);

    if (network-> numBatches > 1 && warmStart == FALSE) {
        deleteMatrix(network->demand, network->numOrigins);
        network->demand = newMatrix(network->batchSize, network->numZones,
                                    double);
    }

}


void writeOBANetwork(network_type *network, char *linkFileName, 
        char *tripFileName) {
    int i, j;
    qsort(network->arcs, network->numArcs, sizeof(arc_type),
            forwardStarOrder); 
    displayMessage(FULL_NOTIFICATIONS, "Arcs sorted.");

    FILE* linkFile = openFile(linkFileName, "w");

    displayMessage(FULL_NOTIFICATIONS, "Opening %s and %s\n", linkFileName,
            tripFileName);
    fprintf(linkFile, "<NUMBER OF NODES> %d\n", network->numNodes);
    fprintf(linkFile, "<NUMBER OF LINKS> %d\n", network->numArcs);
    fprintf(linkFile, "<NUMBER OF ZONES> %d\n", network->numZones);
    fprintf(linkFile, "<FIRST THRU NODE> %d\n", network->firstThroughNode);
    fprintf(linkFile, "<END OF METADATA>\n");
    fprintf(linkFile, "\n");
    fprintf(linkFile, "~\tTail\tHead\tCapacity (veh/h)\tLength (ft)\tFree "
            "Flow Time (min)\tB\tPower\tSpeed (ft/min)\tToll\tType\t;\n");
    for (i = 0; i < network->numArcs; i++) {
        fprintf(linkFile,
                "\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t;\n",
            network->arcs[i].tail,
            network->arcs[i].head,
            network->arcs[i].capacity,
            network->arcs[i].length,
            network->arcs[i].freeFlowTime,
            network->arcs[i].alpha,
            network->arcs[i].beta,
            network->arcs[i].speedLimit,
            network->arcs[i].toll,
            network->arcs[i].linkType);
    }
    fclose(linkFile);
    displayMessage(FULL_NOTIFICATIONS, "Wrote arcs.\n");

    int destinationsWritten;
    FILE* tripFile = openFile(tripFileName, "w");

    fprintf(tripFile, "<NUMBER OF ZONES> %d\n", network->numZones);
    fprintf(tripFile, "<TOTAL OD FLOW> %lf\n", network->totalODFlow);
    fprintf(tripFile, "<END OF METADATA>\n");
    for (i = 0; i < network->numZones; i++) {
        fprintf(tripFile, "\n\nOrigin %d\n", i);
        destinationsWritten = 0;
        for (j = 0; j < network->numZones; j++) {
            if (network->demand[i][j] > 0) {
                fprintf(tripFile, "\t%d :\t%lf;",j,network->demand[i][j]);
                destinationsWritten++;
                if (destinationsWritten % DESTINATIONS_PER_LINE == 0)
                    fprintf(tripFile, "\n");
            }
        }
    }
    fclose(tripFile);
    displayMessage(FULL_NOTIFICATIONS, "Wrote ODs.\n");
}

/*
 * writeNetworkFlows prints the link IDs, flows, and costs in the file
 * format used on the TNTP site.
 */
void writeNetworkFlows(network_type *network, char *outputFileName) {
    FILE *outFile = openFile(outputFileName, "w");
    int ij;

    for (ij = 0; ij < network->numArcs; ij++) {
        fprintf(outFile, "(%d,%d) %f %f\n", network->arcs[ij].tail + 1,
                                            network->arcs[ij].head + 1,
                                            network->arcs[ij].flow,
                                            network->arcs[ij].cost
                                                - network->arcs[ij].fixedCost);
    }
    
    fclose(outFile);
}

/*
 * writeMulticlassFlows prints link IDs, total flows and costs, and then
 * class-specific flows, followed by summary information on total time,
 * distance, and toll paid by class and in total.
 *
 * Also includes a header row for multiclass.  (This is a "richer" version of
 * writeNetworkFlows).  If there is just a single class, it is the same
 * as writeNetworkFlows, but with time/distance/toll printed at bottom.
 */
void writeMulticlassFlows(network_type *network, char *outputFileName) {
    if (network->numClasses == 1) { /* just call regular writeNetworkFlows */
        writeNetworkFlows(network, outputFileName);
        /* But add summary statistics at the end */
        FILE *outFile = openFile(outputFileName, "a");
        fprintf(outFile, "Total time: %f\n", classTravelTime(network, 0));
        fprintf(outFile, "Total distance: %f\n", classDistance(network, 0));
        fprintf(outFile, "Total toll: %f\n", classRevenue(network, 0));
        fprintf(outFile, "Total generalized cost: %f\n", TSTT(network));
        fclose(outFile);
        return;
    }
    /* Usual case */
    FILE *outFile = openFile(outputFileName, "w");
    int c, ij;

    fprintf(outFile, "Link TotalFlow TravelTime");
    for (c = 0; c < network->numClasses; c++) {
        fprintf(outFile, " Class%d", c+1);
    }
    fprintf(outFile, "\n");
    for (ij = 0; ij < network->numArcs; ij++) {
        fprintf(outFile, "(%d,%d) %f %f", network->arcs[ij].tail + 1,
                                          network->arcs[ij].head + 1,
                                          network->arcs[ij].flow,
                                          network->arcs[ij].cost
                                              - network->arcs[ij].fixedCost);
        for (c = 0; c < network->numClasses; c++) {
            fprintf(outFile, " %f", network->arcs[ij].classFlow[c]);
        }
        fprintf(outFile, "\n");
    }

    fprintf(outFile, "\n\nClass Time Distance Toll GeneralizedCost\n");
    float classTime, classDist, classToll, classGC;
    float totalTime = 0, totalDistance = 0, totalToll = 0, totalCost = 0;
    for (c = 0; c < network->numClasses; c++) {
        classTime = classTravelTime(network, c);
        classDist = classDistance(network, c);
        classToll = classRevenue(network, c);
        classGC = classCost(network, c, 1, network->tollFactor[c],
                                           network->distanceFactor[c]);
        totalTime += classTime;
        totalDistance += classDist;
        totalToll += classToll;
        totalCost += classGC;
        fprintf(outFile, "%d %f %f %f %f\n", c+1, classTime, classDist,
                                             classToll, classGC);
    }
    fprintf(outFile, "TOTAL %f %f %f %f\n", totalTime, totalDistance,
                                            totalToll, totalCost);
    fclose(outFile);
}


    

///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length) {
    int i;
    for (i = 0; i < length; i++) string[i] = '\0';
}

int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue) {
    /* metadataTag and metadataValue both need to be of at least STRING_SIZE */
    int i = 0, j = 0;
    inputLine[STRING_SIZE-1] = '\0';
    while (inputLine[i] != '\0' && inputLine[i] != '\n' && inputLine[i] != '\r' 
            && inputLine[i] != '<' && inputLine[i] != '~') i++;
    if (inputLine[i] == '\0' || inputLine[i] == '\n' || inputLine[i] == '\r') 
        return BLANK_LINE;
    if (inputLine[i] == '~') return COMMENT;
    i++;
    while (inputLine[i] != '\0' && inputLine[i] != '>') {
        metadataTag[j++] = toupper(inputLine[i++]);
    }
    metadataTag[j] = '\0';
    if (inputLine[i] == '\0')
        fatalError("Metadata tag not closed: ", metadataTag);
    i++;
    while (inputLine[i] != '\0' 
            && (inputLine[i] == ' ' || inputLine[i] == '\t')) 
        i++;
    j = 0;
    while (inputLine[i] != '\0' && inputLine[i] != '\n' && inputLine[i] != '~')
        metadataValue[j++] = inputLine[i++];
    metadataValue[j] = '\0';
    return SUCCESS;
}

// Checks for comments and blank lines, and removes leading spaces
int parseLine(char* inputLine, char* outputLine) {
    int i = 0, j = 0;
    while (inputLine[i] != '\0' && (inputLine[i] == ' ' 
                || inputLine[i] == '\t')) i++;
    if (inputLine[i] == '~') return COMMENT;
    if (inputLine[i] == '\0' || inputLine[i] == '\n' 
            || inputLine[i] == '\r') return BLANK_LINE;
    while (inputLine[i] != '\0' && i < STRING_SIZE - 1) {
        outputLine[j++] = inputLine[i++];
    }
    outputLine[j] = '\0';
    return SUCCESS;
}

void parseCSV(char field[][STRING_SIZE], char *fullLine, int numFields) {
    int curField, pos;
    char *comma, *position;

    /* Strip trailing newline, if any */
    int len = strlen(fullLine);
    while (len > 0 && isspace(fullLine[len-1])) {
        fullLine[--len] = '\0';
    }

    comma = strchr(fullLine, ',');
    position = fullLine;
    curField = 0;
    while (comma) {
        if (curField >= numFields) fatalError("Too many args to parseCSV!");
        pos = 0;
        while (position < comma && pos <= STRING_SIZE) {
            field[curField][pos] = *position;
            pos++;
            position++;
        }
        field[curField][pos] = '\0';
        position++;
        curField++;
        comma = strchr(position, ',');
    }
    /* Copy last field */
    strncpy(field[curField], position, STRING_SIZE);
}

#ifdef PARALLELISM
/* Code indepdent taken from Stack Overflow (https://stackoverflow.com/a/3006416) */
int getNumCores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if(count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if(count < 1) { count = 1; }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}
#endif
