#include <unistd.h>
#include <stdio.h>
#include "fileio.h"

///////////////////////////
// Reading network files //
///////////////////////////

/* These routines are HARD-CODED for the NCTCOG network.  They cannot read
 * any other network or file.  Dimensions, etc. will be wrong.*/
#define NCTCOG_NODES 32782
#define NCTCOG_DIRECTED_LINKS 85089
#define NCTCOG_UNDIRECTED_LINKS 50295
#define NCTCOG_ZONES 5352
#define NCTCOG_MAX_NODE_ID 77429
#define NCTCOG_NUM_CLASSES 10
void readNCTCOGNetwork(network_type *network, char *networkFileName,
                       char *tripFileName, char *converterFileName) {
    int r, s, c;
    int NCTCOG2SDB[NCTCOG_MAX_NODE_ID+1]; /* Conversion table */

    /* Set up network and data structures */
    network->numNodes = NCTCOG_NODES;
    network->numArcs = NCTCOG_DIRECTED_LINKS;
    network->numZones = NCTCOG_ZONES;
    network->firstThroughNode = NCTCOG_ZONES + 1;
    network->numClasses = NCTCOG_NUM_CLASSES;
    network->numOrigins = network->numZones * network->numClasses;
    network->batchSize = network->numZones;
    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->arc_muts = newVector(network->numArcs, pthread_mutex_t);

    if (tripFileName == NULL) { /* Warm-start = compact matrix */
        network->demand = newMatrix(network->numZones, network->numZones,
                                    double);
    } else { /* Cold start = full matrix */
        network->demand = newMatrix(network->numOrigins, network->numZones,
                                    double);
    }

    network->tollFactor = newVector(network->numClasses, double);
    network->distanceFactor = newVector(network->numClasses, double);

    for (r = 0;
         r < (tripFileName == NULL ? network->numZones : network->numOrigins);
         r++) {
        for (s = 0; s < network->numZones; s++) {
            network->demand[r][s] = 0;
        }
    }

    for (c = 0; c < network->numClasses; c++) network->distanceFactor[c] = 0;
    /* Toll factors are *inverse* of VOT */
    network->tollFactor[SOLO_17] = 60.0 / 10.5; 
    network->tollFactor[SOLO_35] = 60.0 / 21.0;
    network->tollFactor[SOLO_45] = 60.0 / 27.0;
    network->tollFactor[SOLO_90] = 60.0 / 54.0;
    network->tollFactor[HOV_17] = 60.0 / 10.5;
    network->tollFactor[HOV_35] = 60.0 / 21.0;
    network->tollFactor[HOV_45] = 60.0 / 27.0;
    network->tollFactor[HOV_90] = 60.0 / 54.0;
    network->tollFactor[MED_TRUCKS] = 60.0 / 60.0;
    network->tollFactor[HVY_TRUCKS] = 60.0 / 60.0;

    /* Now read files */
    readConverterFile(converterFileName, NCTCOG2SDB, network->numNodes,
                      sizeof(NCTCOG2SDB)/sizeof(NCTCOG2SDB[0]), FALSE);
    readNCTCOGLinks(network, networkFileName, NCTCOG2SDB);
    if (tripFileName != NULL && strcmp("STREAM", tripFileName) == 0) {
        streamNCTCOGTrips(network, NCTCOG2SDB);
    } else if (tripFileName != NULL) {
        /* If you give a NULL tripFileName, we assume a warm start and skip
         * this file. */
        readNCTCOGTrips(network, tripFileName, NCTCOG2SDB);
    }

    finalizeNetwork(network);



    /* If trip file was read, write the batch OD matrices so we don't have
     * to do so again */
//    if (tripFileName != NULL) {
//        writeBinaryMatrices(network);
//        deleteMatrix(network->demand, network->numOrigins);
//        network->demand=newMatrix(network->batchSize,network->numZones,double);
//    }
    displayMessage(FULL_NOTIFICATIONS, "Finished reading NCTCOG Network and trips\n");
}

/*
 * Write OD matrices into binary files, one for each batch.
 * For use in later warm-starts or assignments.
 */
void writeBinaryMatrices(network_type *network) {
    int batch, r, origin;
    char filename[STRING_SIZE];
    FILE* matrixFile;

    for (batch = 0; batch < network->numBatches; batch++) {
        sprintf(filename, "matrix%d.bin", batch);
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


void readConverterFile(char *converterFileName, int *table, int maxValue,
                       int maxKey, bool isInverse) {
    int i;
    int NCTCOG, SDB;

    FILE *converterFile = openFile(converterFileName, "r");
    for (i = 0; i < maxKey; i++) table[i] = IS_MISSING;
    while (fscanf(converterFile, "%d,%d", &SDB, &NCTCOG) != EOF) {
        if (isInverse == TRUE) swap(NCTCOG, SDB);
        if (NCTCOG >= 0 && NCTCOG <= maxKey && SDB >= 0 && SDB <= maxValue) {
            table[NCTCOG] = SDB;
        } else {
            warning(DEBUG,
                    "Converter tables out of range: %d -> %d\n", NCTCOG, SDB);
        }
    }
    fclose(converterFile);
}

int convert(int value, int *table, int maxKey) {
    if (value < 0 || value > maxKey || table[value] == IS_MISSING) {
        fatalError("Invalid table lookup: %d (max %d)\n", value, maxKey);
    }
    return table[value];
}

void readNCTCOGLinks(network_type *network, char *networkFileName, int *table){
    int e, ij;
    char lineData[NUM_NCTCOG_NET_COLUMNS][STRING_SIZE];
    char fullLine[STRING_SIZE];
    FILE *networkFile = openFile(networkFileName, "r");

    /* Skip header row */
    if (fgets(fullLine, STRING_SIZE, networkFile) == NULL)
        fatalError("Cannot read header of network file %s", networkFileName);

    /* Now read for each link */
    ij = 0; /* ij counts directional links, e undirected links */
    for (e = 1; e < NCTCOG_UNDIRECTED_LINKS; e++) {
        if (fgets(fullLine, STRING_SIZE, networkFile) == NULL)
            fatalError("Link file done before all links read.");
        parseCSV(lineData, fullLine, NUM_NCTCOG_NET_COLUMNS);
        if (atof(lineData[PMCAP_AB]) > 0) { /* Create AB link (handles '--')*/
            makeLink(network, ij++, table, lineData[ID], lineData[FROM_NODE],
                lineData[TO_NODE], lineData[PMCAP_AB], lineData[LENGTH],
                lineData[PKFRTIME_AB], lineData[A_PK_CONICAL],
                lineData[VDF_SHIFT], lineData[SPAR_AB], lineData[PMSATFLOW_AB],
                lineData[CA_AB], lineData[CB_AB], lineData[CC_AB],
                lineData[CD_AB], lineData[UNSIG_MINDELAY], lineData[UPAR_AB],
                lineData[OPERCOSTPM_AB], lineData[TOLL_AUTO_DA_PM_AB],
                lineData[TOLL_AUTO_SR_PM_AB],
                lineData[TOLL_MEDIUM_TRUCK_PM_AB],
                lineData[TOLL_HEAVY_TRUCK_PM_AB], lineData[FLAG_EXCLUDE_DA]);
        }
        if (atof(lineData[PMCAP_BA]) > 0) { /* Create BA link */
            makeLink(network, ij++, table, lineData[ID], lineData[TO_NODE],
                lineData[FROM_NODE], lineData[PMCAP_BA], lineData[LENGTH],
                lineData[PKFRTIME_BA], lineData[A_PK_CONICAL],
                lineData[VDF_SHIFT], lineData[SPAR_BA], lineData[PMSATFLOW_BA],
                lineData[CA_BA], lineData[CB_BA], lineData[CC_BA],
                lineData[CD_BA], lineData[UNSIG_MINDELAY], lineData[UPAR_BA],
                lineData[OPERCOSTPM_BA], lineData[TOLL_AUTO_DA_PM_BA],
                lineData[TOLL_AUTO_SR_PM_BA],
                lineData[TOLL_MEDIUM_TRUCK_PM_BA],
                lineData[TOLL_HEAVY_TRUCK_PM_BA], lineData[FLAG_EXCLUDE_DA]);
        }
    }
    fclose(networkFile);
}

void makeLink(network_type *network, int ij, int *table, char *ID, char *from,
        char *to, char *cap, char *len, char *freeFlow, char *conical, char
        *shift, char *sPar, char *satFlow, char *CA, char *CB, char *CC, char
        *CD, char *minDelay, char *uPar, char *operCost, char *daToll, char
        *srToll, char *medToll, char *heavyToll, char *exclude) {
    int c;

    network->arcs[ij].ID = atoi(ID);
    network->arcs[ij].tail = convert(atoi(from), table, NCTCOG_MAX_NODE_ID);
    network->arcs[ij].head = convert(atoi(to), table, NCTCOG_MAX_NODE_ID);
    network->arcs[ij].flow = 0;
    network->arcs[ij].cost = atof(freeFlow);
    network->arcs[ij].der = 0;
    network->arcs[ij].freeFlowTime = atof(freeFlow);
    network->arcs[ij].capacity = atof(cap);
    network->arcs[ij].a = atof(conical);
    network->arcs[ij].e = atof(shift);
    network->arcs[ij].sParam = atof(sPar);
    network->arcs[ij].saturationFlow = atof(satFlow);
    network->arcs[ij].CA = atof(CA);
    network->arcs[ij].CB = atof(CB);
    network->arcs[ij].CC = atof(CC);
    network->arcs[ij].CD = atof(CD);
    network->arcs[ij].m = atof(minDelay);
    network->arcs[ij].u = atof(uPar);
    network->arcs[ij].operCost = atof(operCost);
    network->arcs[ij].toll_DA = atof(daToll);
    network->arcs[ij].toll_SR = atof(srToll);
    network->arcs[ij].toll_Med = atof(medToll);
    network->arcs[ij].toll_Heavy = atof(heavyToll);
    network->arcs[ij].length = atof(len);
    if (strcmp(exclude, "TRUE") == 0) {
        network->arcs[ij].excludeDA = TRUE;
    } else if (strcmp(exclude, "FALSE") == 0) {
        network->arcs[ij].excludeDA = FALSE;
    } else {
        fatalError("Unreadable exclude parameter '%s'", exclude);
    }

    network->arcs[ij].alpha = 0.15;
    network->arcs[ij].beta = 4;
    network->arcs[ij].toll = 0;
    network->arcs[ij].speedLimit = IS_MISSING;
    network->arcs[ij].linkType = IS_MISSING;
    network->arcs[ij].fixedCost = 0;
    if (network->arcs[ij].saturationFlow > 0) { /* Regular link */
        network->arcs[ij].calculateCost = &conicCost;
        network->arcs[ij].calculateDer = &conicDer;
        network->arcs[ij].calculateInt = &conicInt;
    } else { /* Centroid connector */
        network->arcs[ij].alpha = 0;
        network->arcs[ij].beta = 1;
        network->arcs[ij].calculateCost = &linearBPRcost;
        network->arcs[ij].calculateDer = &linearBPRder;
        network->arcs[ij].calculateInt = &linearBPRint;
    }

    network->arcs[ij].classFlow = newVector(network->numClasses, double);
    network->arcs[ij].classCost = newVector(network->numClasses, double);
    network->arcs[ij].classToll = newVector(network->numClasses, double);

    for (c = 0; c < network->numClasses; c++) {
        if (isDA(c)) {
            network->arcs[ij].classToll[c] = network->arcs[ij].toll_DA;
        } else if (isHOV(c)) {
            network->arcs[ij].classToll[c] = network->arcs[ij].toll_SR;
        } else if (c == MED_TRUCKS) {
            network->arcs[ij].classToll[c] = network->arcs[ij].toll_Med;
        } else if (c == HVY_TRUCKS) {
            network->arcs[ij].classToll[c] = network->arcs[ij].toll_Heavy;
        } else {
            fatalError("Class index and NCTCOG_classes typedef not aligned.");
        }
        if (isSolo(c) && network->arcs[ij].excludeDA == TRUE) {
            network->arcs[ij].classCost[c] = ARTIFICIAL;
        } else {
            network->arcs[ij].classCost[c] = network->tollFactor[c] 
            * (network->arcs[ij].operCost + network->arcs[ij].classToll[c]);
        }
        network->arcs[ij].classFlow[c] = 0;
    }
}

void readNCTCOGTrips(network_type *network, char *tripFileName, int *table) {
    int r, s;
    char lineData[NUM_NCTCOG_CLASSES + 2][STRING_SIZE];
    char fullLine[STRING_SIZE];
    FILE *tripFile = openFile(tripFileName, "r");
    while (TRUE) { /* Break in middle when out of lines */
        if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) break;
        if (strstr(fullLine, ",") == NULL) continue;
        parseCSV(lineData, fullLine, NUM_NCTCOG_CLASSES + 2);
        r = convert(atoi(lineData[0]), table, NCTCOG_MAX_NODE_ID);
        s = convert(atoi(lineData[1]), table, NCTCOG_MAX_NODE_ID);
        assignDemand(network, r, s, SOLO_35, lineData[2]);
        assignDemand(network, r, s, SOLO_90, lineData[3]);
        assignDemand(network, r, s, HOV_35, lineData[4]);
        assignDemand(network, r, s, HOV_90, lineData[5]);
        assignDemand(network, r, s, SOLO_17, lineData[6]);
        assignDemand(network, r, s, SOLO_45, lineData[7]);
        assignDemand(network, r, s, HOV_17, lineData[8]);
        assignDemand(network, r, s, HOV_45, lineData[9]);
        assignDemand(network, r, s, MED_TRUCKS, lineData[10]);
        assignDemand(network, r, s, HVY_TRUCKS, lineData[11]);
    } ;
    fclose(tripFile);
}

void streamNCTCOGTrips(network_type *network, int *table) {
    ssize_t n;
    ssize_t off = 0;
    int r, s;
//    char initial[6];
    // char *buffer = (char *) calloc(48, sizeof(char));
    float buffer[12];
    displayMessage(FULL_NOTIFICATIONS, "Starting to read...\n");
    displayMessage(FULL_NOTIFICATIONS, "Printing before the tight loop\n");
    int count = 0;
    do {
        n = read(STDIN_FILENO, ((char*)buffer) + off, 48 - off);  /* Break in middle when out of lines */
        off += n;
        if (n < 0) {
            fatalError("Issue reading from stdin");
        } else if (strcmp((char *) buffer, "Sto") == 0) {
            break;
        } else if (off != 48){
            continue;
        }
//        displayMessage(FULL_NOTIFICATIONS, "read the following bytes %d\n", off);
        off = 0;
//        displayMessage(FULL_NOTIFICATIONS, "Printing before initial conversion\n");
//        displayMessage(FULL_NOTIFICATIONS, "r: %d, s: %d, SOLO_35: %f, SOLO_90: %f, HOV_35: %f, HOV_90: %f, SOLO_17: %f, SOLO_45: %f, HOV_17: %f, HOV_45: %f, MED_TRUCKS: %f, HVY_TRUCKS: %f\n",
//                       ((int *)buffer)[0], ((int *) buffer)[1],
//                       ((float *) buffer)[2], ((float *) buffer)[3],
//                       ((float *) buffer)[4], ((float *) buffer)[5],
//                       ((float *) buffer)[6], ((float *) buffer)[7],
//                       ((float *) buffer)[8], ((float *) buffer)[9],
//                       ((float *) buffer)[10], ((float *) buffer)[11]);
        count += 1;
        r = convert(((int *)buffer)[0], table, NCTCOG_MAX_NODE_ID);
        s = convert(((int *)buffer)[1], table, NCTCOG_MAX_NODE_ID);
        assignStreamedDemand(network, r, s, SOLO_35, ((float *) buffer)[2]);
        assignStreamedDemand(network, r, s, SOLO_90, ((float *) buffer)[3]);
        assignStreamedDemand(network, r, s, HOV_35, ((float *) buffer)[4]);
        assignStreamedDemand(network, r, s, HOV_90, ((float *) buffer)[5]);
        assignStreamedDemand(network, r, s, SOLO_17, ((float *) buffer)[6]);
        assignStreamedDemand(network, r, s, SOLO_45, ((float *) buffer)[7]);
        assignStreamedDemand(network, r, s, HOV_17, ((float *) buffer)[8]);
        assignStreamedDemand(network, r, s, HOV_45, ((float *) buffer)[9]);
        assignStreamedDemand(network, r, s, MED_TRUCKS, ((float *) buffer)[10]);
        assignStreamedDemand(network, r, s, HVY_TRUCKS, ((float *) buffer)[11]);
        } while(1);
    displayMessage(FULL_NOTIFICATIONS, "Finished reading %d OD pairs from buffer\n", count);
}

void assignDemand(network_type *network, int originNode, int destinationNode,
                  int class, char *demandValue) {
    int origin = nodeclass2origin(network, originNode, class);
    network->demand[origin][destinationNode] = atof(demandValue);
}

void assignStreamedDemand(network_type *network, int originNode, int destinationNode,
                  int class, float demandValue) {
    int origin = originNode + class * network->batchSize;
    network->demand[origin][destinationNode] = demandValue;
}

void readSDBNetwork(network_type *network, char *filename, double defaultAlpha,
        double defaultBeta) {
    int i, j, c;
    int numParams;

    FILE *networkFile = openFile(filename, "r");

    /* Read header row and dimension arrays */
    numParams = fscanf(networkFile, "%d,%d,%d,%d",
                        &(network->numNodes),
                        &(network->numArcs),
                        &(network->numZones),
                        &(network->firstThroughNode));
    if (numParams != 4)
        fatalError("Header information incorrect in network file %s",filename);
    displayMessage(MEDIUM_NOTIFICATIONS, "Nodes, arcs, zones, thrunode: "
            "%d %d %d %d\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);

    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->demand = newMatrix(network->numZones, network->numZones, double);
    network->numClasses = 1;
    network->numOrigins = network->numZones * network->numClasses;
    // This file format cannot include distance or toll factors
    network->distanceFactor = newVector(network->numClasses, double);
    network->tollFactor = newVector(network->numClasses, double);
    for (c = 0; c < network->numClasses; c++) {
        network->distanceFactor[c] = 0;
        network->tollFactor[c] = 0;
    }
    network->totalODFlow = 0;

    network->beckmann = INFINITY;
    network->beckmannLB = INFINITY;

    /* Read arc and OD data */
    for(i = 0; i < network->numArcs; i++) {
        numParams = fscanf(networkFile, "%d,%d,%lf,%lf",
                &network->arcs[i].tail, &network->arcs[i].head,
                &network->arcs[i].capacity, &network->arcs[i].freeFlowTime);
        if (numParams != 4)
            fatalError("Unable to read information for arc %d in network "
                       "file %s.", i, filename);
        if (network->arcs[i].tail < 0 ||
                network->arcs[i].tail >= network->numNodes)
            fatalError("Arc %d tail %d out of range in network file %s.", i,
                    network->arcs[i].tail, filename);
        if (network->arcs[i].head < 0 || 
                network->arcs[i].head >= network->numNodes) 
            fatalError("Arc %d head %d out of range in network file %s.", i, 
                    network->arcs[i].head, filename);
        if (network->arcs[i].freeFlowTime < 0) 
            fatalError("Arc free flow time %d negative in network file %s.", i,
                    filename);
        if (network->arcs[i].capacity <= 0) 
            fatalError("Capacity %d nonpositive in network file %s.", i, 
                    filename);
        if (feof(networkFile)) 
            fatalError("network file ended after reading arc %ld.", i);
        network->arcs[i].flow = 0;
        network->arcs[i].cost = network->arcs[i].freeFlowTime;
        network->arcs[i].alpha = defaultAlpha;
        network->arcs[i].beta = defaultBeta;
        /* This file format does not include the following link fields */
        network->arcs[i].toll = 0;
        network->arcs[i].speedLimit = IS_MISSING;
        network->arcs[i].length = 0;
        network->arcs[i].linkType = IS_MISSING;

        network->arcs[i].classFlow = newVector(network->numClasses, double);
        network->arcs[i].classCost = newVector(network->numClasses, double);
        network->arcs[i].classToll = newVector(network->numClasses, double);        
        for (c = 0; c < network->numClasses; c++) {
            /* classFlow and classCost are initialized in finalizeNetwork */
            network->arcs[i].classToll[c] = network->arcs[i].toll;
        }

        
        if (defaultBeta == 1) {
            network->arcs[i].calculateCost = &linearBPRcost;
            network->arcs[i].calculateDer = &linearBPRder;           
            network->arcs[i].calculateInt = &linearBPRint;           
        } else if (defaultBeta == 4) {
            network->arcs[i].calculateCost = &quarticBPRcost;
            network->arcs[i].calculateDer = &quarticBPRder;           
            network->arcs[i].calculateInt = &quarticBPRint;           
        } else {
            network->arcs[i].calculateCost = &generalBPRcost;
            network->arcs[i].calculateDer = &generalBPRder;           
            network->arcs[i].calculateInt = &generalBPRint;           
        }       
    }
    for (i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones; j++) {
            if (fscanf(networkFile, "%lf,", &network->demand[i][j]) != 1)
                fatalError("Error reading OD table entry %d -> %d", i, j);
            if (feof(networkFile) && j < network->numZones 
                    && i < network->numZones)
                fatalError("network file %s ended before OD table is "
                            "complete.", filename);
            network->totalODFlow += network->demand[i][j];
        }
    }
    fclose(networkFile);
    displayMessage(LOW_NOTIFICATIONS,"File read and memory allocated.\n");

    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS,"Forward/reverse stars generated.\n");
}

void readOBANetwork(network_type *network, char *linkFileName, 
                    char **tripFileName, int numClasses,
                    algorithmBParameters_type *parameters) {
    int i, j, r = 0, c;
    int check;
    int numParams, status;
    double demand, totalDemandCheck = 0;
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
                    "in  parameters file %s", metadataTag, linkFileName);
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

    displayMessage(MEDIUM_NOTIFICATIONS, "Nodes, arcs, zones, thrunode: "
            "%ld %ld %ld %ld\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);
    displayMessage(MEDIUM_NOTIFICATIONS, "Distance factor, toll factor: "
            "%lf %lf\n", defaultDistanceFactor, defaultTollFactor);

    network->numOrigins = network->numZones * network->numClasses;
    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->demand = newMatrix(network->numOrigins, network->numZones,double);
    network->tollFactor = newVector(network->numClasses, double);
    network->distanceFactor = newVector(network->numClasses, double);
    network->arc_muts = newVector(network->numArcs, pthread_mutex_t);

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
        tripFile = openFile(tripFileName[c], "r");
        /* Verify trip table metadata */
        endofMetadata = FALSE;
        network->totalODFlow = IS_MISSING;
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
                network->totalODFlow = atof(metadataValue);
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
            demandMultiplier = parameters->demandMultiplier;

        /* Now read trip table */
        while (!feof(tripFile)) {
            if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) break;
            status = parseLine(fullLine, trimmedLine);
            if (status == BLANK_LINE || status == COMMENT) continue;
            if (strstr(trimmedLine, "Origin") != NULL) {
                // i indexes current origin
                sscanf(strstr(trimmedLine, "Origin")+6,"%d", &i);  
                if (i <= 0 || i > network->numNodes) fatalError("Origin %d is"
                        "out of range in trips file %s", i, tripFileName[c]);
                i--;
                r = nodeclass2origin(network, i, c);
                continue;
            }
            token = strtok(trimmedLine , ";");
            while (token != NULL && strlen(token) > 1) {
                numParams = sscanf(token, "%d : %lf", &j, &demand);
                if (numParams < 2) break;
                if (j <= 0 || j > network->numNodes) fatalError("Destination "
                        "%d is out of range in trips file %s\n%s\n%s", j, 
                        tripFileName[c], fullLine, token);
                j--;
                network->demand[r][j] = demand * demandMultiplier;
                if (demand < 0) fatalError("Negative demand from origin %d to "
                        "destination %d in class %d", i, j, c);
                totalDemandCheck += network->demand[r][j];
                token = strtok(NULL, ";");
            }
            blankInputString(trimmedLine, STRING_SIZE);
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

void writeDUENetwork(network_type *network, char *filename) {
    FILE* outfile = openFile(filename, "w");
    displayMessage(FULL_NOTIFICATIONS, "Opening %s\n", filename);
    int i, j;
    fprintf(outfile, "%d,%d,%d,%d\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);
    for (i = 0; i < network->numArcs; i++) {
        fprintf(outfile, "%d,%d,%f,%f\n", network->arcs[i].tail, 
                network->arcs[i].head, network->arcs[i].capacity, 
                network->arcs[i].freeFlowTime);
    }
    displayMessage(FULL_NOTIFICATIONS, "Wrote arcs.\n");
    for (i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones - 1; j++) {
            fprintf(outfile, "%f,", network->demand[i][j]);
        }
        fprintf(outfile, "%f\n", network->demand[i][j]);
    }
    displayMessage(FULL_NOTIFICATIONS, "Wrote ODs.\n");
    fclose(outfile);
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




///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length) {
    int i;
    for (i = 0; i < length; i++) string[i] = '\0';
}

int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue) {
    int i = 0, j = 0;
    while (inputLine[i] != 0 && inputLine[i] != '\n' && inputLine[i] != '\r' 
            && inputLine[i] != '<') i++;
    if (inputLine[i] == 0 || inputLine[i] == '\n' || inputLine[i] == '\r') 
        return BLANK_LINE;
    if (inputLine[i] == '~') return COMMENT;
    i++;
    while (inputLine[i] != 0 && inputLine[i] != '>') {
        metadataTag[j++] = toupper(inputLine[i++]);
    }
    metadataTag[j] = 0;
    if (inputLine[i] == 0) fatalError("Metadata tag not closed in parameters "
            "file - ", metadataTag);
    i++;
    while (inputLine[i] != 0 && (inputLine[i] == ' ' || inputLine[i] == '\t')) 
        i++;
    j = 0;
    while (inputLine[i] != 0 && inputLine[i] != '\n') {
        metadataValue[j++] = inputLine[i++];
    }
    metadataValue[j] = 0;
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
    while (inputLine[i] != '\0') {
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
