#include "fileio.h"

///////////////////////////
// Reading network files //
///////////////////////////

FILE *openFile(const char *filename, const char *access) {
    FILE *handle = fopen(filename, access);
    if (handle == NULL) fatalError("File %s not found", filename);
    return handle;
}

void readSDBNetwork(network_type *network, char *filename, double defaultAlpha,
        double defaultBeta) {
    long i, j;
    int numParams;

    FILE *networkFile = openFile(filename, "r");

    /* Read header row and dimension arrays */
    numParams = fscanf(networkFile, "%ld,%ld,%ld,%ld", &(network->numNodes), &(network->numArcs), &(network->numZones), &(network->firstThroughNode));
    if (numParams != 4)
        fatalError("Header information incorrect in network file %s",filename);
    displayMessage(MEDIUM_NOTIFICATIONS, "Nodes, arcs, zones, thrunode: "
            "%ld %ld %ld %ld\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);

    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->OD = newMatrix(network->numZones, network->numZones, od_type);
    // This file format cannot include distance or toll factors
    network->distanceFactor = 0;
    network->tollFactor = 0;
    network->totalODFlow = 0;

    network->beckmann = INFINITY;
    network->beckmannLB = INFINITY;

    /* Read arc and OD data */
    for(i = 0; i < network->numArcs; i++) {
        numParams = fscanf(networkFile, "%ld,%ld,%lf,%lf",
                &network->arcs[i].tail, &network->arcs[i].head,
                &network->arcs[i].capacity, &network->arcs[i].freeFlowTime);
        if (numParams != 4)
            fatalError("Unable to read information for arc %ld in network "
                       "file %s.", i, filename);
        if (network->arcs[i].tail <= 0 ||
                network->arcs[i].tail > network->numNodes)
            fatalError("Arc tail %d out of range in network file %s.", i,
                    filename);
        if (network->arcs[i].head <= 0 || 
                network->arcs[i].head > network->numNodes) 
            fatalError("Arc head %d out of range in network file %s.", i, 
                    filename);
        if (network->arcs[i].freeFlowTime < 0) 
            fatalError("Arc free flow time %d negative in network file %s.", i,
                    filename);
        if (network->arcs[i].capacity <= 0) 
            fatalError("Capacity %d nonpositive in network file %s.", i, 
                    filename);
        if (feof(networkFile)) 
            fatalError("network file ended after reading arc %ld.", i);
        network->arcs[i].head--;
        network->arcs[i].tail--;
        network->arcs[i].flow = 0;
        network->arcs[i].cost = network->arcs[i].freeFlowTime;
        network->arcs[i].alpha = defaultAlpha;
        network->arcs[i].beta = defaultBeta;
        /* This file format does not include the following link fields */
        network->arcs[i].toll = IS_MISSING;
        network->arcs[i].speedLimit = IS_MISSING;
        network->arcs[i].length = IS_MISSING;
        network->arcs[i].linkType = IS_MISSING;
        
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
    for(i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones; j++) {
            if (fscanf(networkFile, "%lf,", &network->OD[i][j].demand) != 1)
                fatalError("Error reading OD table entry %d -> %d", i, j);
            if (feof(networkFile) && j < network->numZones 
                    && i < network->numZones)
                fatalError("network file %s ended before OD table is "
                            "complete.", filename);
            network->totalODFlow += network->OD[i][j].demand;
        }
    }
    fclose(networkFile);
    displayMessage(LOW_NOTIFICATIONS,"File read and memory allocated.\n");

    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS,"Forward/reverse stars generated.\n");
}

void readOBANetwork(network_type *network, char *linkFileName, 
        char *tripFileName) {
    long i, j;
    long check;
    int numParams, status;
    double demand, demandMultiplier, totalDemandCheck = 0;

    char fullLine[STRING_SIZE], trimmedLine[STRING_SIZE], *token;
    char metadataTag[STRING_SIZE], metadataValue[STRING_SIZE];

    FILE *linkFile = openFile(linkFileName, "r");
    FILE *tripFile = openFile(tripFileName, "r");

    network->numZones = IS_MISSING;
    network->numArcs = IS_MISSING;
    network->numNodes = IS_MISSING;
    network->firstThroughNode = IS_MISSING;
    network->tollFactor = IS_MISSING;
    network->distanceFactor = IS_MISSING;

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
            network->distanceFactor = atof(metadataValue);
        } else if (strcmp(metadataTag, "TOLL FACTOR") == 0) {
            network->tollFactor = atof(metadataValue);
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
    if (network->distanceFactor == IS_MISSING) {
        warning(LOW_NOTIFICATIONS, "Link file %s does not contain distance "
                "factor, setting to 0 as default.\n", linkFileName);
        network->distanceFactor = 0;
    }
    if (network->tollFactor == IS_MISSING) {
        warning(LOW_NOTIFICATIONS, "Link file %s does not contain toll factor,"
                " setting to 0 as default.\n", linkFileName);
        network->tollFactor = 0;
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
            "%lf %lf\n", network->distanceFactor, network->tollFactor);

    network->nodes = newVector(network->numNodes, node_type);
    network->arcs = newVector(network->numArcs, arc_type);
    network->OD = newMatrix(network->numZones, network->numZones, od_type);


    for (i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones; j++) {
            network->OD[i][j].demand = 0;
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
        numParams=sscanf(trimmedLine,"%ld %ld %lf %lf %lf %lf %lf %lf %lf %d",
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
            fatalError("Arc length %d negative in network file %s.\n%s", i,
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
      }           }

    /* Verify trip table metadata */
    endofMetadata = FALSE;
    demandMultiplier = 1;
    network->totalODFlow = IS_MISSING;
    do {
        if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) 
            fatalError("Trip file %s ended (or other I/O error) before "
                    "metadata complete.", tripFileName);
        status = parseMetadata(fullLine, metadataTag, metadataValue);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if         (strcmp(metadataTag, "NUMBER OF ZONES") == 0) {
            check = atoi(metadataValue);
            if (check != network->numZones) fatalError("Number of zones in "
                    "trip and link files do not match.");
        } else if (strcmp(metadataTag, "TOTAL OD FLOW") == 0) {
            network->totalODFlow = atof(metadataValue);
        } else if (strcmp(metadataTag, "DEMAND MULTIPLIER") == 0) {
            demandMultiplier = atof(metadataValue);
        } else if (strcmp(metadataTag, "END OF METADATA") == 0) {
            endofMetadata = TRUE;
        } else {
            warning(MEDIUM_NOTIFICATIONS, "Ignoring unknown metadata tag %s "
                    "in trips file %s", metadataTag, tripFileName);
        }
    } while (endofMetadata == FALSE);

    /* Now read trip table */
    while (!feof(tripFile)) {
        if (fgets(fullLine, STRING_SIZE, tripFile) == NULL) break;
        status = parseLine(fullLine, trimmedLine);
        if (status == BLANK_LINE || status == COMMENT) continue;
        if (strstr(trimmedLine, "Origin") != NULL) {
            // i indexes current origin
            sscanf(strstr(trimmedLine, "Origin")+6,"%ld", &i);  
            if (i <= 0 || i > network->numNodes) fatalError("Origin %d is out "
                    "of range in trips file %s", i, tripFileName);
            i--;
            continue;
        }
        token = strtok(trimmedLine , ";");
        while (token != NULL && strlen(token) > 1) {
            numParams = sscanf(token, "%ld : %lf", &j, &demand);
            if (numParams < 2) break;
            if (j <= 0 || j > network->numNodes) fatalError("Destination %d "
                    "is out of range in trips file %s\n%s\n%s", j, 
                    tripFileName, fullLine, token);
            j--;
            network->OD[i][j].demand = demand * demandMultiplier;
            if (demand < 0) fatalError("Negative demand from origin %d to "
                    "destination %d", i, j);
            totalDemandCheck += network->OD[i][j].demand;
            token = strtok(NULL, ";");
        }
        blankInputString(trimmedLine, STRING_SIZE);
    }

    displayMessage(MEDIUM_NOTIFICATIONS, "%d nodes, %d arcs, and %d zones\n",
            network->numNodes, network->numArcs, network->numZones);
    displayMessage(FULL_NOTIFICATIONS, "Total demand %f compared to metadata "
            "%f\n", totalDemandCheck, network->totalODFlow);
    //network->totalODFlow = totalDemandCheck;
    fclose(linkFile);
    fclose(tripFile);
    displayMessage(LOW_NOTIFICATIONS, "File read and memory allocated.\n");

    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS, "Forward and reverse star lists "
            "generated.\n");

}

void writeDUENetwork(network_type *network, char *filename) {
    FILE* outfile = openFile(filename, "w");
    displayMessage(FULL_NOTIFICATIONS, "Opening %s\n", filename);
    long i, j;
    fprintf(outfile, "%ld,%ld,%ld,%ld\n", network->numNodes, network->numArcs,
            network->numZones, network->firstThroughNode);
    for (i = 0; i < network->numArcs; i++) {
        fprintf(outfile, "%ld,%ld,%f,%f\n", network->arcs[i].tail, 
                network->arcs[i].head, network->arcs[i].capacity, 
                network->arcs[i].freeFlowTime);
    }
    displayMessage(FULL_NOTIFICATIONS, "Wrote arcs.\n");
    for (i = 0; i < network->numZones; i++) {
        for (j = 0; j < network->numZones - 1; j++) {
            fprintf(outfile, "%f,", network->OD[i][j].demand);
        }
        fprintf(outfile, "%f\n", network->OD[i][j].demand);
    }
    displayMessage(FULL_NOTIFICATIONS, "Wrote ODs.\n");
    fclose(outfile);
}


void writeOBANetwork(network_type *network, char *linkFileName, 
        char *tripFileName) {
    long i, j;
    qsort(network->arcs, network->numArcs, sizeof(arc_type),
            forwardStarOrder); 
    displayMessage(FULL_NOTIFICATIONS, "Arcs sorted.");

    FILE* linkFile = openFile(linkFileName, "w");

    displayMessage(FULL_NOTIFICATIONS, "Opening %s and %s\n", linkFileName,
            tripFileName);
    fprintf(linkFile, "<NUMBER OF NODES> %ld\n", network->numNodes);
    fprintf(linkFile, "<NUMBER OF LINKS> %ld\n", network->numArcs);
    fprintf(linkFile, "<NUMBER OF ZONES> %ld\n", network->numZones);
    fprintf(linkFile, "<FIRST THRU NODE> %ld\n", network->firstThroughNode);
    fprintf(linkFile, "<DISTANCE FACTOR> %lf\n", network->distanceFactor);
    fprintf(linkFile, "<TOLL FACTOR> %lf\n", network->tollFactor);
    fprintf(linkFile, "<END OF METADATA>\n");
    fprintf(linkFile, "\n");
    fprintf(linkFile, "~\tTail\tHead\tCapacity (veh/h)\tLength (ft)\tFree "
            "Flow Time (min)\tB\tPower\tSpeed (ft/min)\tToll\tType\t;\n");
    for (i = 0; i < network->numArcs; i++) {
        fprintf(linkFile,
                "\t%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t;\n",
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

    long destinationsWritten;
    FILE* tripFile = openFile(tripFileName, "w");

    fprintf(tripFile, "<NUMBER OF ZONES> %ld\n", network->numZones);
    fprintf(tripFile, "<TOTAL OD FLOW> %lf\n", network->totalODFlow);
    fprintf(tripFile, "<END OF METADATA>\n");
    for (i = 0; i < network->numZones; i++) {
        fprintf(tripFile, "\n\nOrigin %ld\n", i);
        destinationsWritten = 0;
        for (j = 0; j < network->numZones; j++) {
            if (network->OD[i][j].demand > 0) {
                fprintf(tripFile, "\t%ld :\t%lf;",j,network->OD[i][j].demand);
                destinationsWritten++;
                if (destinationsWritten % DESTINATIONS_PER_LINE == 0)
                    fprintf(tripFile, "\n");
            }
        }
    }
    fclose(tripFile);
    displayMessage(FULL_NOTIFICATIONS, "Wrote ODs.\n");
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


