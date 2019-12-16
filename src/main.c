#include "main.h"

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_NOTIFICATIONS;
    #ifdef DEBUG_MODE
        debugFile = openFile("full_log.txt", "w");
    #endif
 
    /* Uncomment one of these, depending on whether you want to read a network
       in the TNTP format or the NCTCOG network specifically */
    main_TNTP(argc, argv);
    //main_NCTCOG(argc, argv);                    

    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif

    return EXIT_SUCCESS;
}


void main_TNTP(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();
    readOBANetwork(network, argv[1], argv + 2, argc - 2, &Bparameters);    

    if (argc < 3) 
         fatalError("Must specify at least two arguments\n\nUsage:\n "
                    "networkfile demandfile    for 1-class TNTP files\n"
                    "networkfile demandfiles... for multiclass TNTP\n");
                        
    /* Default: one batch */
    setBatches(network, network->numOrigins);
    
    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-14; 
    Bparameters.maxIterations = 99; 
    Bparameters.maxTime = 3600;
    Bparameters.warmStart = FALSE;
    Bparameters.gapFunction = RELATIVE_GAP_1;

    AlgorithmB(network, &Bparameters);
    writeNetworkFlows(network, Bparameters.flowsFile);
    deleteNetwork(network); 
}

void main_NCTCOG(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();
    
    if (argc != 3)
        fatalError("Must specify three arguments for NCTCOG:\n\n"
                   "networkfile triptable convertertable\n\n"
                   "-network file has link data\n"
                   "-triptable has the OD matrix in CSV form\n"
                   "-convertertable translates wrap IDs to TAP-B\n");
        
    
    /* Uncomment the following line to read demand file afresh (rather than
     * from the pre-read binary matrices */
    readNCTCOGNetwork(network, argv[1], argv[2], argv[3]);
    
    /* Uncomment the following line to read archived binary OD matrices */
    /* readNCTCOGNetwork(network, argv[1], NULL, argv[3]); */
    
    /* Default batching for NCTCOG: one per *class* */
    setBatches(network, network->numZones);        
    
    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-14; 
    Bparameters.maxIterations = 5; 
    Bparameters.maxTime = 3600;
    
    Bparameters.warmStart = FALSE;
    Bparameters.calculateBeckmann = FALSE; /* Expensive with conic functions */
    Bparameters.gapFunction = RELATIVE_GAP_1;
    
    AlgorithmB(network, &Bparameters);
    writeNetworkFlows(network, Bparameters.flowsFile);
    deleteNetwork(network);     
}

/*
setBatches: Re-partitions the network into origin batches of the given
size.  This should NEVER be called during the middle of a run, or
unpredictable things may happen. 

As a result, we write the binary matrices HERE.  
*/
void setBatches(network_type *network, int batchSize) {
    network->batchSize = batchSize;
    network->numBatches = (network->numOrigins - 1) / batchSize + 1;
    network->curBatch = 0;

    writeBinaryMatrices(network);    
}

/*
 * CURRENTLY NOT USED since we can now read full NCTCOG.
 *
 * Full version of NCTCOG network not currently available.  So make it from
 * a previously-aggregated version:
 *  1. Multiply # origins by 10 (for 10 classes)
 *  2. Split each OD pair into 10 even parts
 *  3. Assign each class a particular VOT
 *  4. Assign each link a toll equal to its free-flow time.
 *
 * We assume that the NCTCOG network was read in SDB format, not TNTP.
 */
#define NCTCOG_CLASSES 10
#define MAX_VOT 1.0
void inferNCTCOGNetwork(network_type *network) {
    int ij, r, s, c, origin;
    double **newDemand;

    /* 1. Create more classes */
    network->numOrigins = network->numZones * NCTCOG_CLASSES;
    network->numClasses = NCTCOG_CLASSES;
    newDemand = newMatrix(network->numOrigins, network->numZones, double);

    /* 2. Split OD matrix */
    for (r = 0; r < network->numZones; r++) {
        for (s = 0; s < network->numZones; s++) {
            for (c = 0; c < NCTCOG_CLASSES; c++) {
                origin = nodeclass2origin(network, r, c);
                newDemand[origin][s] = network->demand[r][s] / NCTCOG_CLASSES;
            }
        }
    }
    free(network->demand);
    network->demand = newDemand;
    
    /* 3. Set class VOTs */
    network->tollFactor = newVector(NCTCOG_CLASSES, double);
    network->distanceFactor = newVector(NCTCOG_CLASSES, double);
    for (c = 0; c < network->numClasses; c++) {
        network->tollFactor[c] = MAX_VOT * c / NCTCOG_CLASSES;
        network->distanceFactor[c] = 0;
    }

    /* 4.Set link toll to free-flow time and do other needful link things */
    for (ij = 0; ij < network->numArcs; ij++) {
        network->arcs[ij].toll = network->arcs[ij].freeFlowTime;
        network->arcs[ij].length = 0;
        free(network->arcs[ij].classFlow);
        network->arcs[ij].classFlow = newVector(NCTCOG_CLASSES, double);
        for (c = 0; c < network->numClasses; c++) {
            network->arcs[ij].classFlow[c] = 0;
        }
    }


}
