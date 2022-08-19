#include "main.h"

#ifdef PARALLELISM
    #ifdef _WIN32
        #include <windows.h>
    #elif MACOS
        #include <sys/param.h>
        #include <sys/sysctl.h>
    #else
        #include <unistd.h>
    #endif
#endif

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_NOTIFICATIONS;
#ifdef DEBUG_MODE
    debugFile = openFile("full_log.txt", "w");
#endif

    network_type *network = newScalar(network_type);
    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();

#ifdef PARALLELISM
    int numOfThreads = 0;
#endif

    if (argc < 5)
        fatalError("Must specify at least four arguments\n\nUsage: tap gap "
                   "num_classes networkfile demandfiles [num_threads]\n");
    if (argc == atoi(argv[2]) + 4) {
#ifdef PARALLELISM
        displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will "
                       "define the num of threads based on the number of "
                       "available cores.\n");
        displayMessage(FULL_NOTIFICATIONS, "Number of available cores: %d\n",
                       getNumCores());
        numOfThreads = getNumCores();
#endif
    } else {
        if (argc != atoi(argv[2]) + 5)
            fatalError("Number of classes and input demand files must match.");
#ifdef PARALLELISM
        numOfThreads = atoi(argv[argc-1]);
#endif
    }
#ifdef PARALLELISM
    Bparameters.numThreads = numOfThreads;
    if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
        fatalError("Invalid number of threads: %d must be between 1 and 64",
                   Bparameters.numThreads);
    } else {
      displayMessage(FULL_NOTIFICATIONS,"Number of threads: %d\n",numOfThreads);
    }
#endif
    readOBANetwork(network, argv[3], argv + 4, atoi(argv[2]),
                   Bparameters.demandMultiplier);
    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = atof(argv[1]);
    Bparameters.maxIterations = 200;
    Bparameters.maxTime = 10000;
    // Bparameters.storeBushes = TRUE; // Uncomment if you want to save bushes for future warm start use
    Bparameters.warmStart = FALSE; //Set to true if you want to warm start. Batch size must be set to the size used when first storing the bush
    Bparameters.gapFunction = RELATIVE_GAP_1;
    Bparameters.calculateBeckmann = TRUE; /* Expensive with conic functions */

    /* Default: one batch */
    setBatches(network, network->numOrigins, atoi(argv[2]) == 0);

    AlgorithmB(network, &Bparameters);
    writeNetworkFlows(network, Bparameters.flowsFile);
    if (network->numClasses > 1) {
        for (int c = 0; c < network->numClasses; c++) {
           displayMessage(FULL_NOTIFICATIONS,
                          "Class %d TSTT: %f\n", c+1, classGeneralizedCost(network, c));
        }
    }
    displayMessage(FULL_NOTIFICATIONS, "Aggregate TSTT: %f\n", TSTT(network));
    deleteNetwork(network);

#ifdef DEBUG_MODE
    fclose(debugFile);
#endif

    return EXIT_SUCCESS;
}


/*
setBatches: Re-partitions the network into origin batches of the given
size.  This should NEVER be called during the middle of a run, or
unpredictable things may happen.

As a result, we write the binary matrices HERE.
*/
void setBatches(network_type *network, int batchSize, bool warmStart) {
    network->batchSize = batchSize;
    if (network->numOrigins%batchSize != 0) {
        fatalError("Number of Origins (%d) must be divisible by the batch size(%d)\n", network->numOrigins, batchSize);
    }
    network->numBatches = (network->numOrigins - 1) / batchSize + 1;
    network->curBatch = 0;

    if(!warmStart)
      writeBinaryMatrices(network);

    if(network-> numBatches > 1 && !warmStart) {

        deleteMatrix(network->demand, network->numOrigins);
        network->demand = newMatrix(network->batchSize, network->numZones, double);
    }

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
