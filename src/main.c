#include "main.h"
#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

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

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_NOTIFICATIONS;
    #ifdef DEBUG_MODE
        debugFile = openFile("full_log.txt", "w");
    #endif
    /* Uncomment one of these, depending on whether you want to read a network
       in the TNTP format or the NCTCOG network specifically */
    main_NCTCOG(argc, argv);
    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif

    return EXIT_SUCCESS;
}

void main_FWtest(int argc, char* argv[]) {
    verbosity = DEBUG;
    
    network_type *network = newScalar(network_type);
    CCparameters_type parameters = initializeCCparameters(CONJUGATE_FRANK_WOLFE);

#if PARALLELISM
    int numOfThreads = 0;
   if (argc != 4) {
       displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will define the num of threads based on the number of available cores.\n");
       displayMessage(FULL_NOTIFICATIONS, "Number of available cores: %d\n", getNumCores());
       numOfThreads = getNumCores();
   } else {
       numOfThreads = atoi(argv[argc - 1]);
   }
#else
    if (argc != 3)
        fatalError("Must specify two arguments\n\nUsage: tap "
                   "networkfile demandfile\n");

#endif

#if PARALLELISM
    parameters.numThreads = numOfThreads;
   if(parameters.numThreads < 1 || parameters.numThreads > 64) {
       fatalError("Invalid number of threads: %d must be between 1 and 64", parameters.numThreads);
   }
   if (argc != 4) {
       readOBANetwork(network, argv[1], argv + 2, argc - 2,
                      parameters.demandMultiplier);
   } else {
       readOBANetwork(network, argv[1], argv + 2, argc - 3,
                      parameters.demandMultiplier);
   }
#else
    readOBANetwork(network, argv[1], argv + 2, argc - 2,
                   parameters.demandMultiplier);
#endif
    parameters.convergenceGap = 1e-4;
    parameters.maxLineSearchIterations = 1;
    displayMessage(FULL_NOTIFICATIONS, "Starting Frank-Wolfe Variation Algorithm...\n");

    convexCombinations(network, &parameters);

    writeNetworkFlows(network, parameters.flowsFile);
    deleteNetwork(network);
}

void main_TNTP(int argc, char* argv[]) {
   network_type *network = newScalar(network_type);
   algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();


#if PARALLELISM
   int numOfThreads = 0;
  if (argc < 5)
    fatalError("Must specify at least four arguments\n\nUsage: tap gap num_classes "
               "networkfile demandfiles [num_threads]\n");
   if (argc == atoi(argv[2]) + 4) {
       displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will define the num of threads based on the number of available cores.\n");
       displayMessage(FULL_NOTIFICATIONS, "Number of available cores: %d\n", getNumCores());
       numOfThreads = getNumCores();
     if (argc-4 != atoi(argv[2]))
        fatalError("Number of classes must match number of input demand files");
   } else {
      if (argc-5 != atoi(argv[2]))
        fatalError("Number of classes must match number of input demand files");
      numOfThreads = atoi(argv[argc-1]);
   }
   Bparameters.numThreads = numOfThreads;
   if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
       fatalError("Invalid number of threads: %d must be between 1 and 64", Bparameters.numThreads);
   }
   else {
      displayMessage(FULL_NOTIFICATIONS, "Number of threads: %d\n", numOfThreads);
   }
#else
   if (argc < 5)
        fatalError("Must specify at least four arguments\n\nUsage: tap gap num_classes "
                   "networkfile demandfile\n");
   if (argc-4 != atoi(argv[2]))
        fatalError("Number of classes must match number of input demand files");
#endif

   readOBANetwork(network, argv[3], argv + 4, atoi(argv[2]), Bparameters.demandMultiplier);
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
   deleteNetwork(network);
}

void main_NCTCOG(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();

#if PARALLELISM
   int numOfThreads = 0;
   if (argc != 5) {
       displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s\n", argv[1], argv[2], argv[3]);
       displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will define the num of threads based on the number of available cores.\n");
       displayMessage(FULL_NOTIFICATIONS, "Number of available cores: %d\n", getNumCores());
       numOfThreads = getNumCores();
   } else {
       displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s, arg4: %s\n", argv[1], argv[2], argv[3], argv[4]);
       numOfThreads = atoi(argv[argc - 1]);
   }

#if PARALLELISM
    Bparameters.numThreads = numOfThreads;
   if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
       fatalError("Invalid number of threads: %d must be between 1 and 64", Bparameters.numThreads);
   }
    /* Uncomment the following line to read demand file afresh (rather than
     * from the pre-read binary matrices */
    displayMessage(FULL_NOTIFICATIONS, "Reading NCTCOG Network...%s\n", argv[1]);
    if (strcmp("", argv[2]) == 0) {
      argv[2] = NULL;
      displayMessage(FULL_NOTIFICATIONS, "No Input Matrix path provided, assuming warm start\n");

    }
    readNCTCOGNetwork(network, argv[1], argv[2], argv[3]);
        
#else
    /* Uncomment the following line to read demand file afresh (rather than
     * from the pre-read binary matrices */
    displayMessage(FULL_NOTIFICATIONS, "Reading NCTCOG Network...%s\n", argv[1]);
    if (strcmp("", argv[2]) == 0) {
      argv[2] = NULL;
      displayMessage(FULL_NOTIFICATIONS, "No Input Matrix path provided, assuming warm start\n");

    }
    readNCTCOGNetwork(network, argv[1], argv[2], argv[3]);
        
#endif

    /* Default batching for NCTCOG: one per *class* */

    setBatches(network, network->numZones, FALSE);
    

    displayMessage(FULL_NOTIFICATIONS, "Total demand: %f\n", network->totalODFlow);
    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-3;
    Bparameters.maxIterations = 5000;
    Bparameters.maxTime = 3600 * 24 * 7;

    Bparameters.warmStart = TRUE;
    Bparameters.reuseFirstBush = TRUE;    
    Bparameters.calculateBeckmann = FALSE; /* Expensive with conic functions */
    Bparameters.gapFunction = RELATIVE_GAP_1;

    AlgorithmB(network, &Bparameters);
    writeNCTCOGFlows(network, Bparameters.flowsFile);
    deleteNetwork(network);
}
void main_NCTCOGFW(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    CCparameters_type parameters = initializeCCparameters(BICONJUGATE_FRANK_WOLFE);
    verbosity = DEBUG;
    
#if PARALLELISM
   int numOfThreads = 0;
   if (argc != 5) {
       displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s\n", argv[1], argv[2], argv[3]);
       displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will define the num of threads based on the number of available cores.\n");
       displayMessage(FULL_NOTIFICATIONS, "Number of available cores: %d\n", getNumCores());
       numOfThreads = getNumCores();
   } else {
       displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s, arg4: %s\n", argv[1], argv[2], argv[3], argv[4]);
       numOfThreads = atoi(argv[argc - 1]);
   }
#else
    displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s\n", argv[1], argv[2], argv[3]);
    if (argc != 4)
        fatalError("Must specify three arguments for NCTCOG:\n\n"
                       "networkfile triptable convertertable\n\n"
                       "-network file has link data\n"
                       "-triptable has the OD matrix in CSV form\n"
                       "-convertertable translates wrap IDs to TAP-B\n");

#endif

#if PARALLELISM
    parameters.numThreads = numOfThreads;
   if(parameters.numThreads < 1 || parameters.numThreads > 64) {
       fatalError("Invalid number of threads: %d must be between 1 and 64", parameters.numThreads);
   }
    /* Uncomment the following line to read demand file afresh (rather than
     * from the pre-read binary matrices */
    displayMessage(FULL_NOTIFICATIONS, "Reading NCTCOG Network...%s\n", argv[2]);
    if (strcmp("", argv[2]) == 0) {
      argv[2] = NULL;
      displayMessage(FULL_NOTIFICATIONS, "Here\n");

    }
    readNCTCOGNetwork(network, argv[1], argv[2], argv[3]);

    /* Uncomment the following line to read archived binary OD matrices */
    /* readNCTCOGNetwork(network, argv[1], NULL, argv[3]); */

#else
    /* Uncomment the following line to read demand file afresh (rather than
     * from the pre-read binary matrices */
    displayMessage(FULL_NOTIFICATIONS, "Reading NCTCOG Network...\n");
    if (strcmp("", argv[2]) == 0) {
      argv[2] = NULL;
      displayMessage(FULL_NOTIFICATIONS, "Here\n");

    }
    readNCTCOGNetwork(network, argv[1], argv[2], argv[3]);

    /* Uncomment the following line to read archived binary OD matrices */
    /* readNCTCOGNetwork(network, argv[1], NULL, argv[3]); */
#endif

    /* Default batching for NCTCOG: one per *class* */

//    setBatches(network, network->numZones, argv[2] == NULL);


    displayMessage(FULL_NOTIFICATIONS, "Starting Frank-Wolfe Variation Algorithm...\n");
    parameters.convergenceGap = 1e-4;
    parameters.maxIterations = 5000;
    parameters.maxTime = 3600 * 24 * 7;

    parameters.warmStart = FALSE;
    parameters.calculateBeckmann = FALSE; /* Expensive with conic functions */

    parameters.maxLineSearchIterations = 1;

    convexCombinations(network, &parameters);

    writeNetworkFlows(network, parameters.flowsFile);
    deleteNetwork(network);
}
#endif

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

    if(network-> numBatches > 1) {
        deleteMatrix(network->demand, network->numOrigins);
        network->demand = newMatrix(network->batchSize, network->numZones, double);
    }

}
