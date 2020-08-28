#include "main.h"
#include "stdlib.h"
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

#if NCTCOG_ENABLED
//    main_NCTCOG(argc, argv);
    // main_NCTCOGFW(argc, argv);
    main_NCTCOGFW_conic_calculator(argc, argv);
#else
     main_TNTP(argc, argv);
//    main_FWtest(argc, argv);
#endif


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

#if NCTCOG_ENABLED
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
    Bparameters.numThreads = numOfThreads;
   if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
       fatalError("Invalid number of threads: %d must be between 1 and 64", Bparameters.numThreads);
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

    setBatches(network, network->numZones, argv[2] == NULL);


    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-8;
    Bparameters.maxIterations = 5000;
    Bparameters.maxTime = 3600 * 24 * 7;

    Bparameters.warmStart = FALSE;
    Bparameters.calculateBeckmann = FALSE; /* Expensive with conic functions */
    Bparameters.gapFunction = RELATIVE_GAP_1;

    AlgorithmB(network, &Bparameters);
    writeNetworkFlows(network, Bparameters.flowsFile);
    deleteNetwork(network);
}

void main_NCTCOGFW(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    CCparameters_type parameters = initializeCCparameters(CONJUGATE_FRANK_WOLFE);

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

struct NCTCOG_Translator {
    int AB_index;
    int BA_index;
};

void main_NCTCOGFW_conic_calculator(int argc, char* argv[]) {
    network_type *network = newScalar(network_type);
    CCparameters_type parameters = initializeCCparameters(CONJUGATE_FRANK_WOLFE);

    displayMessage(FULL_NOTIFICATIONS, "arg1: %s, arg2: %s, arg3: %s\n", argv[1], argv[2], argv[3]);
    if (argc != 5)
        fatalError("Must specify three arguments for NCTCOG:\n\n"
                       "networkfile flowsfile convertertable\n\n"
                       "-networkfile has link data\n"
                       "-flowsfile has flows to calculate times for\n"
                       "-convertertable translates wrap IDs to TAP-B\n"
                       "-translator translates nctcog link index to TAP-B indices");

    displayMessage(FULL_NOTIFICATIONS, "Reading NCTCOG Network...\n");
    readNCTCOGNetwork(network, argv[1], NULL, argv[3]);


    struct NCTCOG_Translator* translator = newVector(89421, struct NCTCOG_Translator);
    for(int i = 0; i < 89421; ++i) {
        translator[i].AB_index = -1;
        translator[i].BA_index = -1;
    }

    char* fname = "nctcog/translator.csv";
    char lineData[5][STRING_SIZE];
    char fullLine[STRING_SIZE];
    FILE *translator_file = openFile(argv[4], "r");
    if (fgets(fullLine, STRING_SIZE, translator_file) == NULL)
        fatalError("Cannot read header of translator_file %s", translator_file);
    while (TRUE) { /* Break in middle when out of lines */
        if (fgets(fullLine, STRING_SIZE, translator_file) == NULL) break;
        if (strstr(fullLine, ",") == NULL) continue;
        parseCSV(lineData, fullLine, 3);
        int idx = atoi(lineData[1]);
        if (strcmp(lineData[2], "True") == 0) {
            translator[idx].AB_index = atoi(lineData[0]);
        } else {
            translator[idx].BA_index = atoi(lineData[0]);
        }
    }
    fclose(translator_file);

    double* times = newVector(network->numArcs, double);
    for (int ij = 0; ij < network->numArcs; ij++) {
        times[ij] = 0.0;
    }
    FILE *flows_file = openFile(argv[2], "r");
    if (fgets(fullLine, STRING_SIZE, flows_file) == NULL)
        fatalError("Cannot read header of flows_file %s", flows_file);
    while (TRUE) { /* Break in middle when out of lines */
        if (fgets(fullLine, STRING_SIZE, flows_file) == NULL) break;
        if (strstr(fullLine, ",") == NULL) continue;
        parseCSV(lineData, fullLine, 5);
        int idx = atoi(lineData[0]);
        int ab_arc_idx = translator[idx].AB_index;
        int ba_arc_idx = translator[idx].BA_index;
        if (ab_arc_idx != -1) {
            network->arcs[ab_arc_idx].flow = atof(lineData[3]);
            times[ab_arc_idx] = conicCost(&network->arcs[ab_arc_idx]);
        }
        if (ba_arc_idx != -1) {
            network->arcs[ba_arc_idx].flow = atof(lineData[4]);
            times[ba_arc_idx] = conicCost(&network->arcs[ba_arc_idx]);
        }
    }
    fclose(flows_file);

    FILE *outFile = openFile("calculated_file.csv", "w");
    int ij;

    for (ij = 0; ij < network->numArcs; ij++) {
        if (network->arcs[ij].capacity == ARTIFICIAL) continue;
        fprintf(outFile, "%f\n", times[ij]);
    }

    fclose(outFile);
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
#if NCTCOG_ENABLED
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
#endif
