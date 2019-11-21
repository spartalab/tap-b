#include "main.h"

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_NOTIFICATIONS; 
    #ifdef DEBUG_MODE
        debugFile = openFile("full_log.txt", "w");
    #endif

#if PARALLELISM
    if (argc != 4)
        fatalError("Must specify three arguments\n\nUsage: tap "
                   "networkfile demandfile numthreads\n");
#else
    if (argc != 3)
         fatalError("Must specify two arguments\n\nUsage: tap "
                    "networkfile demandfile\n");

#endif
    network_type *network = newScalar(network_type);

    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();

    readOBANetwork(network, argv[1], argv[2]);

    makeStronglyConnectedNetwork(network); /* Check connectivity */

    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-14;
    Bparameters.maxIterations = 1000;
    Bparameters.maxTime = 3000;
    Bparameters.gapFunction = RELATIVE_GAP_1;
#if PARALLELISM
    Bparameters.numThreads = atoi(argv[3]);
    if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
        fatalError("Invalid number of threads: %d must be between 1 and 64", Bparameters.numThreads);
    }
#endif
    AlgorithmB(network, &Bparameters);

    deleteNetwork(network);

    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif

    return EXIT_SUCCESS;
}
