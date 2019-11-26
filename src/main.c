#include "main.h"

int main(int argc, char* argv[]) {
    /* verbosity is a global variable controlling how much output to produce,
     * see utils.h for possible values*/
    verbosity = FULL_NOTIFICATIONS; 
    #ifdef DEBUG_MODE
        debugFile = openFile("full_log.txt", "w");
    #endif

    if (argc != 4)
         fatalError("Must specify three  arguments\n\nUsage: tap "
                    "networkfile demandfile num_threads\n");

    network_type *network = newScalar(network_type);

    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();

    readOBANetwork(network, argv[1], argv[2]);

    makeStronglyConnectedNetwork(network); /* Check connectivity */

    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-14;
    Bparameters.maxIterations = 1000;
    Bparameters.maxTime = 3600;
    Bparameters.gapFunction = RELATIVE_GAP_1;
    Bparameters.numThreads = atoi(argv[3]);
    
    AlgorithmB(network, &Bparameters);

    deleteNetwork(network);

    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif

    return EXIT_SUCCESS;
}
