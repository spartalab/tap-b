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
    verbosity = FULL_DEBUG;
    //verbosity = FULL_NOTIFICATIONS;
    network_type *network;
#ifdef DEBUG_MODE
    debugFile = openFile("full_log.txt", "w");
#endif

    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();
    if (argc != 2)
        fatalError("Must specify exactly one parameter (parameters file).\n");
    network = readParametersFile(&Bparameters, argv[1]);
    AlgorithmB(network, &Bparameters);
    writeNetworkFlows(network, Bparameters.flowsFile);
    if (network->numClasses > 1) {
        for (int c = 0; c < network->numClasses; c++) {
           displayMessage(FULL_NOTIFICATIONS,
                          "Class %d TSTT: %f\n", c+1,
                          classGeneralizedCost(network, c));
        }
    }
    displayMessage(FULL_NOTIFICATIONS, "Aggregate TSTT: %f\n", TSTT(network));
    deleteNetwork(network);

#ifdef DEBUG_MODE
    fclose(debugFile);
#endif

    return EXIT_SUCCESS;
}
