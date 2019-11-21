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
    int numOfThreads = 0;
#if PARALLELISM
    if (argc != 4) {
        displayMessage(FULL_NOTIFICATIONS, "Threads were not defined, we will define the num of threads based on the number of available cores.\n");
        displayMessage(FULL_NOTIFICATIONS, "Number of available codes: %d\n", getNumCores);
        numOfThreads = getNumCores();
    } else {
        numOfThreads = atoi(argv[3]);
    }
#else
    if (argc != 3)
         fatalError("Must specify two arguments\n\nUsage: tap "
                    "networkfile demandfile\n");

#endif
    network_type *network = newScalar(network_type);

    algorithmBParameters_type Bparameters = initializeAlgorithmBParameters();

#if PARALLELISM
    Bparameters.numThreads = numOfThreads;
    if(Bparameters.numThreads < 1 || Bparameters.numThreads > 64) {
        fatalError("Invalid number of threads: %d must be between 1 and 64", Bparameters.numThreads);
    }
#endif
    readOBANetwork(network, argv[1], argv[2]);

    makeStronglyConnectedNetwork(network); /* Check connectivity */

    displayMessage(FULL_NOTIFICATIONS, "Starting Algorithm B...\n");
    Bparameters.convergenceGap = 1e-14;
    Bparameters.maxIterations = 1000;
    Bparameters.maxTime = 3000;
    Bparameters.gapFunction = RELATIVE_GAP_1;

    AlgorithmB(network, &Bparameters);

    deleteNetwork(network);

    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif

    return EXIT_SUCCESS;
}
