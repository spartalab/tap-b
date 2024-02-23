/*
 * fileio.h -- This is the header for file reading and writing.  String
 * processing routines also go here.
 */

#ifndef FILEIO_H
#define FILEIO_H

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bush.h"
#include "networks.h"
#include "tap.h"
#include "utils.h"

#define STRING_SIZE 9999
#define DESTINATIONS_PER_LINE 5

#define PAUSE_ON_ERROR FALSE
#define PAUSE_ON_WARNING FALSE

#ifdef DEBUG_MODE
char debugFileName[STRING_SIZE];
FILE *debugFile;
#endif

int screenVerbosity, logVerbosity;

enum { // Return codes for metadata parsing
    SUCCESS,
    BLANK_LINE,
    COMMENT
};

typedef struct stringList_type { // Stack of strings
    char string[STRING_SIZE];
    struct stringList_type *prev;
} stringList_type;

///////////////////////////
// Reading network files //
///////////////////////////

network_type *readParametersFile(algorithmBParameters_type *thisRun,
                                 char *filename);
void writeBinaryMatrices(network_type *network, char *matrixStem);
void readOBANetwork(network_type *network, char *linkFileName,
                    char **tripFileName, int numClasses,
                    double defaultDemandMultiplier);
void setBatches(network_type *network, int batchSize, bool warmStart,
                bool storeMatrices, char *matrixStem);
void writeOBANetwork(network_type *Network, char *linkfile, char *tripfile);
void writeNetworkFlows(network_type *network, char *outputFileName);
void writeMulticlassFlows(network_type *network, char *outputFileName);

///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length);
int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue);
int parseLine(char* inputLine, char* outputLine);
void parseCSV(char field[][STRING_SIZE], char *fullLine, int numFields);

#ifdef PARALLELISM
int getNumCores();
#endif

#endif
