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


///////////////////////////
// Reading network files //
///////////////////////////

void writeBinaryMatrices(network_type *network);
void readOBANetwork(network_type *network, char *linkFileName,
                    char **tripFileName, int numClasses,
                    double defaultDemandMultiplier);
void writeOBANetwork(network_type *Network, char *linkfile, char *tripfile);
void writeNetworkFlows(network_type *network, char *outputFileName);

///////////////////////
// String processing //
///////////////////////

void blankInputString(char *string, int length);
int parseMetadata(char* inputLine, char* metadataTag, char* metadataValue);
int parseLine(char* inputLine, char* outputLine);
void parseCSV(char field[][STRING_SIZE], char *fullLine, int numFields);

#endif
