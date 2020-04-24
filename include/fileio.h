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

typedef enum {
    ID,
    FROM_NODE,
    TO_NODE,
    PMCAP_AB,
    PMCAP_BA,
    LENGTH,
    PKFRTIME_AB,
    PKFRTIME_BA,
    A_PK_CONICAL,
    VDF_SHIFT,
    SPAR_AB,
    SPAR_BA,
    PMSATFLOW_AB,
    PMSATFLOW_BA,
    CA_AB,
    CB_AB,
    CC_AB,
    CD_AB,
    CA_BA,
    CB_BA,
    CC_BA,
    CD_BA,
    UNSIG_MINDELAY,
    UPAR_AB,
    UPAR_BA,
    OPERCOSTPM_AB,
    OPERCOSTPM_BA,
    TOLL_AUTO_DA_PM_AB,
    TOLL_AUTO_DA_PM_BA,
    TOLL_AUTO_SR_PM_AB,
    TOLL_AUTO_SR_PM_BA,
    TOLL_MEDIUM_TRUCK_PM_AB,
    TOLL_MEDIUM_TRUCK_PM_BA,
    TOLL_HEAVY_TRUCK_PM_AB,
    TOLL_HEAVY_TRUCK_PM_BA,
    FLAG_EXCLUDE_DA,
    NUM_NCTCOG_NET_COLUMNS
} NCTCOG_net_columns;

void readNCTCOGNetwork(network_type *network, char *networkFileName,
                       char *tripFileName, char *converterFileName);
void writeBinaryMatrices(network_type *network);
void readConverterFile(char *converterFileName, int *table, int maxValue,
                       int maxKey, bool isInverse);
int convert(int value, int *table, int maxKey);
void readNCTCOGLinks(network_type *network, char *networkFileName, int *table);
void makeLink(network_type *network, int ij, int *table, char *ID, char *from,
        char *to, char *cap, char *len, char *freeFlow, char *conical, char
        *shift, char *sPar, char *satFlow, char *CA ,char *CB, char *CC, char
        *CD, char *minDelay, char *uPar, char *operCost, char *daToll, char
        *srToll, char *medToll, char *heavyToll, char *exclude);
void readNCTCOGTrips(network_type *network, char *tripFileName, int *table);
void streamNCTCOGTrips(network_type *network, int *table);
void assignDemand(network_type *network, int originNode, int destinationNode,
                  int class, char *demandValue);
void assignStreamedDemand(network_type *network, int originNode, int destinationNode,
                          int class, float demandValue);
void readOBANetwork(network_type *network, char *linkFileName,
                    char **tripFileName, int numClasses,
                    algorithmBParameters_type *parameters);
void readSDBNetwork(network_type *network, char *filename, double defaultAlpha,
                    double defaultBeta);
void writeDUENetwork(network_type *Network, char *filename);
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
