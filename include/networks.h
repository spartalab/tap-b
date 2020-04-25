/*
 * networks.h -- header file for general network management.  Defines the data
 * structures for storing networks and their components (arcs, nodes, OD
 * pairs).  Also contains implementations of standard network algorithms
 * (shortest path, topological order) and routines for managing the network
 * data structures.
 */

#ifndef NETWORKS_H
#define NETWORKS_H

#include <math.h>
#include <string.h>
#include "datastructures.h"
#include <pthread.h>

#define NO_PATH_EXISTS -1
#define ARTIFICIAL 99999 /* Value used for costs, etc. on artificial links
                            generated to ensure strong connectivity */
#define NOT_ALLOWED ARTIFICIAL /* Use same value if a class cannot use links*/

typedef enum {
    SOLO_17,
    SOLO_35,
    SOLO_45,
    SOLO_90,
    HOV_17,
    HOV_35,
    HOV_45,
    HOV_90,
    MED_TRUCKS,
    HVY_TRUCKS,
    NUM_NCTCOG_CLASSES
} NCTCOG_classes;

#define isHOV(c) ((c)==HOV_17 || (c)==HOV_35 || (c)==HOV_45 || (c)==HOV_90)
#define isSolo(c) (!isHOV(c))
#define isDA(c) ((c)==SOLO_17 || (c)==SOLO_35 || (c)==SOLO_45 || (c)==SOLO_90)


typedef enum {
    FORWARD,
    REVERSE
} direction_type;

/*
 * arc_type -- struct for storing all information for arcs.  Most of the
 * members are self-explanatory (tail, head, flow, cost, etc.) and store data
 * from TNTP files.  
 *
 * NOTE: This function is specialized now for the NCTCOG network.
 *
 * Of note are the the three function pointers at the end -- calculateCost,
 * calculateDer, and calculateInt.  These respectively calculate the cost on a
 * link, the derivative of this cost at the current flow values, and the
 * integral of the cost from 0 up to the current flow values (for use in the
 * Beckmann function).
 *
 */
typedef struct arc_type {
    int    tail;
    int    head;
    double  flow;
    double* classFlow; /* [class] */
    double* classCost; /* [class]... fixed costs per class */
    double* classToll; /* [class]... separate from classCost to get revenue  */
    double  cost;
    double  der;

    int     ID; /* NCTCOG ID */
    double  freeFlowTime;
    double  capacity;
    double  a; /* Conical parameter */
    double  e; /* VDF shift */
    double  sParam; /* Signal parameter */
    double  saturationFlow;
    double  CA;
    double  CB;
    double  CC;
    double  CD; /* Not sure what these are */
    double  m; /* Unsignalized minimum delay */
    double  u; /* Unsignalized parameter */
    double  operCost; /* Operating cost */
    double  toll_DA; /* Tolls for drive-alone, shared-ride, medium, heavy */
    double  toll_SR;
    double  toll_Med;
    double  toll_Heavy;
    double  length;
    bool    excludeDA; /* Are drive-alone prohibited? */

    double  b; /* Experimental constants to reduce run time */
    double  h0;
    double oldRoot;

    double  alpha; /* Vestigial BPR values */
    double  beta;
    double  toll;
    double  speedLimit;
    int     linkType;

    double  fixedCost; /* Reflects toll and distance */
    double  (*calculateCost)(struct arc_type *arc);
    double  (*calculateDer)(struct arc_type *der);
    double  (*calculateInt)(struct arc_type *der, bool);
} arc_type;


/* Data structures for linked lists of arcs (used for forward/reverse stars) */
#define arcListElt struct AL
arcListElt {
    arc_type    *arc;
    arcListElt  *prev;
    arcListElt  *next;
};

typedef struct {
    arcListElt  *head;
    arcListElt  *tail;
    int         size;
} arcList;

typedef struct {
    arcList *arcs;
    double  cost;
    double  der;
} path_type;

/* Linked lists of paths (not used in Algorithm B) */
#define pathSetElt struct PSE
pathSetElt {
    path_type   *path;
    pathSetElt  *prev;
    pathSetElt  *next;
};

typedef struct {
    pathSetElt *head;
    pathSetElt *tail;
    int       numPaths;
} pathSet;

/* node_type -- data structure for nodes.  Only contains lists of arcs entering
 * and leaving the node */
typedef struct {
    arcList forwardStar;
    arcList reverseStar;
} node_type;


/* network_type -- data structure for the entire network, including arrays of
 * nodes, arcs, and OD pairs, and network size information.  The beckmann and
 * beckmannLB members are used in certain gap calculations. 
 *
 * Specialized to NCTCOG by not including class-specific VOT/VOD... for that
 * see the individual link classCosts
 * */
typedef struct {
    node_type*  nodes;
    arc_type*   arcs;
    pthread_mutex_t* arc_muts;
    double**    demand;
    int    numNodes;
    int    numArcs;
    int    numZones; /* Number of physical zones */
    int    numOrigins; /* Number of origins (zones * classes) */
    int    numClasses;
    int    firstThroughNode;
    double  totalODFlow;
    double  beckmann;
    double  beckmannLB;
    double *tollFactor; /* per class */
    double *distanceFactor; /* per class */
    int    batchSize; /* Indicates number of batches per origin */
    int    numBatches;
    int    curBatch; /* Indicates which origin batch is current */
} network_type;

/* Used to detect when we are at the end of all batches and must reset */
#define END_OF_ORIGINS_SENTINEL -2


void BellmanFord(int origin, double *label, int *backnode,
                 network_type *network, queueDiscipline q);
void arcBellmanFord(int origin, double *label, arc_type **backarc,
                    network_type *network, queueDiscipline q);
void arcIndexBellmanFord(int origin, double *label, int *backarc,
                         network_type *network, queueDiscipline q);
void BellmanFord_NoLabel(int origin, double *label, network_type *network,
                         queueDiscipline q, double *labelGuess, int *order);
void heapDijkstra(int origin, double *label, int *backnode,
                  network_type *network);


void changeFixedCosts(network_type *network, int class);

void finalizeNetwork(network_type *network);
void quicksortDestinations(int *nodes, double *costs, int elements);

void search(int origin, int* order, int *backnode, network_type *network,
            queueDiscipline q, direction_type d);
void topologicalOrder(network_type *network, int* sequence, int* order);

void findPrimaryLink(network_type *network, arc_type *arc);
int forwardStarOrder(const void *arc1, const void *arc2);
int ptr2arc(network_type *network, arc_type *arcptr);
int origin2node(network_type *network, int origin);
int origin2class(network_type *network, int origin);
int nodeclass2origin(network_type *network, int originNode, int class);
bool checkIfCurrentNodeClass(network_type *network, int originNode, int class);
bool checkIfCurrentOrigin(network_type *network, int origin);
bool outOfOrigins(network_type *network, int origin);

/* Not needed for NCTCOG */
/* void setWeights(network_type *network, int class, double *timeFactor,
                double *tollFactor, double *distanceFactor);
*/

/////////////////////////////////
// Custom linked lists for TAP //
/////////////////////////////////

void displayNetwork(int minVerbosity, network_type *network);
void deleteNetwork(network_type *network);

arcList *createArcList();
void initializeArcList(arcList *list);
arcListElt *insertArcList(arcList *list, arc_type *value, arcListElt *after);
void clearArcList(arcList *list);
void deleteArcList(arcList *list);
void deleteArcListElt(arcList *list, arcListElt *elt);
void displayArcList(arcList *list);
path_type *createPath();
void displayPath(path_type *path);
void displayPathCompact(int minVerbosity, path_type *path);
bool comparePaths(path_type *path1, path_type *path2);

pathSet *createPathSet();
void initializePathSet(pathSet *list);
bool pathsEqual(arcList *path1, arcList *path2);
pathSetElt *insertPathSet(pathSet *list, path_type *value, pathSetElt *after);
void clearPathSet(pathSet *list);
void deletePathSet(pathSet *list);
void deletePathSetElt(pathSet *list, pathSetElt *elt);
void displayPathSet(pathSet *list);

#endif

