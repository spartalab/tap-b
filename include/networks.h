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


typedef enum {
    FORWARD,
    REVERSE
} direction_type;

/*
 * arc_type -- struct for storing all information for arcs.  Most of the
 * members are self-explanatory (tail, head, flow, cost, etc.) and store data
 * from TNTP files.  
 *
 * Of note are the the three function pointers at the end -- calculateCost,
 * calculateDer, and calculateInt.  These respectively calculate the cost on a
 * link, the derivative of this cost at the current flow values, and the
 * integral of the cost from 0 up to the current flow values (for use in the
 * Beckmann function).
 *
 * When the TNTP file is read, appropriate functions are chosen based on the
 * BPR function parameters.  Special functions are provided in case the BPR
 * function is linear (beta = 1) or quartic (beta = 4) which are considerably
 * faster than the generic functions used as a fallback -- the pow operation is
 * expensive, and is best avoided if possible.
 */
typedef struct arc_type {
    long    tail;
    long    head;
    double  flow;
    double  cost;
    double  der;
    double  freeFlowTime;
    double  capacity;
    double  alpha;
    double  beta;
    double  toll;
    double  speedLimit;
    double  length;
    double  fixedCost; /* Reflects toll and distance */
    int        linkType;
    double  (*calculateCost)(struct arc_type *arc);
    double  (*calculateDer)(struct arc_type *der);
    double  (*calculateInt)(struct arc_type *der);
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
    long       numPaths;
} pathSet;

/* node_type -- data structure for nodes.  Only contains lists of arcs entering
 * and leaving the node */
typedef struct {
    arcList forwardStar;
    arcList reverseStar;
} node_type;

/* od_type -- data structure for OD pairs, containing the value in the OD
 * matrix and the min-cost travel time for each OD pair.  (To save memory, the
 * latter can be deleted -- any place that refers to this value [like gap
 * calculations] can recalculate shortest paths on the fly.)
 */
typedef struct {
    double demand;
    double cost;
} od_type;

/* network_type -- data structure for the entire network, including arrays of
 * nodes, arcs, and OD pairs, and network size information.  The beckmann and
 * beckmannLB members are used in certain gap calculations. */
typedef struct {
    node_type*  nodes;
    arc_type*   arcs;
    pthread_mutex_t* arc_muts;
    od_type**   OD;
    long    numNodes;
    long    numArcs;
    long    numZones;
    long    firstThroughNode;
    double  totalODFlow;
    double  tollFactor;
    double  distanceFactor;
    double  beckmann;
    double  beckmannLB;
} network_type;

void BellmanFord(long origin, double *label, long *backnode,
                 network_type *network, queueDiscipline q);
void arcBellmanFord(long origin, double *label, arc_type **backarc,
                    network_type *network, queueDiscipline q);
void arcIndexBellmanFord(long origin, double *label, int *backarc,
                         network_type *network, queueDiscipline q);
void heapDijkstra(long origin, double *label, long *backnode,
                  network_type *network);
void finalizeNetwork(network_type *network);
void quicksortDestinations(long *nodes, double *costs, int elements);
void search(long origin, long* order, long *backnode, network_type *network,
            queueDiscipline q, direction_type d);
void topologicalOrder(network_type *network, long* sequence, long* order);
void findPrimaryLink(network_type *network, arc_type *arc);
int forwardStarOrder(const void *arc1, const void *arc2);
int ptr2arc(network_type *network, arc_type *arcptr);


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

