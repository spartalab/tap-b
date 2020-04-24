/*
 * bush.h -- This is the header file for the bush and merge data structures.
 * These structures are the key for handling large-scale networks efficiently.
 *
 * merge_type is a struct containing extra information for "merge" nodes (that
 * is, those with more than one incoming link in a bush).  Empirically, most
 * nodes in a bush have only a single incoming link.  It is redundant to store
 * flow data for such links, since flow conservation from the outgoing link
 * uniquely determines flows on the incoming link.  It is only at merge nodes,
 * with multiple approaches, that we need to specify how flow is distributed.
 *
 * bushes_type is a struct containing information for ALL network bushes.  Some
 * of its elements store information for each bush separately -- bushOrder
 * (topological order), the bush topology (pred), the last merge node, and the
 * merge data structures.  Other elements are defined only for the current bush
 * being operated on -- longest and shortest path costs, node and link flows.
 * This information is re-created when processing each bush, rather than
 * stored.
 *
 * algorithmBParameters_type is the struct storing algorithm parameters,
 * including convergence criteria, numerical tolerances, flags indicating how
 * often different steps are performed, and pointers to functions for
 * initializing bushes, shifting flows, and finding topological orders.  This
 * code has special implementations for some of these that are faster than
 * naive ones.
 *
 * More details on these data structures are in the comments below.  Quantities
 * in [brackets] in these comments show the interpretation of array indices.
 */


#ifndef BUSH_H
#define BUSH_H

#include <limits.h>
#include <math.h>
#include "tap.h"
#include "networks.h"
#include "datastructures.h"
#include "utils.h"

#define NEW_LINK -1
//#define PARALLELISM 1

/*
 * merge_type: Extra data for merge nodes
 *
 *  numApproaches -- the number of bush links entering the merge node
 *  approach -- an array of link IDs for links entering the merge node
 *  approachFlow -- an array of bush flows for the links entering the merge
 *                  node.
 *  LPlink, SPlink -- these values indicate the *approach IDs* for the longest
 *                    and shortest path backlinks at this node.  Note that this
 *                    is NOT the network link IDs, but the index of the
 *                    approaches for this merge.  LPlink and SPlink thus can be
 *                    directly used as indices for the approach and
 *                    approachFlow arrays.
 *  divergenceNode -- indicates the ID of a node further upstream where flow
 *                    shifts for this merge will start.
 */
typedef struct merge_type {
   int      numApproaches;
   int      *approach; /* [approachIndex] */ 
   double   *approachFlow; /* [approachIndex] */
   int      LPlink;
   int      SPlink;
   int      divergenceNode; 
} merge_type;

/*
 * bushes_type: Stores the data for ALL bushes associated with the network
 *
 * The following members are *shared* across all bushes, and are associated
 * with whatever bush is currently being operated on.  To save memory, this
 * information is overwritten when we move to another bush.
 *  LPcost, SPcost -- arrays of longest used and shortest path costs, indexed
 *                    by node ID
 *  flow -- array of bush flows, indexed by link ID.
 *  nodeFlow -- array of total flow through each node in the bush, indexed by
 *              node ID.
 *
 * The following members are stored *separately* for each bush, and contain
 * information about bushes which is persistent even when other bushes are
 * being operated on.
 *  bushOrder -- stores the inverse topological order for each bush, that is,
 *               bushOrder[origin][index] gives the link ID of the node which
 *               is in the index-th position topologically, for bush 'origin'.
 *  pred -- 2D array storing the bush topology, indexed by bush origin and
 *          the node ID.  To test whether the node is a merge or not, use the
 *          isMergeNode() function which returns TRUE if there is more than one
 *          incoming node in the bush, and FALSE otherwise.  This function
 *          hides the following implementation of bush structure: If bushOrder
 *          is a NON-NEGATIVE value, it indicates that this node has only a
 *          single predecessor link (it is not a merge node), and
 *          pred[origin][node] gives the ID of this link.  If bushOrder is a
 *          NEGATIVE value, it indicates that the node is a merge node, and
 *          -pred[origin][node] gives the index of the corresponding merge_type
 *          struct in the bush merges array.
 *  lastMerge -- array storing topological labels for the merge nodes with
 *               highest topological order.  Many bush processing algorithms
 *               can start or stop at this point, since all nodes with higher
 *               topological order are leaf nodes and may be amenable to
 *               streamlined treatment.
 *  numMerges -- array storing the number of merge nodes associated with each 
 *               bush.
 *  merge -- 2D array storing pointers to the merge_type data structures
 *           associated with each bush.  The first index is to the bush, the
 *           second to the index of a particular merge within the bush (see
 *           numMerges above for the size of this array).
 */
typedef struct bushes_type {
   double   *LPcost; /* [node] */
   double   *SPcost; /* [node] */
   double   *flow; /* [link] */
   double   *nodeFlow; /* [node] */
   bool     *updateBush; /* [origin]... process or skip this origin? */
   int      **bushOrder; /* [origin][nodeOrder] */
   int      **pred; /* [origin][node] */
   int      *lastMerge; /* [origin] */
   int      *numMerges; /* [origin] */
   merge_type ***merges; /* [origin][merge]*/
   double   **LPcost_par; /* [thread][node] */
   double   **SPcost_par; /* [thread][node] */
   double   **flow_par; /* [thread][node] */
   double   **nodeFlow_par; /* [thread][node] */
   double   **nodeFlowshift_par; /* [origin][shift] */


   network_type *network; /* Points back to the corresponding network */
} bushes_type;

/*
 * algorithmBparameters_type: Stores miscellaneous parameters associated with
 * Algorithm B.
 *
 * Termination criteria (BE SURE TO SET AT LEAST ONE OF convergenceGap,
 * maxTime, OR maxIterations, OR ELSE THE ALGORITHM WILL NEVER TERMINATE.):
 *  gapFunction -- indicates what gap function to use (e.g., average excess
 *                 cost, or relative gap).  See tap.h for the possible values
 *                 this can take.  Default value is RELATIVE_GAP_1
 *  convergenceGap -- algorithm terminates after the gap is below this value.
 *                    Default value is 0 (i.e., no termination based on gap.)
 *  maxTime -- algorithm terminates after this many seconds of run time.
 *             Default value is INFINITY (i.e., no termination based on time.)
 *  maxIterations -- algorithm stops after this many iterations.  Default value
 *                   is int_MAX (i.e., no termination based on iterartions.)
 *
 * Algorithm parameters:
 *  innerIterations -- How many times to shift flow on bushes before updating
 *                     the bush topology.  Default value is 20.
 *  shiftReps -- How many times to shift flow on one bush before moving to the
 *               next.  Default value is 1.
 *  rescanAfterShift -- Boolean indicating whether to update shortest paths and
 *                      longest paths after each iteration of flow shifting.
 *                      This will not do anything if shiftReps == 1.  Default
 *                      value is FALSE.
 *  thresholdGap -- Allows you to skip over bushes that are close to
 *                  equilibrium.  If the maximum difference between a longest
 *                  used and shortest path label is less than thresholdGap, no
 *                  flow shifts are performed.  Default value is 0.
 *  minCostDifference -- Allows you to skip over nodes that are close to
 *                       equilibrium.  If the maximum difference between the
 *                       longest used and shortest path to a merge node is
 *                       less than minCostDifference, no flow shift is
 *                       performed.  Default value is 0.
 *  minLinkFlowShift -- Allows you to skip over nodes that have little flow
 *                      passing through them.  If the flow through a node is
 *                      less than minLinkFlowShift, no flow shift is performed.
 *                      Default value is 0.
 *  minLinkFlow -- Parameter to guard against numerical errors in link flows;
 *                 any link flow smaller than this value is assumed to be zero.
 *                 Default value is 1e-14.
 *  minDerivative -- Parameter to avoid dividing by zero in Newton's Method.
 *                   A zero denominator is replaced by minDerivative.
 *                   Default value is 1e-6.
 *  newtonStep -- Step size in Newton's method.  Default value is 1.
 *  numNewtonShifts -- Indicates how many iterations of Newton's method to
 *                     perform per merge node.  Default value is 1.
 *  SPQueueDiscipline -- How to manage the scan eligible list in case the
 *                       label correcting shortest path is run (e.g., during
 *                       initialization).  See datastructures.h for options;
 *                       default value is DEQUE.
 *  createInitialBush -- Function pointer for how bushes are initially set up.
 *                       Default value is initialBushShortestPath (setting bush
 *                       to the one-to-all shortest path tree at free flow.
 *  topologicalOrder -- Function pointer for calculating topological order.  
 *                      Default value is genericTopologicalOrder, the classic
 *                      algorithm for topological ordering.  More specialized
 *                      algorithms are possible (e.g., trying to minimize the
 *                      topological order of the last merge node.  cf.
 *                      lastMerge in bushes_type).
 *  linkShiftB -- Function pointer for adjusting flow on a link.  Default value
 *                is exactCostUpdate, which recalculates the link performance
 *                functions and derivatives whenever flow is shifted.
 *                Alternatives are linearCostUpdate (which approximates the
 *                cost change using the derivative, and leaves the derivative
 *                unchanged) and noCostUpdate (which only updates flow, but not
 *                cost).  If you switch to one of these methods, be sure that
 *                link costs are eventually updated somewhere.
 *
 *  storeMatrices -- boolean.  If TRUE, OD matrices when read are converted
 *                       to binary form (see matrixStem below).  Otherwise,
 *                       the entire OD matrix is retained in memory.  Default
 *                       behavior is FALSE.  Binary matrices are always used if
 *                       if numBatches > 1 (overriding this parameter).
 *  storeBushes -- Force writing of intermediate bushes at each iteration.  
 *                 Default behavior is FALSE.  Intermediate bushes are always
 *                 written if numBatches > 1 (overriding this parameter).
 *  reuseFirstBush -- Faster initialization by copying first bush structure for
 *                    other batches.  This only works if the origins in each
 *                    batch coincide (e.g., if each batch represents one user
 *                    class.  Default behavior is FALSE.  Be careful with using
 *                    this setting, the code does not (yet) check that the
 *                    origins actually do coincide, and strange behavior will
 *                    result if they do not.  
 *
 *
 *  includeGapTime -- include gap calculation time in run times?  Default TRUE
 *
 *  batchStem -- prefix for files storing batches of bushes in binary format.
 *               Default = "batch", so files are batch0.bin, batch1.bin, etc.
 *  matrixStem -- prefix for files storing binary OD matrices for each batch.
 *                Default = "matrix", so files are matrix0.bin, etc.
 *  flowsFile -- name for file to write flows, default is "flows.txt"
 */
typedef struct algorithmBParameters_type{
   gap_type gapFunction;
   double   convergenceGap;
   double   maxTime;
   int      maxIterations;
   int      innerIterations;
   int      shiftReps;
   bool     rescanAfterShift;
   double   demandMultiplier;
   double   thresholdGap;
   double   thresholdAEC;
   double   minCostDifference;
   double   minLinkFlowShift;
   double   minLinkFlow;
   double   minDerivative;
   double   newtonStep;
   int      numNewtonShifts;
   int      numFlowShifts; /* Counter to measure algorithm performance */
   bool     warmStart;
   bool     calculateBeckmann;
   queueDiscipline SPQueueDiscipline;
   void     (*createInitialBush)(int, network_type *, bushes_type *,
                                 struct algorithmBParameters_type *);
   void     (*topologicalOrder)(int, network_type *, bushes_type *,
                                struct algorithmBParameters_type *);
   void     (*linkShiftB)(int, double, network_type *);

#if PARALLELISM
   int numThreads;
#endif
   bool     storeMatrices;
   bool     storeBushes;
   bool     reuseFirstBush;
   bool     includeGapTime;
   char     batchStem[STRING_SIZE-20]; /* Leave space for index */
   char     matrixStem[STRING_SIZE-20];
   char     flowsFile[STRING_SIZE];
} algorithmBParameters_type;

typedef enum scan_type {
    LONGEST_BUSH_PATH,
    LONGEST_USED_PATH,
    LONGEST_USED_OR_SP,
    NO_LONGEST_PATH
} scan_type;

/* Master routine and parameters */
void AlgorithmB(network_type *network, algorithmBParameters_type *parameters);
algorithmBParameters_type initializeAlgorithmBParameters();

/* Main Algorithm B helper functions */
void initializeAlgorithmB(network_type *network, bushes_type **bushes,
                          algorithmBParameters_type *parameters);
void updateBatchBushes(network_type *network, bushes_type *bushes,
                       int *lastClass, algorithmBParameters_type *parameters);
void updateBatchFlows(network_type *network, bushes_type *bushes,
                      int *lastClass, algorithmBParameters_type *parameters);
void loadBatch(int batch, network_type *network, bushes_type **bushes,
               algorithmBParameters_type *parameters);
void storeBatch(int batch, network_type *network, bushes_type *bushes,
               algorithmBParameters_type *parameters);                                         
void initializeBushesB(network_type *network, bushes_type *bushes,
                       struct algorithmBParameters_type *parameters);
void updateBushB(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters);
bool updateFlowsB(int origin, network_type *network, bushes_type *bushes,
                  algorithmBParameters_type *parameters);

/* Custom gap routines using bushes */
double bushSPTT(network_type *network, bushes_type *bushes, 
              algorithmBParameters_type *parameters);
double bushTSTT(network_type *network, bushes_type *bushes);
double bushRelativeGap(network_type *network, bushes_type *bushes,
              algorithmBParameters_type *parameters);
double bushAEC(network_type *network, bushes_type *bushes,
              algorithmBParameters_type *parameters);
double bushMEC(network_type *network, bushes_type *bushes,
              algorithmBParameters_type *parameters);

/* Basic bush manipulations */
bushes_type *createBushes(network_type *network);
void deleteBushes(network_type *network, bushes_type *bushes);
void initialBushShortestPath(int origin, network_type *network,
                             bushes_type *bushes,
                             algorithmBParameters_type *parameters);
void initialBushBFS(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters);
void genericTopologicalOrder(int origin, network_type *network,
                             bushes_type *bushes,
                             algorithmBParameters_type *parameters);
void mergeFirstTopologicalOrder(int origin, network_type *network,
                                bushes_type *bushes,
                                algorithmBParameters_type *parameters);
void scanBushes(int origin, network_type *network, bushes_type *bushes,
                algorithmBParameters_type *parameters, scan_type LPrule);
void reconstructMerges(int origin, network_type *network, bushes_type *bushes);
void findDivergenceNodes(int origin, network_type *network,
                         bushes_type *bushes);
bool rescanAndCheck(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters);

void updateFlowPass(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters);

void calculateBushFlows(int origin, network_type *network,
                        bushes_type *bushes);
void pushBackFlowSimple(int j, int origin, network_type *network,
                        bushes_type *bushes);

void pushBackFlowMerge(merge_type *merge, network_type *network,
                       bushes_type *bushes);

void rectifyMerge(int j, merge_type *merge, bushes_type *bushes);
void rectifyBushFlows(int origin,network_type *network,bushes_type *bushes);
void newtonFlowShift(int j, merge_type *merge, int origin,
                     network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters);

/* Utility functions */
bool isInBush(int origin, int ij, network_type *network, bushes_type *bushes);
bool isMergeNode(int origin, int i, bushes_type *bushes);
int pred2merge(int ij);
int merge2pred(int m);
void exactCostUpdate(int ij, double shift, network_type *network);
void linearCostUpdate(int ij, double shift, network_type *network);
void noCostUpdate(int ij, double shift, network_type *network);
void checkFlows(network_type *network, bushes_type *bushes);

/**** Merges and merge-doubly linked lists ****/

merge_type *createMerge(int numApproaches);
void deleteMerge(merge_type *merge);
void displayMerge(int minVerbosity, merge_type *merge, network_type *network);

typedef struct mergeDLLelt_s {
    int node;
    merge_type *merge;
    struct mergeDLLelt_s *next;
    struct mergeDLLelt_s *prev;
} mergeDLLelt;

typedef struct {
    mergeDLLelt *head;
    mergeDLLelt *tail;
    int size;
} mergeDLL;

mergeDLL *createMergeDLL();
mergeDLLelt *insertMergeDLL(mergeDLL *list, merge_type *merge, int i,
                            mergeDLLelt *after);
void deleteMergeDLL(mergeDLL *list);
void deleteMergeDLLelt(mergeDLL *list, mergeDLLelt *elt);
void displayMergeDLL(int minVerbosity, mergeDLL *list);

/**** Batching ****/
void writeBushes(network_type *network, bushes_type *bushes, char *filename);
void readBushes(network_type *network, bushes_type **bushes, char *filename);
void readBinaryMatrix(network_type *network,
                      algorithmBParameters_type *parameters);
#endif
