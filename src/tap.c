/*
   This file contains routines applicable across many types of algorithms for
   traffic assignment: calculation of Beckmann function, gap functions
   (relative gap, average excess cost), evaluation of BPR functions, total
   system travel time, and so forth.
*/

#include "tap.h"
#if PARALLELISM
    #include "pthread.h"
#endif

//////////////////////////////////
// General network calculations //
//////////////////////////////////

/*
arcNumber converts a pointer to an arc to the index number for that arc; to do
this, the network struct needs to be passed along with the arc pointer.
*/
long arcNumber(network_type *network, arc_type *arc) {
    return arc - network->arcs;
}

/*
BeckmannFunction calculates the Beckmann function for a network given its
current flows.
*/
double BeckmannFunction(network_type *network) {
    double Beckmann = 0;
    long i;
    for (i = 0; i < network->numArcs; i++) {
        Beckmann += network->arcs[i].calculateInt(&network->arcs[i]);
    }
    return Beckmann;
}

/*
calculateGap is a wrapper function for the various gap calculation functions.
Argument gapFunction indicates which function to call.  RELATIVE_GAP_1 is what
I taught you in class.  RELATIVE_GAP_2 is an alternative definition of relative
gap which appears sometimes in the literature.  AEC is average excess cost.
AEC_OB and MEC are other gap functions which I have not yet implemented.
*/
double calculateGap(network_type *network, gap_type gapFunction) {
    updateAllCosts(network);
    switch (gapFunction) {
        case RELATIVE_GAP_1: return relativeGap1(network); break;
        case RELATIVE_GAP_2: return relativeGap2(network); break;
        case AEC: return averageExcessCost(network); break;
        case AEC_OB: fatalError("Gap type AEC_OB not yet implemented."); break;
        case MEC: fatalError("Gap type MEC not yet implemented."); break;
        case PASS: return INFINITY; break;
        default: fatalError("Unknown gap type %d\n", gapFunction);
    }
    return IS_MISSING; /* This should never be reached; provided to eliminate
                          compiler warnings */
}

/*

*****
NOTE: linkCost is a legacy function not used in the current implementation.
*****

linkCost performs a link-specific calculation, specified by the argument
costFunction.  Most commonly used to evaluate a link performance function such
as the BPR function, but passing different costFunctions allows you to
calculate other quantities with the same function.  The following costFunctions
have been implemented thus far:
   BPR -- evaluates a BPR function using the arc'link s alpha and beta values,
          its capacity, and its current flow
   BPR_DER -- evaluates the *derivative* of the BPR function at the current
              flow (useful for Newton's method)
   BPR_INT -- evaluates the *integral* of the BPR function, between zero and
              the current flow (useful for Beckmann function!)
   SO_BPR -- evaluates the marginal cost of the link (flow + externality).  By
             using this instead of the BPR function, the system optimal
             assignment can be solved with the existing code and no other
             changes.
*/
double linkCost(arc_type *arc, cost_type costFunction) {
   if (arc->capacity <= 0)
       return INFINITY; /* Protect against division by zero */
    switch (costFunction) {
    case BPR:
      if (arc->flow <= 0) 
          return arc->freeFlowTime + arc->fixedCost; 
          // Protect against negative flow values and 0^0 errors
        return arc->fixedCost + arc->freeFlowTime * 
            (1 + arc->alpha * pow(arc->flow / arc->capacity, arc->beta));
    case BPR_DER:
      if (arc->flow <= 0) { /* Protect against negative flow values and 0^0 
                               errors */
         if (arc->beta != 1)
            return 0;
         else
            return arc->freeFlowTime * arc->alpha / arc->capacity;
      }
        return arc->freeFlowTime * arc->alpha * arc->beta / arc->capacity
            * pow(arc->flow / arc->capacity, arc->beta - 1);
    case BPR_INT:
      if (arc->flow <= 0) return 0; /* Protect against negative flow values and
                                       0^0 errors */
      return arc->flow * (arc->fixedCost + arc->freeFlowTime * 
              (1 + arc->alpha / (arc->beta + 1) * 
               pow(arc->flow / arc->capacity, arc->beta)));
    default:
        fatalError("Unknown cost function.");
        return IS_MISSING;
    }
}

/*
 * generalBPRcost -- Evaluates the BPR function for an arbitrary polynomial.
 */
double generalBPRcost(struct arc_type *arc) {
   if (arc->flow <= 0) 
   // Protect against negative flow values and 0^0 errors
       return arc->freeFlowTime + arc->fixedCost;

   return arc->fixedCost + arc->freeFlowTime * 
       (1 + arc->alpha * pow(arc->flow / arc->capacity, arc->beta));
}

/*
 * generalBPRder -- Evaluate derivative of an arbitrary polynomial BPR
 * function.
 */
double generalBPRder(struct arc_type *arc) {
   if (arc->flow <= 0) { /* Protect against negative flow values and 0^0 
                            errors */
      if (arc->beta != 1)
         return 0;
      else
         return arc->freeFlowTime * arc->alpha / arc->capacity;
   }
    return arc->freeFlowTime * arc->alpha * arc->beta / arc->capacity
              * pow(arc->flow / arc->capacity, arc->beta - 1);
}

/* 
 * generalBPRint -- Evaluate integral of an arbitrary polynomial BPR function.
 */
double generalBPRint(struct arc_type *arc) {
   if (arc->flow <= 0) return 0; /* Protect against negative flow values and 
                                    0^0 errors */
   return arc->flow * (arc->fixedCost + arc->freeFlowTime * 
           (1 + arc->alpha / (arc->beta + 1) * 
            pow(arc->flow / arc->capacity, arc->beta)));
}

/* linearBPRcost/der/int -- Faster implementation for linear BPR functions. */
double linearBPRcost(struct arc_type *arc) {
   return arc->fixedCost + arc->freeFlowTime * 
       (1 + arc->alpha * arc->flow / arc->capacity);
}

double linearBPRder(struct arc_type *arc) {
   return arc->freeFlowTime * arc->alpha / arc->capacity;
}

double linearBPRint(struct arc_type *arc) {
   return arc->flow * (arc->fixedCost + arc->freeFlowTime * 
           (1 + arc->alpha / arc->capacity / 2));
}

/* quarticBPRcost/der/int -- Faster implementation for 4th-power BPR functions
 */
double quarticBPRcost(struct arc_type *arc) {
   double y = arc->flow / arc->capacity;
   y *= y;
   y *= y;
   return arc->fixedCost + arc->freeFlowTime * (1 + arc->alpha * y);
}

double quarticBPRder(struct arc_type *arc) {
   double y = arc->flow / arc->capacity / arc->capacity;
   y *= y;
   y *= arc->flow;
   return 4 * arc->freeFlowTime * arc->alpha * y;
   
}

double quarticBPRint(struct arc_type *arc) {
   double y = arc->flow / arc->capacity;
   y *= y;
   y *= y;
   return arc->flow * (arc->fixedCost + arc->freeFlowTime * 
           (1 + arc->alpha * y / 5));
}

/*
SPTT calculates the shortest-path travel time on the network, that is, the
total system travel time if everyone could be loaded on the current shortest
paths without changing the travel times.  This function actually re-solves
shortest paths for each origin, rather than using the OD cost field in the
network struct -- this allows SPTT to be calculated without changing any values
in the network, at the expense of a little more run time.
*/
double SPTT(network_type *network) {
    long i, j;
    double sptt = 0;
    declareVector(double, SPcosts, network->numNodes);
    declareVector(long, backnode, network->numNodes);
    declareVector(double, oldCosts, network->numArcs);
    for (i = 0; i < network->numArcs; i++) { 
    // Save old costs (we will be using new ones based on costFunction)
      oldCosts[i] = network->arcs[i].cost;
      network->arcs[i].cost = network->arcs[i].calculateCost(&network->arcs[i]);
   }
    
    for (i = 0; i < network->numZones; i++) {
        BellmanFord(i, SPcosts, backnode, network, DEQUE);
        for (j = 0; j < network->numZones; j++) {
         if (SPcosts[j] > 9e9 && network->OD[i][j].demand > 0)
             displayMessage(LOW_NOTIFICATIONS, "%d -> %d: %f %f\n", i, j, 
                     SPcosts[j], network->OD[i][j].demand);
            sptt += network->OD[i][j].demand * SPcosts[j];
        }
    }
    for (i = 0; i < network->numArcs; i++) // Restore old costs
      network->arcs[i].cost = oldCosts[i];
    deleteVector(SPcosts);
    deleteVector(backnode);
    deleteVector(oldCosts);
    return sptt;
}

#if PARALLELISM

struct thread_args {
    int id;
    int start;
    int num_points;
    network_type *network;
	long* backnode;
	double* SPcosts;
	double* threadSPTT;
};

void BellmanFord_par(struct thread_args* args) {
    int start = args->start;
    int end = args->num_points + start;
    int thread_id = args->id;
    double* SPcosts = args->SPcosts;
    long* backnode = args->backnode;
    double* threadSPTT = args->threadSPTT;
    network_type *network = args->network;

    for (; start < end && start < network->numZones; start++) {
        BellmanFord(start, SPcosts, backnode, network, DEQUE);
        for (int j = 0; j < network->numZones; j++) {
            if (SPcosts[j] > 9e9 && network->OD[start][j].demand > 0)
                displayMessage(LOW_NOTIFICATIONS, "%d -> %d: %f %f\n", start, j, 
                        SPcosts[j], network->OD[start][j].demand);
            threadSPTT[thread_id] += network->OD[start][j].demand * SPcosts[j];
        }
    }
}

/*
SPTT_par calculates the shortest-path travel time on the network, that is, the
total system travel time if everyone could be loaded on the current shortest
paths without changing the travel times.  This function actually re-solves
shortest paths for each origin, rather than using the OD cost field in the
network struct -- this allows SPTT to be calculated without changing any values
in the network, at the expense of a little more run time.
*/
double SPTT_par(network_type *network) {
    long i;
    double sptt = 0;
    declareMatrix(double, SPcosts, 8, network->numNodes);
    declareVector(long, backnode, network->numNodes);
    declareVector(double, oldCosts, network->numArcs);
    declareVector(double, threadSPTT, 8);

    for (i = 0; i < network->numArcs; i++) { 
    // Save old costs (we will be using new ones based on costFunction)
      oldCosts[i] = network->arcs[i].cost;
      network->arcs[i].cost = network->arcs[i].calculateCost(&network->arcs[i]);
   }
    
    pthread_t handles[8];
    struct thread_args args[8];
    int pointsPerThread = network->numZones/8;

    for (int j = 0; j < 8; ++j) {
        args[j].id = j;
        args[j].start = pointsPerThread * j;
        args[j].num_points =  pointsPerThread;
        args[j].network = network;
        args[j].backnode = backnode;
        args[j].SPcosts = SPcosts[j];
        args[j].threadSPTT = threadSPTT;
    }

    for (int j = 0; j < 8; ++j) {
        pthread_create(&handles[j], NULL, (void* (*)(void*)) BellmanFord_par, &args[j]);
    }
    for(int i = 0; i < 8; i++) {
        pthread_join(handles[i], NULL);
    }

    for (int i = 0; i < 8; i++) {
        sptt += threadSPTT[i];
    }

    for (i = 0; i < network->numArcs; i++) // Restore old costs
      network->arcs[i].cost = oldCosts[i];
    deleteMatrix(SPcosts);
    deleteVector(backnode);
    deleteVector(oldCosts);
    deleteVector(threadSPTT)
    return sptt;
}
#endif

/*
TSTT calculates the total system travel time on the network, that is, the dot
product of link flow and travel time.
*/
double TSTT(network_type *network) {
    double sum = 0;
    long i;
    for (i = 0; i < network->numArcs; i++)
        sum += network->arcs[i].flow * 
            network->arcs[i].calculateCost(&network->arcs[i]);
   if (isnan(sum)) displayNetwork(DEBUG, network); /* Oops.  Indicates some 
                   kind of numerical error (likely division by zero or 0^0) */
    return sum;
}

/*
updateAllCosts recalculates and updates all link costs in the network,
according to costFunction
*/
void updateAllCosts(network_type *network) {
    long i;
    for (i = 0; i < network->numArcs; i++)
      network->arcs[i].cost=network->arcs[i].calculateCost(&network->arcs[i]);
}

/*
updateAllCostDers recalculates and updates all link derivatives in the network,
according to derFunction
*/
void updateAllCostDers(network_type *network) {
    long i;
    for (i = 0; i < network->numArcs; i++)
      network->arcs[i].der = network->arcs[i].calculateDer(&network->arcs[i]);
}

///////////////////
// GAP FUNCTIONS //
///////////////////

/*
averageExcessCost calculates the difference between TSTT and SPTT, normalized
by total demand in the network.
*/
double averageExcessCost(network_type *network) {
    double sptt = SPTT_par(network), tstt = TSTT(network);
    if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and SPTT "
                                                "are %f %f\n", tstt, sptt);
    return ((tstt - sptt) / network->totalODFlow);
}

/*
relativeGap1 calculates the ratio between (TSTT - SPTT) and SPTT.
*/
double relativeGap1(network_type *network) {
    double sptt = SPTT_par(network), tstt = TSTT(network);
    displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest "
                          "path TSTT: %f\n", tstt, sptt);
    if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and "
                                              "denom are %f %f\n", tstt, sptt);
    return (tstt / sptt - 1);
}

/*
relativeGap2 calculates the ratio between the current value of the Beckmann
function and a lower bound on the optimal value of the Beckmann function, which
is calculated using convexity: the tangent plane at any point always
underestimates the true value of the Beckmann function, and one such value is
given by the current value of the Beckmann function, minus (TSTT - SPTT).
(This is because the gradient of the Beckmann function is just the current link
travel times; evaluating the tangent plane approximation at the target link
flow solution x* gives this formula.)  Since any lower bound will do, this
function stores the best lower bound seen thus far.
*/
double relativeGap2(network_type *network) {
    double sptt = SPTT_par(network), tstt = TSTT(network);

    // Warning: This is a hack and will give incorrect values if relativeGap2
    // is used with a non-BPR function.  TODO: Fix later
    network->beckmann = BeckmannFunction(network); 
    network->beckmannLB = min(network->beckmannLB, network->beckmann + sptt - 
                                                   tstt);
    displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest "
            "path TSTT: %f\n", tstt, sptt);
    return (network->beckmann / network->beckmannLB - 1);
}


/*
makeStronglyConnectedNetwork creates artificial links, if needed, to ensure
that every link in the network is reachable from every other.  It does this by
calling 'search' in both the forward and reverse directions from an arbitrarily
chosen node (the one with the highest index), call this node *.  Artificial
links are created to and from * and any nodes which cannot be found from these
searches.  As a result, the network becomes strongly connected: any two nodes
are at least reachable through the path i -> * -> j.  The freeFlowTime (and
various other parameters) on these links are set equal to the symbolic constant
ARTIFICIAL, which should be set high enough that the links will not be used.
*/
void makeStronglyConnectedNetwork(network_type *network) {
    long i, j;
    declareVector(long, order, network->numNodes);
    declareVector(long, backnode, network->numNodes);
    declareVector(long, forwardnode, network->numNodes);

   /* Run forward and reverse searches to see what links must be created */
    search(network->numNodes - 1, order, backnode, network, FIFO, FORWARD);
    search(network->numNodes - 1, order, forwardnode, network, FIFO, REVERSE);

    long newArcs = 0;
    for (i = 0; i < network->numNodes; i++) {
        if (backnode[i] == NO_PATH_EXISTS) newArcs++;
        if (forwardnode[i] == NO_PATH_EXISTS) newArcs++;
    }
    if (newArcs == 0) { /* Nothing to do, network is already strongly 
                           connected, so clean up/return */
        deleteVector(order);
        deleteVector(backnode);
        deleteVector(forwardnode);
        return;
    }

   /* Create new arc vector, with all the old ones plus the new artificial 
    * links */
    displayMessage(FULL_NOTIFICATIONS, "Warning: Creating %d artifical arcs "
            "to ensure strong connectivity.\n", newArcs);
    declareVector(arc_type, newArcVector, network->numArcs + newArcs);
    for (i = 0; i < network->numArcs; i++) {
        newArcVector[i] = network->arcs[i];
    }
    for (j = 0; j < network->numNodes; j++) {
        if (backnode[j] == NO_PATH_EXISTS) {
            displayMessage(DEBUG, "Creating (%d,%d)\n", network->numNodes, 
                           j + 1);
            newArcVector[i].tail = network->numNodes - 1;
            newArcVector[i].head = j;
            newArcVector[i].alpha = 0;
            newArcVector[i].beta = 1;
            newArcVector[i].flow = 0;
            newArcVector[i].capacity = ARTIFICIAL;
            newArcVector[i].length = ARTIFICIAL;
            newArcVector[i].toll = ARTIFICIAL;
            newArcVector[i].freeFlowTime = ARTIFICIAL;
           newArcVector[i].calculateCost = &linearBPRcost;
           newArcVector[i].calculateDer = &linearBPRder;           
           newArcVector[i].calculateInt = &linearBPRint;        
           newArcVector[i].cost = newArcVector[i].freeFlowTime; 
            i++;
        }
        if (forwardnode[j] == NO_PATH_EXISTS) {
            displayMessage(DEBUG, "Creating (%d,%d)\n", j + 1,
                           network->numNodes);
            newArcVector[i].tail = j;
            newArcVector[i].head = network->numNodes - 1;
            newArcVector[i].alpha = 0;
            newArcVector[i].beta = 1;
            newArcVector[i].flow = 0;
            newArcVector[i].capacity = ARTIFICIAL;
            newArcVector[i].length = ARTIFICIAL;
            newArcVector[i].toll = ARTIFICIAL;
            newArcVector[i].freeFlowTime = ARTIFICIAL;
           newArcVector[i].calculateCost = &linearBPRcost;
           newArcVector[i].calculateDer = &linearBPRder;           
           newArcVector[i].calculateInt = &linearBPRint;         
           newArcVector[i].cost = newArcVector[i].freeFlowTime;           
            i++;
        }
    }
    network->numArcs += newArcs;
    deleteVector(network->arcs);
    network->arcs = newArcVector;

   /* Regenerate forward/reverse star lists */
    for (j = 0; j < network->numNodes; j++) {
        clearArcList(&(network->nodes[j].forwardStar));
        clearArcList(&(network->nodes[j].reverseStar));
    }
    finalizeNetwork(network);

   /* Clean up and return */
    deleteVector(order);
    deleteVector(backnode);
    deleteVector(forwardnode);
}

