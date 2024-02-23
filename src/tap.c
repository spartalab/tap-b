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
this, the network struct needs to be passed aint with the arc pointer.
*/
int arcNumber(network_type *network, arc_type *arc) {
    return arc - network->arcs;
}

/*
BeckmannFunction calculates the Beckmann function for a network given its
current class-flows.  
*/
double BeckmannFunction(network_type *network) {
    double Beckmann = 0;
    int ij, c;
    for (ij = 0; ij < network->numArcs; ij++) {
        Beckmann += network->arcs[ij].calculateInt(&network->arcs[ij], FALSE);
        for (c = 0; c < network->numClasses; c++) {
            Beckmann += network->arcs[ij].classFlow[c]
                       * network->arcs[ij].classCost[c];
        }
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

/* All these functions are commented out until the conic delay function has
 * been implemented. */

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
double generalBPRint(struct arc_type *arc, bool includeFixedCost) {
   if (arc->flow <= 0) return 0; /* Protect against negative flow values and 
                                    0^0 errors */
   return arc->flow * (includeFixedCost == TRUE ? arc->fixedCost : 0
           + arc->freeFlowTime * 
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

double linearBPRint(struct arc_type *arc, bool includeFixedCost) {
   return arc->flow * (includeFixedCost == TRUE ? arc->fixedCost : 0
           + arc->freeFlowTime * (1 + arc->flow*arc->alpha/arc->capacity/2));
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

double quarticBPRint(struct arc_type *arc, bool includeFixedCost) {
   double y = arc->flow / arc->capacity;
   y *= y;
   y *= y;
   return arc->flow * (includeFixedCost == TRUE ? arc->fixedCost : 0
           + arc->freeFlowTime * (1 + arc->alpha * y / 5));
}

double conicCost(struct arc_type *arc) {
    double time = arc->freeFlowTime + arc->fixedCost;
        
    arc->oldRoot = sqrt(arc->b * arc->b + arc->a * arc->a
                        * (1 - arc->flow/arc->capacity + arc->e)
                        * (1 - arc->flow/arc->capacity + arc->e));
    /* Add conic delay */
    time += arc->freeFlowTime *
            (arc->oldRoot -arc->h0
              + arc->a * arc->flow/arc->capacity );
    
    /* Add signalized delay */
    if (arc->flow/arc->saturationFlow <= 0.875) {
        time += arc->sParam / (1 - arc->flow/arc->saturationFlow);
    } else if (arc->flow/arc->saturationFlow < 0.925) {
        time += arc->CD + arc->flow/arc->saturationFlow *
                    (arc->CC + arc->flow/arc->saturationFlow *
                        (arc->CB + arc->flow/arc->saturationFlow * arc->CA));
    } else {
        time += arc->sParam / 0.1;
    }

    /* Add unsignalized delay */
        time += arc->m + arc->u * arc->flow/arc->capacity;
    
    return time;
}

/* Note: To save on computation time, the derivative re-uses the root last
 * calculated in the objective. */
double conicDer(struct arc_type *arc) {
    double der = 0;
    
    /* Start with conic part... use saved root */
    der += arc->a * arc->freeFlowTime / arc->capacity
           * (1 - arc->a*(1 - arc->flow/arc->capacity + arc->e) / arc->oldRoot);

    /* Add signalized delay */
    if (arc->flow/arc->saturationFlow <= 0.875) {
        der += arc->sParam / (arc->saturationFlow
                * (1 - arc->flow/arc->saturationFlow) 
                * (1 - arc->flow/arc->saturationFlow));
    } else if (arc->flow/arc->saturationFlow < 0.925) {
        der += (arc->CC + arc->flow/arc->saturationFlow *
                 (2 * arc->CB + arc->flow/arc->saturationFlow * 3 *arc->CA))
                 / arc->saturationFlow;
    }

    /* Add unsignalized delay */
    der += arc->u /arc->capacity;

    return der;
}

double conicInt(struct arc_type *arc, bool includeFixedCost) {
    fatalError("Conic integrals not yet implemented; set calculateBeckmann "
               "to FALSE in parameters.");
    /* Suppress compiler warning about unused arguments */
    displayMessage(FULL_DEBUG, "%p%d", arc, (int) includeFixedCost);
    return IS_MISSING;
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
    int r, j, c, originNode;
    double sptt = 0;
    declareVector(double, SPcosts, network->numNodes);
    for (r = 0; r < network->numOrigins; r++) {
        originNode = origin2node(network, r);
        c = origin2class(network, r);
        changeFixedCosts(network, c);
        BellmanFord_NoLabel(originNode, SPcosts, network, DEQUE, NULL, NULL);
        for (j = 0; j < network->numZones; j++) {
            sptt += network->demand[r][j] * SPcosts[j];
        }
    }
    deleteVector(SPcosts);
    return sptt;
}

/*
TSTT calculates the total generalized cost on the network, that is, the dot
product of link cost and travel time.
*/
double TSTT(network_type *network) {
    double sum = 0;
    int c;
    for (c = 0; c < network->numClasses; c++) {
        sum += classGeneralizedCost(network, c);
    }
    if (isnan(sum)) displayNetwork(DEBUG, network); /* Oops.  Indicates some 
                    kind of numerical error (likely division by zero or 0^0) */
    return sum;
}

/*
classCost calculates cost, but ONLY for one specified class.  Note that the
given time/toll/distanceFactors do NOT need to be the same as that for the
class.  This gives flexibility in calculating total toll paid, total travel
time, total travel cost, etc.
*/
double classCost(network_type *network, int class, double timeFactor,
                 double tollFactor, double distanceFactor) {
    double sum = 0;
    int ij;
    
    for (ij = 0; ij < network->numArcs; ij++) {
        sum += network->arcs[ij].classFlow[class]
               * (timeFactor *
                    (network->arcs[ij].cost - network->arcs[ij].fixedCost)
                  + distanceFactor * network->arcs[ij].length
                  + tollFactor * network->arcs[ij].classToll[class]);
    }
    return sum;
}

/* 
The following are specific functions for common classCost calls
*/
double classRevenue(network_type *network, int class) {
    return classCost(network, class, 0, 1, 0);
}

double classDistance(network_type *network, int class) {
    return classCost(network, class, 0, 0, 1);
}

double classTravelTime(network_type *network, int class) {
    return classCost(network, class, 1, 0, 0);
}

double classGeneralizedCost(network_type *network, int class) {
    int ij;
    double sum = 0;
    changeFixedCosts(network, class);
    for (ij = 0; ij < network->numArcs; ij++) {
        sum += network->arcs[ij].classFlow[class] * network->arcs[ij].cost;
    }
    return sum;
}


/*
updateAllCosts recalculates and updates all link costs in the network,
according to costFunction
*/
void updateAllCosts(network_type *network) {
    int i;
    for (i = 0; i < network->numArcs; i++)
      network->arcs[i].cost=network->arcs[i].calculateCost(&network->arcs[i]);
}

/*
updateAllCostDers recalculates and updates all link derivatives in the network,
according to derFunction
*/
void updateAllCostDers(network_type *network) {
    int i;
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
    double sptt = SPTT(network);
    double tstt = TSTT(network);
    if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and SPTT "
                                                "are %f %f\n", tstt, sptt);
    if (network->totalODFlow == 0) {
        warning(LOW_NOTIFICATIONS, "No flow or demand on network\n");
        return 0;
    }
    return ((tstt - sptt) / network->totalODFlow);
}

/*
relativeGap1 calculates the ratio between (TSTT - SPTT) and SPTT.
*/
double relativeGap1(network_type *network) {
    double sptt = SPTT(network);
    double tstt = TSTT(network);
    displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest "
                          "path TSTT: %f\n", tstt, sptt);
    if (tstt < sptt) warning(LOW_NOTIFICATIONS, "Negative gap.  TSTT and "
                                              "denom are %f %f\n", tstt, sptt);
    if (sptt == 0 && tstt == 0) {
        warning(LOW_NOTIFICATIONS, "No flow or demand on network\n");
        return  0;
    } else if(sptt == 0) {
        warning(LOW_NOTIFICATIONS, "SPTT is zero\n");
    }
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
    double sptt = SPTT(network);
    double tstt = TSTT(network);

    // Warning: This is a hack and will give incorrect values if relativeGap2
    // is used with a non-BPR function.  TODO: Fix later
    network->beckmann = BeckmannFunction(network);
    network->beckmannLB = min(network->beckmannLB, network->beckmann + sptt -
                                                   tstt);
    displayMessage(DEBUG, "Current relative gap:\nCurrent TSTT: %f\nShortest "
            "path TSTT: %f\n", tstt, sptt);
    if(network->beckmannLB == 0) {
        warning(LOW_NOTIFICATIONS, "BeckmannLB is zero\n");
    }
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
    int c, i, j;
    declareVector(int, order, network->numNodes);
    declareVector(int, backnode, network->numNodes);
    declareVector(int, forwardnode, network->numNodes);

   /* Run forward and reverse searches to see what links must be created */
    search(network->numNodes - 1, order, backnode, network, FIFO, FORWARD);
    search(network->numNodes - 1, order, forwardnode, network, FIFO, REVERSE);

    int newArcs = 0;
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
            /*displayMessage(DEBUG, "Creating (%d,%d)\n", network->numNodes, 
                           j + 1);*/
            newArcVector[i].tail = network->numNodes - 1;
            newArcVector[i].head = j;
            newArcVector[i].alpha = 0;
            newArcVector[i].beta = 1;
            newArcVector[i].flow = 0;
            newArcVector[i].capacity = ARTIFICIAL;
            newArcVector[i].length = ARTIFICIAL;
            newArcVector[i].freeFlowTime = ARTIFICIAL;
            newArcVector[i].calculateCost = &linearBPRcost;
            newArcVector[i].calculateDer = &linearBPRder;           
            newArcVector[i].calculateInt = &linearBPRint;        
            newArcVector[i].cost = newArcVector[i].freeFlowTime; 
            newArcVector[i].classFlow = newVector(network->numClasses, double);
            newArcVector[i].classCost = newVector(network->numClasses, double);
            newArcVector[i].classToll = newVector(network->numClasses, double);
            for (c = 0; c < network->numClasses; c++) {
                newArcVector[i].classFlow[c] = 0;
                newArcVector[i].classCost[c] = ARTIFICIAL; /* Ban trips... */
                newArcVector[i].classToll[c] = 0; /*...but protect revenue. */
            }
            i++;
        }
        if (forwardnode[j] == NO_PATH_EXISTS) {
            /*displayMessage(DEBUG, "Creating (%d,%d)\n", j + 1,
                           network->numNodes);*/
            newArcVector[i].tail = j;
            newArcVector[i].head = network->numNodes - 1;
            newArcVector[i].alpha = 0;
            newArcVector[i].beta = 1;
            newArcVector[i].flow = 0;
            newArcVector[i].capacity = ARTIFICIAL;
            newArcVector[i].length = ARTIFICIAL;
            newArcVector[i].freeFlowTime = ARTIFICIAL;
            newArcVector[i].calculateCost = &linearBPRcost;
            newArcVector[i].calculateDer = &linearBPRder;           
            newArcVector[i].calculateInt = &linearBPRint;         
            newArcVector[i].cost = newArcVector[i].freeFlowTime;           
            newArcVector[i].classFlow = newVector(network->numClasses, double);
            newArcVector[i].classCost = newVector(network->numClasses, double);
            newArcVector[i].classToll = newVector(network->numClasses, double);
            for (c = 0; c < network->numClasses; c++) {
                newArcVector[i].classFlow[c] = 0;
                newArcVector[i].classCost[c] = ARTIFICIAL; /* Ban trips... */
                newArcVector[i].classToll[c] = 0; /*...but protect revenue. */
            }
            i++;
        }
    }
    deleteVector(network->arcs);
#ifdef PARALLELISM
    for(i = 0; i < network->numArcs; i++) {
        pthread_mutex_destroy(&network->arc_muts[i]);
    }
    deleteVector(network->arc_muts);
#endif
    network->numArcs += newArcs;
    network->arcs = newArcVector;
#ifdef PARALLELISM
    network->arc_muts = newVector(network->numArcs, pthread_mutex_t);
#endif
   /* Regenerate forward/reverse star lists */
    for (j = 0; j < network->numNodes; j++) {
        clearArcList(&(network->nodes[j].forwardStar));
        clearArcList(&(network->nodes[j].reverseStar));
    }
    displayMessage(FULL_NOTIFICATIONS, "Finalizing network setup\n");
    finalizeNetwork(network);
    displayMessage(FULL_NOTIFICATIONS, "Finalized network setup\n");

   /* Clean up and return */
    deleteVector(order);
    deleteVector(backnode);
    deleteVector(forwardnode);
}

