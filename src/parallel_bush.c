#include "parallel_bush.h"

/*
 * scanBushes_par -- simultaneously find shortest and longest paths for a given
 * bush by a thread.  If the longestUsed argument is TRUE, the longest path search will
 * restrict attention to longest *used* paths.  If FALSE, longest bush paths
 * will be calculated regardless of whether there is flow on them.
*/
void scanBushes_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters, scan_type LPrule) {
//    displayMessage(FULL_NOTIFICATIONS, "Top of scan %d\n", t_id);
    int h, i, hi, m, curnode, curarc;
    double tempcost;
    merge_type *merge;

    for (i = 0; i < network->numNodes; i++) {
        bushes->LPcost_par[origin][i] = -INFINITY;
        bushes->SPcost_par[origin][i] = INFINITY;
    }

    /* Ensure costs are up to date */
    if (parameters->linkShiftB != &exactCostUpdate) updateAllCosts(network);

    bushes->SPcost_par[origin][origin2node(network, origin)] = 0;
    bushes->LPcost_par[origin][origin2node(network, origin)] = 0;
    for (curnode = 1; curnode < network->numNodes; curnode++) {
        i = bushes->bushOrder[origin][curnode];
        /* Iterate over incoming links */
        if (isMergeNode(origin, i, bushes) == TRUE) {
            m = pred2merge(bushes->pred[origin][i]);
            merge = bushes->merges[origin][m];
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->SPcost_par[origin][h] + network->arcs[hi].cost;
                if (tempcost < bushes->SPcost_par[origin][i]) {
                    bushes->SPcost_par[origin][i] = tempcost;
                    merge->SPlink = curarc;
                }
            }

            /* Find longest (need to separate out depending on LPrule) */
            if (LPrule == NO_LONGEST_PATH) continue;
            for (curarc = 0; curarc < merge->numApproaches; curarc++) {
                hi = merge->approach[curarc];
                h = network->arcs[hi].tail;
                tempcost = bushes->LPcost_par[origin][h] + network->arcs[hi].cost;     
                if (tempcost > bushes->LPcost_par[origin][i]
                  && (LPrule == LONGEST_BUSH_PATH
                    || (LPrule == LONGEST_USED_PATH
                      && merge->approachFlow[curarc] > parameters->minLinkFlow)
                    || (LPrule == LONGEST_USED_OR_SP
                      &&(merge->approachFlow[curarc] > parameters->minLinkFlow
                      || merge->SPlink == curarc)))) {
                    
                    bushes->LPcost_par[origin][i] = tempcost;
                    merge->LPlink = curarc;
                }
            }
        } else { /* Only one incoming bush link, not much to do */
            hi = bushes->pred[origin][i];
            h = network->arcs[hi].tail;
            bushes->LPcost_par[origin][i] = bushes->LPcost_par[origin][h] + network->arcs[hi].cost;
            bushes->SPcost_par[origin][i] = bushes->SPcost_par[origin][h] + network->arcs[hi].cost;
        }
    }
}

/*
 * updateBushB_par -- updates the topology of a bush: deleting unused links, and
 * adding links which can serve as shortcuts.  Follows Nie (2009) for how this
 * is done, first checking a stronger criterion for adding bush links, then as
 * a fallback checking a weaker one.
 *
 * New links are marked using the NEW_LINK constant, for use in
 * reconstructMerges (the function which rebuilds the merge data structures
 * after bush updating.)
 */
void updateBushB_par(int origin, network_type *network, bushes_type *bushes,
                 algorithmBParameters_type *parameters) {

    int ij, i, j, newArcs = 0;

    /* First update labels... ignoring longest unused paths since those will be
     * removed in the next step. */
    scanBushes_par(origin, network, bushes, parameters, LONGEST_USED_OR_SP);
    calculateBushFlows_par(origin, network, bushes);

    /* Make a first pass... */
    for (ij = 0; ij < network->numArcs; ij++) {
        /* Mark links with near-zero contribution for removal by setting to
         * exactly zero */
        if (bushes->flow_par[origin][ij] < parameters->minLinkFlow) {
            if (isInBush(origin, ij, network, bushes) == TRUE)
                displayMessage(FULL_DEBUG, "Attempting to delete (%d,%d)\n", 
                           network->arcs[ij].tail+1, network->arcs[ij].head+1);
            bushes->flow_par[origin][ij] = 0;
        }
        if (bushes->flow_par[origin][ij] > 0) continue; /* Link is already in the bush, no
                                               need to add it again */
        /* See if the link provides a shortcut using the strict criterion */
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        if (bushes->LPcost_par[origin][i] == -INFINITY && bushes->LPcost_par[origin][j] > -INFINITY)
            continue; /* No path to extend */
        if (bushes->SPcost_par[origin][i] + network->arcs[ij].cost < bushes->SPcost_par[origin][j]
            && bushes->LPcost_par[origin][i] < bushes->LPcost_par[origin][j]
            && (network->arcs[ij].tail == origin2node(network, origin)
                || network->arcs[ij].tail >= network->firstThroughNode))
        {
            bushes->flow_par[origin][ij] = NEW_LINK;
            newArcs++;
            /* Never delete shortest path tree... should be OK with floating point
             * comparison since this is how SPcost is calculated */
        } else if (bushes->SPcost_par[origin][i]+network->arcs[ij].cost==bushes->SPcost_par[origin][j]
                   && bushes->flow_par[origin][ij] == 0
                   && isInBush(origin, ij, network, bushes) == TRUE) {
            bushes->flow_par[origin][ij] = NEW_LINK;
        }
    }
//
//    /* If strict criterion fails, try a looser one */
//    if (newArcs == 0) {
//        for (ij = 0; ij < network->numArcs; ij++) {
//            i = network->arcs[ij].tail;
//            j = network->arcs[ij].head;
//            if (bushes->LPcost_par[origin][i]==-INFINITY && bushes->LPcost_par[origin][j]>-INFINITY)
//                continue; /* No path to extend */
//            if (bushes->flow_par[origin][ij] == 0 && bushes->LPcost_par[origin][i] < bushes->LPcost_par[origin][j]
//                && (network->arcs[ij].tail == origin2node(network, origin)
//                    || network->arcs[ij].tail >= network->firstThroughNode))
//            {
//                bushes->flow_par[origin][ij] = NEW_LINK;
//            }
//        }
//    }

    /* Finally update bush data structures: delete/add merges, find a new
     * topological order, rectify approach proportions */
    reconstructMerges_par(origin, network, bushes);
    parameters->topologicalOrder(origin, network, bushes, parameters);
}

/*
 * reconstructMerges_par -- update the merge data structures after a bush is
 * updated by adding and removing links.  In particular, for each node we do
 * the following...
 *
 * 1. Normalize approach proportions for incoming links -- if there are zero
 *    such links for a non-origin, throw an error.
 * 2. Set pred if there is just one incoming link, otherwise create a merge.
 *
 * Stores merges in a linked list at first, so they can be created in one pass.
 * Then transfer into an array for fast indexing.
 */
void reconstructMerges_par(int origin, network_type *network, bushes_type *bushes){
//    displayMessage(FULL_NOTIFICATIONS, "Top of update reconstructmerg %d\n", t_id);

    int i, hi, lastApproach, m, arc, numApproaches;
    arcListElt *curArc;
    merge_type *merge;
    mergeDLL *mergeList = createMergeDLL();
    mergeDLLelt *curMerge;


    /* Create necessary merges */
    for (i = 0; i < network->numNodes; i++) {
        if (i == origin2node(network, origin)) continue;
        numApproaches = 0;
        for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
             curArc = curArc->next) {
            hi = ptr2arc(network, curArc->arc);
            if (bushes->flow_par[origin][hi] > 0 || bushes->flow_par[origin][hi] == NEW_LINK) {
                numApproaches++;
                lastApproach = hi;
            }
        }
        if (numApproaches == 0)
            fatalError("Cannot have non-origin node %d in bush %d without"
                       "incoming contributing links", i, origin);
        if (numApproaches == 1) { /* No merge */
            bushes->pred[origin][i] = lastApproach;
        } else { /* Must create a merge */
            merge = createMerge(numApproaches);
            arc = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                 curArc = curArc->next) {
                hi = ptr2arc(network, curArc->arc);
                if (bushes->flow_par[origin][hi] > 0 || bushes->flow_par[origin][hi] == NEW_LINK) {
                    if (bushes->flow_par[origin][hi] == NEW_LINK) bushes->flow_par[origin][hi] = 0;
                    merge->approach[arc] = hi;
                    merge->approachFlow[arc] = bushes->flow_par[origin][hi];
                    arc++;
                }
            }
            insertMergeDLL(mergeList, merge, i, mergeList->tail);
        }
    }

    /* Now transfer to array, deleting old merges and replacing with new ones*/
    for (m = 0; m < bushes->numMerges[origin]; m++) {
        deleteMerge(bushes->merges[origin][m]);
    }
    deleteVector(bushes->merges[origin]);
    bushes->numMerges[origin] = mergeList->size;
    bushes->merges[origin] = newVector(mergeList->size, merge_type *);
    m = 0;
    for (curMerge = mergeList->head; curMerge != NULL;
         curMerge = curMerge->next) {
        bushes->merges[origin][m] = curMerge->merge;
        bushes->pred[origin][curMerge->node] = merge2pred(m);
        m++;
    }

    deleteMergeDLL(mergeList);
}

/*
 * updateFlowsB -- shift flows on a bush to move closer to equilibrium.
 * This function is fairly short, as it mainly just manages a few other
 * functions that do the actual work.
 * - calculateBushFlows gives current flows (remember, these are not stored
 *   persistently).
 * - rescanAndCheck both updates longest/shortest path labels, but also checks
 *   whether these are different enough to actually do anything -- depending
 *   on parameters we can skip bushes which are closer to equilibrium, to
 *   spend more time on bushes further away.
 * - updateFlowPass then does the actual flow shifting, inside of a loop in
 *   case we want to do multiple shifts per origin.
 */
bool updateFlowsB_par(int origin, network_type *network, bushes_type *bushes,
                  algorithmBParameters_type *parameters) {
    int i;

    /* Recompute bush flows for this origin */
    calculateBushFlows_par(origin, network, bushes);


    /* Update longest/shortest paths, check whether there is work to do */
    if (rescanAndCheck_par(origin, network, bushes, parameters) == FALSE) {
        bushes->updateBush[origin] = FALSE;
        displayMessage(DEBUG, "bailing out\n");
        return FALSE;
    };

    /* Now do (possibly multiple) shifts per origin */
    for (i = 0; i < parameters->shiftReps; i++) {
        updateFlowPass_par(origin, network, bushes, parameters);

        // Too lazy to make this atomic rn
        pthread_mutex_lock(&shift_lock);
        parameters->numFlowShifts++;
        pthread_mutex_unlock(&shift_lock);

        /* Uncomment next line for extra validation checking */
//        checkFlows_par(network, bushes, t_id);
        if (parameters->rescanAfterShift == TRUE
            && i + 1 < parameters->shiftReps) {
            if (rescanAndCheck_par(origin, network, bushes, parameters) == FALSE)
                return TRUE;
        }
    }

    return TRUE;

}

/*
 * rescanAndCheck -- calls functions updating the longest/shortest path labels
 * on the bush, along with divergence nodes; and then checks whether the bush
 * is close enough to equilibrium to skip for now (in which case the function
 * returns FALSE).  If we need to shift flows, it returns TRUE.
 */
bool rescanAndCheck_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    int i, j, ij;
    double maxgap = 0;
    double bushSPTT = 0, bushExcess = 0;

    scanBushes_par(origin, network, bushes, parameters, LONGEST_USED_PATH);
    findDivergenceNodes(origin, network, bushes);
    for (i = 0; i < network->numZones; i++) {
        bushSPTT += network->demand[origin][i] * bushes->SPcost_par[origin][i];
        maxgap = max(maxgap, fabs(bushes->LPcost_par[origin][i] -bushes->SPcost_par[origin][i]));
    }
    for (; i < network->numNodes; i++) {
        maxgap = max(maxgap, fabs(bushes->LPcost_par[origin][i] -bushes->SPcost_par[origin][i]));
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        i = network->arcs[ij].tail;
        j = network->arcs[ij].head;
        bushExcess += bushes->flow_par[origin][ij] * (network->arcs[ij].cost +
                                          bushes->SPcost_par[origin][i] - bushes->SPcost_par[origin][j]);
    }
    displayMessage(DEBUG, "Scanning %d, gap is %f\n", origin,
                   bushExcess / bushSPTT);
    displayMessage(DEBUG, "Max gap, threshold, threshold AEC: %f %f %f\n",
                   maxgap, parameters->thresholdGap, parameters->thresholdAEC);
    if (maxgap < parameters->thresholdGap) return FALSE;
    if (bushExcess / bushSPTT < parameters->thresholdAEC) return FALSE;

    return TRUE;
}

/*
 * updateFlowPass_par -- does the actual work of shifting flows from longest to
 * shortest paths, through a pass in descending topological order.
 */
void updateFlowPass_par(int origin, network_type *network, bushes_type *bushes,
                    algorithmBParameters_type *parameters) {
    int i, j, m, node;
    merge_type *merge;


    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++)
        bushes->nodeFlow_par[origin][i] = network->demand[origin][i];

    for (; i < network->numNodes; i++) bushes->nodeFlow_par[origin][i] = 0;


    /* Descending pass for flow shifts */
    for (node = network->numNodes - 1; node > 0; node--) {
        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple_par(j, origin, network, bushes);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            if (merge->LPlink == merge->SPlink
                || (fabs(bushes->LPcost_par[origin][j] - bushes->SPcost_par[origin][j])
                        < parameters->minCostDifference)
                || bushes->nodeFlow_par[origin][j] < parameters->minLinkFlowShift
                || merge->LPlink == NO_PATH_EXISTS)
            { /* Also nothing to do */
                pushBackFlowMerge_par(merge, network, bushes, origin);
            } else {
                newtonFlowShift_par(j,merge,origin,network,bushes,parameters);
            }
        }
    }
}

/*
 * calculateBushFlows_par -- a key feature of this implementation is that bush
 * flows are not persistently stored, but generated as needed.  This saves
 * memory at the expense of a little extra computation.  calculateBushFlows
 * does this by making a pass in descending topological order, loading flows
 * from the OD matrix, and splitting flows as needed.
 */
void calculateBushFlows_par(int origin,network_type *network,bushes_type *bushes) {

    int i, j, ij, m, node;
    merge_type *merge;


    /* Initialize node flows with OD matrix */
    for (i = 0; i < network->numZones; i++) {
        bushes->nodeFlow_par[origin][i] = network->demand[origin][i];
    }
    for (; i < network->numNodes; i++) {
        bushes->nodeFlow_par[origin][i] = 0;
    }
    for (ij = 0; ij < network->numArcs; ij++) {
        bushes->flow_par[origin][ij] = 0;
    }

    /* Descending pass for flow calculations  */
    for (node = network->numNodes - 1; node > 0; node--) {

        j = bushes->bushOrder[origin][node];
        if (isMergeNode(origin, j, bushes) == FALSE) { /* Nothing to do */
            pushBackFlowSimple_par(j, origin, network, bushes);
        } else {
            m = pred2merge(bushes->pred[origin][j]);
            merge = bushes->merges[origin][m];
            pushBackFlowMerge_par(merge, network, bushes, origin);
        }
    }

}

/*
 * pushBackFlowSimple_par -- the easy way to "split" flow.  For a non-merge node,
 * just push flow onto the predecessor.
 */
inline void pushBackFlowSimple_par(int j, int origin, network_type *network,
                        bushes_type *bushes) {

    int i, ij;

    ij = bushes->pred[origin][j];
    i = network->arcs[ij].tail;
    bushes->flow_par[origin][ij] = bushes->nodeFlow_par[origin][j];
    bushes->nodeFlow_par[origin][i] += bushes->nodeFlow_par[origin][j];
}

/*
 * pushBackFlowMerge_par -- the harder way to split flow, when there are multiple
 * approaches to a merge node.
 */
inline void pushBackFlowMerge_par(merge_type *merge, network_type *network,
                       bushes_type *bushes, int t_id) {

//    displayMessage(FULL_NOTIFICATIONS, "top pushback merge\n");

    int i, ij, arc;
    double flow;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        ij = merge->approach[arc];
        i = network->arcs[ij].tail;
        flow = merge->approachFlow[arc];
        bushes->flow_par[t_id][ij] = flow;
        bushes->nodeFlow_par[t_id][i] += flow;
    }
}

/*
 * rectifyMerge_par -- to guard against numerical errors at a particular merge
 * node, recalculates the approach proportions.  In the degenerate case of zero
 * node flow, push everything onto the shortest path.
 */
void rectifyMerge_par(int j, merge_type *merge, bushes_type *bushes, int t_id) {

    int arc;
    double totalFlow = 0;

    for (arc = 0; arc < merge->numApproaches; arc++) {
        totalFlow += merge->approachFlow[arc];
    }

    if (totalFlow > 0) {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] *= bushes->nodeFlow_par[t_id][j] / totalFlow;
        }
    } else {
        for (arc = 0; arc < merge->numApproaches; arc++) {
            merge->approachFlow[arc] = 0;
        }
        merge->approachFlow[merge->SPlink] = bushes->nodeFlow_par[t_id][j];
    }

}

/*
 * newtonFlowShift -- For a given merge node, use Newton's method to shift flow
 * from the longer segment to the shorter one.
 */
void newtonFlowShift_par(int j, merge_type *merge, int origin,
                     network_type *network, bushes_type *bushes,
                     algorithmBParameters_type *parameters) {
    double flow1, flow2, cost1, cost2, der1, der2, shift;
    int i, k, hi, m, c;
    merge_type *segmentMerge;


    c = origin2class(network, origin);
    /* Calculate information for longest path segment */
    i = j;
    cost1 = 0;
    der1 = 0;
    flow1 = merge->approachFlow[merge->LPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            flow1=min(flow1,segmentMerge->approachFlow[segmentMerge->LPlink]);
        }
        cost1 += network->arcs[hi].cost;
        der1 += network->arcs[hi].der;
        i = network->arcs[hi].tail;
    }

    /* Do same for shortest path segment */
    i = j;
    cost2 = 0;
    der2 = 0;
    flow2 = merge->approachFlow[merge->SPlink];
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            flow2=min(flow2,segmentMerge->approachFlow[segmentMerge->SPlink]);
        }
        cost2 += network->arcs[hi].cost;
        der2 += network->arcs[hi].der;
        i = network->arcs[hi].tail;
    }

    /* Determine Newton shift, truncating if need be to ensure feasibility */
    if (der1 + der2 == 0) der1 = parameters->minDerivative;
    shift = parameters->newtonStep * (cost1 - cost2) / (der1 + der2);
    shift = min(shift, flow1);
    shift = max(shift, -flow2);

    /* Now apply shift */
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->LPlink];
            segmentMerge->approachFlow[segmentMerge->LPlink] -= shift;
        }
        bushes->flow_par[origin][hi] -= shift;
        if (!exactCostUpdate_par(hi, -shift, network, c)) {
            bushes->flow_par[origin][hi] += shift;
                if(isMergeNode(origin, i, bushes)) {
                    segmentMerge->approachFlow[segmentMerge->LPlink] += shift;
                }
                k = j;
                while (k != i) {
//                    displayMessage(FULL_NOTIFICATIONS, "Exact update failed %d: %f\n",k, shift );
                    if (isMergeNode(origin, k, bushes) == FALSE) {
                        hi = bushes->pred[origin][k];
                    } else {
                        m = pred2merge(bushes->pred[origin][k]);
                        segmentMerge = bushes->merges[origin][m];
                        hi = segmentMerge->approach[segmentMerge->LPlink];
                        segmentMerge->approachFlow[segmentMerge->LPlink] += shift;
                    }
                    exactCostUpdate_par(hi, shift, network, c);
                    bushes->flow_par[origin][hi] += shift;
                    k = network->arcs[hi].tail;
                }
                shift /= 2;
                continue;
            }
        i = network->arcs[hi].tail;
    }
    i = j;
    while (i != merge->divergenceNode) {
        if (isMergeNode(origin, i, bushes) == FALSE) {
            hi = bushes->pred[origin][i];
        } else {
            m = pred2merge(bushes->pred[origin][i]);
            segmentMerge = bushes->merges[origin][m];
            hi = segmentMerge->approach[segmentMerge->SPlink];
            segmentMerge->approachFlow[segmentMerge->SPlink] += shift;
        }
        bushes->flow_par[origin][hi] += shift;
        exactCostUpdate_par(hi, shift, network,c);
        i = network->arcs[hi].tail;
    }
}

/*
 * checkFlows_par -- an error-checking routine which can be invoked to see if flow
 * conservation is properly maintained on the bushes.  In the release
 * implementation, all calls to this function are commented out.
 */
#define FLOW_TOLERANCE 1e-5
void checkFlows_par(network_type *network, bushes_type *bushes, int t_id) {
    int r, i, ij, kl;
    arcListElt *curArc;
    double balance;

    declareVector(double, flowCheck, network->numArcs);
    for (ij = 0; ij < network->numArcs; ij++) {
        flowCheck[ij] = 0;
    }

    for (r = 0; r < network->numZones; r++) {
        calculateBushFlows_par(r, network, bushes);
        /* First check bush flow consistency */
        for (ij = 0; ij < network->numArcs; ij++) {
            flowCheck[ij] += bushes->flow_par[t_id][ij];
            if (bushes->flow_par[t_id][ij] < 0) {
                for (kl = 0; kl < network->numArcs; kl++) {
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                   network->arcs[kl].tail + 1,
                                   network->arcs[kl].head + 1,
                                   bushes->flow[kl]);
                }
                fatalError("Flow validation failed: origin %d, link (%ld,%ld)"
                           "has negative flow.", r, network->arcs[ij].tail + 1,
                           network->arcs[ij].head + 1);
            }
        }
        for (i = 0; i < network->numNodes; i++) {
            if (i == r) continue;
            balance = 0;
            for (curArc = network->nodes[i].reverseStar.head; curArc != NULL;
                 curArc = curArc->next) {
                balance += bushes->flow_par[t_id][ptr2arc(network, curArc->arc)];
            }
            for (curArc = network->nodes[i].forwardStar.head; curArc != NULL;
                 curArc = curArc->next) {
                balance -= bushes->flow_par[t_id][ptr2arc(network, curArc->arc)];
            }
            if (i < network->numZones) balance -= network->demand[r][i];
            if (fabs(balance) > FLOW_TOLERANCE) {
                for (kl = 0; kl < network->numArcs; kl++) {
                    displayMessage(DEBUG, "(%ld,%ld) %f\n",
                                   network->arcs[kl].tail + 1,
                                   network->arcs[kl].head + 1,
                                   bushes->flow_par[t_id][kl]);
                }
                fatalError("Flow validation failed: origin %d, node %ld"
                           "violates conservation.", r, i + 1);

            }
        }

    }

    /* Then check overall link flow consistency */
    for (ij = 0; ij < network->numArcs; ij++) {
        if (fabs(flowCheck[ij] - network->arcs[ij].flow) > FLOW_TOLERANCE) {
            fatalError("Flow validation failed: link (%ld,%ld) has aggregate"
                       "flow %f, but %f when disaggregated.",
                       network->arcs[ij].tail + 1,
                       network->arcs[ij].head + 1,
                       network->arcs[ij].flow, flowCheck[ij]);
        }
    }

    deleteVector(flowCheck);
}

bool exactCostUpdate_par(int ij, double shift, network_type *network, int class) {
    pthread_mutex_lock(&network->arc_muts[ij]);
    if (network->arcs[ij].flow + shift < -FLOW_TOLERANCE) {
        pthread_mutex_unlock(&network->arc_muts[ij]);
        return FALSE;
    }
    network->arcs[ij].classFlow[class] += shift;
    network->arcs[ij].flow += shift;
    network->arcs[ij].cost=network->arcs[ij].calculateCost(&network->arcs[ij]);
    network->arcs[ij].der = network->arcs[ij].calculateDer(&network->arcs[ij]);
    pthread_mutex_unlock(&network->arc_muts[ij]);
    return TRUE;
}

void classUpdate_par(int hi, int class, double shift,  network_type *network) {
    pthread_mutex_lock(&network->arc_muts[hi]);
    network->arcs[hi].classFlow[class] -= shift;
    pthread_mutex_lock(&network->arc_muts[hi]);
}

