#include "datastructures.h"

/*
This file contains implementation for commonly-used data structures, including
singly and doubly linked lists, binary heaps, queues, as well as memory
allocation and deallocation.

You probably don't need to dig into this unless you are particularly curious.
*/

/******************
 ** Linked lists **
 ******************/

/**** Singly linked lists ****/

linkedList *createLinkedList() {
    declareScalar(linkedList, newll);
    newll->head = NULL;
    newll->tail = NULL;
    newll->size = 0;
    return newll;
}

linkedListElt *insertLinkedList(linkedList *list, int value,
                                linkedListElt *after) {
    declareScalar(linkedListElt, newNode);
    newNode->value = value;
    if (after != NULL) { /* Not inserting at head */
        newNode->next = after->next;
        if (list->tail == after) list->tail = newNode;
        after->next = newNode;
    } else { /* Inserting at head */
        newNode->next = list->head;
        if (list->tail == after) list->tail = newNode;
        list->head = newNode;
    }
    list->size++;
    return newNode;
}

void deleteLinkedList(linkedList *list) {
   clearLinkedList(list);
    killScalar(list);
}

void clearLinkedList(linkedList *list) {
    linkedListElt *savenode, *curnode = list->head;
    while (curnode != NULL) {
        savenode = curnode->next;
        killScalar(curnode);
        curnode = savenode;
    }
    list->head = NULL;
    list->tail = NULL;
    list->size = 0;
}

void displayLinkedList(int minVerbosity, linkedList *list) {
    linkedListElt *curnode = list->head;
    displayMessage(minVerbosity,"Start of the list: %p\n", (void *)list->head);
    while (curnode != NULL) {
        displayMessage(minVerbosity, "%p: %d -> %p\n", (void *)curnode,
                       curnode->value, (void *)curnode->next);
        curnode = curnode->next;
    }
    displayMessage(minVerbosity, "End of the list: %p\n", (void *)list->tail);
}

/**** Doubly linked lists ****/

doublyLinkedList *createDoublyLinkedList() {
    declareScalar(doublyLinkedList, newdll);
    newdll->head = NULL;
    newdll->tail = NULL;
    newdll->size = 0;
    return newdll;
}

doublyLinkedListElt *insertDoublyLinkedList(doublyLinkedList *list,
                                            double value,
                                            doublyLinkedListElt *after) {
    declareScalar(doublyLinkedListElt, newNode);
    newNode->value = value;
    if (after != NULL) {
        newNode->prev = after;
        newNode->next = after->next;
        if (list->tail != after)
            newNode->next->prev = newNode;
        else
            list->tail = newNode;
        after->next = newNode;
    } else {
        newNode->prev = NULL;
        newNode->next = list->head;
        if (list->tail != after)
            newNode->next->prev = newNode;
        else
            list->tail = newNode;
        list->head = newNode;
    }
    list->size++;
    return newNode;
}

void deleteDoublyLinkedList(doublyLinkedList *list) {
    while (list->head != NULL)
        deleteDoublyLinkedListElt(list, list->tail);
    killScalar(list);
}

void deleteDoublyLinkedListElt(doublyLinkedList *list,
                               doublyLinkedListElt *elt) {
    if (list->tail != elt) {
        if (list->head != elt)
            elt->prev->next = elt->next;
        else
            list->head = elt->next;
        elt->next->prev = elt->prev;
    } else {
        list->tail = elt->prev;
        if (list->head != elt)
            elt->prev->next = elt->next;
        else
            list->head = elt->next;
    }
    list->size--;
    killScalar(elt);
}

void displayDoublyLinkedList(int minVerbosity, doublyLinkedList *list) {
    doublyLinkedListElt *curnode = list->head;
    displayMessage(minVerbosity,"Start of the list: %p\n", (void *)list->head);
    while (curnode != NULL) {
        displayMessage(minVerbosity, "%p %f %p %p\n", (void *)curnode,
                curnode->value, (void *)curnode->prev, (void *)curnode->next);
        curnode = (*curnode).next;
    }
    displayMessage(minVerbosity, "End of the list: %p\n", (void *)list->tail);
}

/************
 ** Queues **
 ************/

/**** Standard queue with memory ****/

queue_type createQueue(int size, long eltsize) {
    int i;

    queue_type queue;
    queue.node = newVector(size, int);
    queue.history = newVector(eltsize, char);
    queue.readptr = 0;
    queue.writeptr = 0;
    queue.size = size;
    queue.curelts = 0;

    for (i = 0; i < eltsize; i++) queue.history[i] = NEVER_IN_QUEUE;
    for (i = 0; i < size; i++) queue.node[i] = 0;
    return queue;
}

void deleteQueue(queue_type *queue) {
    deleteVector(queue->node);
    deleteVector(queue->history);
}

void enQueue(queue_type *queue, int elt) {
    if (queue->history[elt] == IN_QUEUE) return;
    if (queue->curelts == queue->size) fatalError("Queue not large enough!");
    queue->curelts++;
    queue->node[queue->writeptr] = elt;
    queue->writeptr++;
    if (queue->writeptr == queue->size) queue->writeptr = 0;
    queue->history[elt] = IN_QUEUE;
}

void frontQueue(queue_type *queue, int elt) {
    if (queue->history[elt] == IN_QUEUE) return;
    if (queue->readptr == 0)
        queue->readptr = queue->size - 1;
    else
        queue->readptr--;
    if (queue->curelts == queue->size) fatalError("Queue not large enough!");
    queue->curelts++;
    queue->node[queue->readptr] = elt;
    queue->history[elt] = IN_QUEUE;
}

int deQueue(queue_type *queue) {
    int val = queue->node[queue->readptr];
    queue->history[queue->node[queue->readptr]] = WAS_IN_QUEUE;
    queue->readptr++;
    queue->curelts--;
    if (queue->readptr >= queue->size) queue->readptr = 0;
    return val;
}

void displayQueue(int minVerbosity, queue_type *queue) {
    int i;
    for (i = 0; i < queue->size; i++) {
        displayMessage(minVerbosity, "%ld ", queue->node[i]);
        if (i == queue->readptr)
            displayMessage(minVerbosity, "R");
        else
            displayMessage(minVerbosity, " ");
        if (i == queue->writeptr)
            displayMessage(minVerbosity, "W");
        else
            displayMessage(minVerbosity, " ");
        displayMessage(minVerbosity, "   %ld %d", i, queue->history[i]);
        displayMessage(minVerbosity, "\n");
    }
}


/******************
 ** Binary heaps **
 ******************/

heap_type *createHeap(int heapsize, int eltsize) {

    int i;
    declareScalar(heap_type, newHeap);
    newHeap->node = newVector(heapsize, int);
    newHeap->nodeNDX = newVector(eltsize, int);
    newHeap->last = NOT_IN_HEAP;
    newHeap->valueFn = newVector(eltsize, double);
    newHeap->maxsize = heapsize;
    newHeap->maxelts = eltsize;

    for(i = 0; i < eltsize; i++) newHeap->nodeNDX[i] = NOT_IN_HEAP;
    return newHeap;
}

void insertHeap(heap_type *heap, int key, int value) {
    int elt = ++(heap->last);
    if (heap->last >= heap->maxsize) fatalError("Heap not big enough.");
    heap->node[heap->last] = key;
    heap->nodeNDX[key] = heap->last;
    heap->valueFn[key] = value;
    siftUp(heap, elt);
}

int findMinHeap(heap_type *heap) {
    return heap->node[0];
}

void deleteMinHeap(heap_type *heap) {
    if (heap->last < 0) fatalError("Negative heap size!");
    heap->nodeNDX[heap->node[heap->last]] = 0;
    heap->nodeNDX[heap->node[0]] = NOT_IN_HEAP;
    if (heap->last > 0) swap(heap->node[0], heap->node[heap->last]);
    heap->last--;
    if (heap->last >= 0) siftDown(heap, 1);
}

void deleteHeap(heap_type *heap) {
    deleteVector(heap->node);
    deleteVector(heap->nodeNDX);
    deleteVector(heap->valueFn);
    deleteScalar(heap);
}

void decreaseKey(heap_type *heap, int elt, int value) {
    heap->valueFn[heap->node[heap->nodeNDX[elt]]] = value;
    siftUp(heap, heap->nodeNDX[elt]);
}

void increaseKey(heap_type *heap, int elt, int value) {
    heap->valueFn[heap->node[heap->nodeNDX[elt]]] = value;
    siftDown(heap, heap->nodeNDX[elt]);
}

void siftUp(heap_type *heap, int elt) {
    while (elt > 0 && heap->valueFn[heap->node[elt]]
                      < heap->valueFn[heap->node[heapPred(elt)]]) {
        swap(heap->nodeNDX[heap->node[elt]],
             heap->nodeNDX[heap->node[heapPred(elt)]]);
        swap(heap->node[elt], heap->node[heapPred(elt)]);
        elt = heapPred(elt);
    }
}

void siftDown(heap_type *heap, int elt) {
    int tmp;
    while (heapSucc(elt) <= heap->last
           && heap->valueFn[heap->node[elt]]
                > heap->valueFn[heap->node[minChild(heap, elt)]]) {
        tmp = minChild(heap, elt);
        swap(heap->nodeNDX[heap->node[elt]], heap->nodeNDX[heap->node[tmp]]);
        swap(heap->node[elt], heap->node[tmp]);
        elt = tmp;
    }
}

int heapPred(int elt) {
    return (elt - 1) / 2;
}

int heapSucc(int elt) {
    return (elt * 2) + 1;
}

int minChild(heap_type *heap, int elt) {
    if (heapSucc(elt) == heap->last)
        return heap->last;
    if (heap->valueFn[heap->node[heapSucc(elt)]]
            <= heap->valueFn[heap->node[heapSucc(elt) + 1]])
        return heapSucc(elt);
    return heapSucc(elt) + 1;
}

void heapify(heap_type *heap) {
    int i;
    for (i = heapPred(heap->last); i >= 0; i--)
        siftDown(heap, i);
}

void displayHeap(int minVerbosity, heap_type *heap) {
    int i;
    displayMessage(minVerbosity, "HEAP current size, capacity, number of "
            "elements: %d %d %d\n", heap->last, heap->maxsize, heap->maxelts);
    displayMessage(minVerbosity, "HEAP STATUS");
    for (i = 0; i < heap->maxsize; i++) {
        displayMessage(minVerbosity, "\n%d %d", i, heap->node[i]);
        if (i == heap->last) displayMessage(minVerbosity, " LAST");
    }
    displayMessage(minVerbosity, "\nNODE NDX");
    for (i = 0; i < heap->maxelts; i++)
        displayMessage(minVerbosity, "\n%d %d", i, heap->nodeNDX[i]);
    displayMessage(minVerbosity, "\nVALUE FUNCTION");
    for (i = 0; i < heap->maxelts; i++)
        displayMessage(minVerbosity, "\n%d %f", i, heap->valueFn[i]);
}

/***************************
 ** Memory (de)allocation **
 ***************************/


void *allocateScalar(size_t size) {
    void *scalar = malloc(size);
    if (scalar == NULL) fatalError("Unable to allocate memory for a scalar.");
    #ifdef MEMCHECK
    memcheck_numScalars++;
    if (memcheck_numScalars > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG, "Now have %ld scalars.\n", memcheck_numScalars);
    }
    #endif
    return scalar;
}

void *allocateVector(int u, size_t size) {
    void *vector = malloc(u * size);
    if (vector == NULL)
        fatalError("Unable to allocate memory for vector of size %ld.", u);
    #ifdef MEMCHECK
    memcheck_numVectors++;
    if (memcheck_numVectors > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG, "Now have %ld vectors.\n", memcheck_numVectors);
    }
    #endif
    return vector;
}

void **allocateMatrix(int u1, long u2, size_t size) {
    int i;
    void **matrix = malloc(u1 * sizeof(void *));
    if (matrix == NULL) fatalError("Unable to allocate memory for matrix of "
                                   "size %ld x %ld.", u1, u2);
    for (i = 0; i < u1; i++) {
        matrix[i] = malloc(u2 * size);
        if (matrix[i] == NULL)
            fatalError("Unable to allocate memory for matrix of size "
                       "%ld x %ld.", u1, u2);
    }
    #ifdef MEMCHECK
    memcheck_numMatrices++;
    if (memcheck_numMatrices > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG,"Now have %ld matrices.\n", memcheck_numMatrices);
    }
    #endif
    return matrix;
}

void ***allocate3DArray(int u1, long u2, long u3, size_t size) {
    int i, j;
    void ***matrix = malloc(u1 * sizeof(void **));
    if (matrix == NULL)
        fatalError("Unable to allocate 3D array of size "
                   "%ld x %ld x %ld.", u1, u2, u3);
    for (i = 0; i < u1; i++) {
        matrix[i] = malloc(u2 * sizeof(void *));
        if (matrix[i] == NULL)
            fatalError("Unable to allocate 3D array of size "
                       "%ld x %ld x %ld.", u1, u2, u3);
        for (j = 0; j < u2; j++) {
            matrix[i][j] = malloc(u3 * size);
            if (matrix[i][j] == NULL)
                fatalError("Unable to allocate 3D array of size "
                           "%ld x %ld x %ld.", u1, u2, u3);
        }
    }
    #ifdef MEMCHECK
    memcheck_num3DArrays++;
    if (memcheck_num3DArrays > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG,"Now have %ld 3D arrays.\n",memcheck_num3DArrays);
    }
    #endif
    return matrix;
}

void killScalar(void *scalar) {
    free(scalar);
    #ifdef MEMCHECK
    memcheck_numScalars--;
    if (memcheck_numScalars > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG, "Now have %ld scalars.\n", memcheck_numScalars);
    }
    #endif
}

void killVector(void *vector) {
    free(vector);
    #ifdef MEMCHECK
    memcheck_numVectors--;
    if (memcheck_numVectors > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG, "Now have %ld vectors.\n", memcheck_numVectors);
    }
    #endif
}

void killMatrix(void **matrix, int u1) {
    int i;
    for (i = 0; i < u1; i++) free(matrix[i]);
    free(matrix);
    #ifdef MEMCHECK
    memcheck_numMatrices--;
    if (memcheck_numMatrices > MEMCHECK_THRESHOLD) {
        displayMessage(DEBUG,"Now have %ld matrices.\n", memcheck_numMatrices);
    }
    #endif
}

void kill3DArray(void ***array, int u1, long u2) {
    int i, j;
    for (i = 0; i < u1; i++) {
        for (j = 0; j < u2; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
    #ifdef MEMCHECK
    memcheck_num3DArrays--;
    if (memcheck_num3DArrays > MEMCHECK_THRESHOLD) { 
        displayMessage(DEBUG,"Now have %ld 3D arrays.\n",memcheck_num3DArrays);
    }
    #endif
}
