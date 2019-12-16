/*
 * datastructures.h -- This is the header file for implementations of common
 * data structures: singly and doubly linked lists, queues, and binary heaps.
 * Custom memory allocation/deallocation functions are also provided here
 * (wrappers to malloc and free).
 */

#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

/******************
 ** Linked lists **
 ******************/

/***** Singly linked lists *****/

typedef struct linkedListElt_s {
    int value;
    struct linkedListElt_s *next;
} linkedListElt;

typedef struct {
    linkedListElt *head;
    linkedListElt *tail;
    int size;
} linkedList;

linkedList *createLinkedList();
linkedListElt *insertLinkedList(linkedList *list, int value,
                                linkedListElt *after);
void deleteLinkedList(linkedList *list);
void clearLinkedList(linkedList *list);
void displayLinkedList(int minVerbosity, linkedList *list);

/**** Doubly linked lists ****/

typedef struct doublyLinkedListElt_s {
    double value;
    struct doublyLinkedListElt_s *next;
    struct doublyLinkedListElt_s *prev;
} doublyLinkedListElt;

typedef struct {
    doublyLinkedListElt *head;
    doublyLinkedListElt *tail;
    int size;
} doublyLinkedList;

doublyLinkedList *createDoublyLinkedList();
doublyLinkedListElt *insertDoublyLinkedList(doublyLinkedList *list,
                                            double value,
                                            doublyLinkedListElt *after);
void deleteDoublyLinkedList(doublyLinkedList *list);
void deleteDoublyLinkedListElt(doublyLinkedList *list,
                               doublyLinkedListElt *elt);
void displayDoublyLinkedList(int minVerbosity, doublyLinkedList *list);


/************
 ** Queues **
 ************/

/**** Standard queue with memory ****/

enum {
    IN_QUEUE,
    WAS_IN_QUEUE,
    NEVER_IN_QUEUE
};

typedef enum {
    DEQUE,
    FIFO,
    LIFO
} queueDiscipline;

typedef struct {
    int* node;
    char* history;
    int readptr;
    int writeptr;
    int size;
    int curelts;
} queue_type;

queue_type createQueue(int size, long eltsize);
void deleteQueue(queue_type *queue);
void enQueue(queue_type *queue, int elt);
void frontQueue(queue_type *queue, int elt);
int deQueue(queue_type *queue);
void displayQueue(int minVerbosity, queue_type *queue);

/******************
 ** Binary heaps **
 ******************/

#define NOT_IN_HEAP -1

typedef struct {
    int* node;
    int last;
    double* valueFn;
    int* nodeNDX;
    int maxsize;
    int maxelts;
} heap_type;

heap_type *createHeap(int heapsize, int eltsize);
void insertHeap(heap_type *heap, int key, int value) ;
int  findMinHeap(heap_type *heap);
void deleteMinHeap(heap_type *heap);
void deleteHeap(heap_type *heap);
void decreaseKey(heap_type *heap, int elt, int value);
void increaseKey(heap_type *heap, int elt, int value);
void siftUp(heap_type *heap, int elt);
void siftDown(heap_type *heap, int elt);
int  heapPred(int elt);
int  heapSucc(int elt);
int  minChild(heap_type *heap, int elt);
void heapify(heap_type *heap);
void displayHeap(int minVerbosity, heap_type *heap);

/***************************
 ** Memory (de)allocation **
 ***************************/


/* Uncomment this line to enable memory leak checking */
/* #define MEMCHECK */
#define MEMCHECK_THRESHOLD 1000 /* Threshold before reporting data structure
                                   counts for memory leak checking */

void *allocateScalar(size_t size);
void *allocateVector(int u, size_t size);
void **allocateMatrix(int u1, long u2, size_t size);
void ***allocate3DArray(int u1, long u2, long u3, size_t size);
void killScalar(void *scalar);
void killVector(void *vector);
void killMatrix(void **matrix, int u1);
void kill3DArray(void ***array, int u1, long u2);

#define newScalar(y)            (y *)allocateScalar(sizeof(y))
#define newVector(u,y)          (y *)allocateVector(u,sizeof(y))
#define newMatrix(u1,u2,y)      (y **)allocateMatrix(u1,u2,sizeof(y))
#define new3DArray(u1,u2,u3,y)  (y ***)allocate3DArray(u1,u2,u3,sizeof(y))

#define declareScalar(y,S)           y *S = newScalar(y)
#define declareVector(y,V,u)         y *V = newVector(u,y)
#define declareMatrix(y,M,u1,u2)     y **M = newMatrix(u1,u2,y)
#define declare3DArray(y,A,u1,u2,u3) y ***A = new3DArray(u1,u2,u3,y)

#define deleteScalar(y)            killScalar(y)
#define deleteVector(y)            killVector(y)
#define deleteMatrix(y,u1)        killMatrix((void **)y,u1)
#define delete3DArray(y,u1,u2)    kill3DArray((void ***)y,u1,u2)

#ifdef MEMCHECK
int memcheck_numScalars, memcheck_numVectors, memcheck_numMatrices;
int memcheck_num3DArrays;
#endif

#endif
