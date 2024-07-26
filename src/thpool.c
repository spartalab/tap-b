/* ********************************
 * Author:       Johan Hanssen Seferidis
 * License:	     MIT
 * Description:  Library providing a threading pool where you can add
 *               work. For usage, check the thpool.h file or README.md
 *
 *//** @file thpool.h *//*
 * Karthik & Rishabh: Updated 11/22/2019 to clean up cast warnings and implicit declaration of functions
 * Alex Vermeulen: Changed the file to work for Windows based platforms. this explicitly breaks linux/macos functionals for now
 ********************************/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <process.h>
#include <signal.h>
#include <time.h>
#include <wchar.h>
#include <stringapiset.h>

#include "thpool.h"

#ifdef THPOOL_DEBUG
#define THPOOL_DEBUG 1
#else
#define THPOOL_DEBUG 0
#endif

#if !defined(DISABLE_PRINT) || defined(THPOOL_DEBUG)
#define err(str) fprintf(stderr, str)
#else
#define err(str)
#endif

static volatile int threads_keepalive;
static volatile int threads_on_hold;

/* ========================== STRUCTURES ============================ */

/* Binary semaphore */
typedef struct bsem {
    CRITICAL_SECTION cs;
    HANDLE sem;
    int v;
} bsem;

/* Job */
typedef struct job {
    struct job* prev;
    void (*function)(void* arg);
    void* arg;
} job;

/* Job queue */
typedef struct jobqueue {
    CRITICAL_SECTION rwmutex;
    job* front;
    job* rear;
    bsem* has_jobs;
    int len;
} jobqueue;

/* Thread */
typedef struct thread {
    int id;
    HANDLE thread_handle;
    HANDLE thread_event;
    struct thpool_* thpool_p;
} thread;

/* Threadpool */
typedef struct thpool_ {
    thread** threads;
    volatile int num_threads_alive;
    volatile int num_threads_working;
    CRITICAL_SECTION thcount_lock;
    HANDLE threads_all_idle;
    jobqueue jobqueue;
} thpool_;

/* ========================== PROTOTYPES ============================ */

static int thread_init(thpool_* thpool_p, thread** thread_p, int id);
static unsigned __stdcall thread_do(void* param);
static void thread_destroy(thread* thread_p);

static int jobqueue_init(jobqueue* jobqueue_p);
static void jobqueue_clear(jobqueue* jobqueue_p);
static void jobqueue_push(jobqueue* jobqueue_p, job* newjob_p);
static job* jobqueue_pull(jobqueue* jobqueue_p);
static void jobqueue_destroy(jobqueue* jobqueue_p);

static void bsem_init(bsem* bsem_p, int value);
static void bsem_post(bsem* bsem_p);
static void bsem_post_all(bsem* bsem_p);
static void bsem_wait(bsem* bsem_p);

/* ========================== THREADPOOL ============================ */

threadpool thpool_init(int num_threads) {
    threads_on_hold = 0;
    threads_keepalive = 1;

    if (num_threads < 0) {
        num_threads = 0;
    }

    /* Make new thread pool */
    thpool_* thpool_p = (thpool_*)malloc(sizeof(thpool_));
    if (thpool_p == NULL) {
        err("thpool_init(): Could not allocate memory for thread pool\n");
        return NULL;
    }
    thpool_p->num_threads_alive = 0;
    thpool_p->num_threads_working = 0;

    /* Initialise the job queue */
    if (jobqueue_init(&thpool_p->jobqueue) == -1) {
        err("thpool_init(): Could not allocate memory for job queue\n");
        free(thpool_p);
        return NULL;
    }

    /* Make threads in pool */
    thpool_p->threads = (thread**)malloc(num_threads * sizeof(thread*));
    if (thpool_p->threads == NULL) {
        err("thpool_init(): Could not allocate memory for threads\n");
        jobqueue_destroy(&thpool_p->jobqueue);
        free(thpool_p);
        return NULL;
    }

    /* Thread init */
    InitializeCriticalSection(&thpool_p->thcount_lock);
    thpool_p->threads_all_idle = CreateEvent(NULL, TRUE, FALSE, NULL);

    for (int n = 0; n < num_threads; n++) {
        if (thread_init(thpool_p, &thpool_p->threads[n], n) == -1) {
            // Cleanup and return on error
            thpool_destroy(thpool_p);
            return NULL;
        }
#if THPOOL_DEBUG
        printf("THPOOL_DEBUG: Created thread %d in pool \n", n);
#endif
    }
    
    /* Wait for threads to initialize */
    while (thpool_p->num_threads_alive != num_threads) {
        Sleep(1);
    }

    return thpool_p;
}

/* Add work to the thread pool */
int thpool_add_work(threadpool thpool_p, void (*function_p)(void*), void* arg_p) {
    job* newjob = (job*)malloc(sizeof(job));
    if (newjob == NULL) {
        err("thpool_add_work(): Could not allocate memory for new job\n");
        return -1;
    }
    /* add function and argument */
    newjob->function = function_p;
    newjob->arg = arg_p;
    /* add job to queue */
    jobqueue_push(&thpool_p->jobqueue, newjob);

    return 0;
}

/* Wait until all jobs have finished */
void thpool_wait(threadpool thpool_p) {
    EnterCriticalSection(&thpool_p->thcount_lock);
    while (thpool_p->jobqueue.len || thpool_p->num_threads_working) {
        WaitForSingleObject(thpool_p->threads_all_idle, INFINITE);
    }
    LeaveCriticalSection(&thpool_p->thcount_lock);
}

/* Destroy the threadpool */
void thpool_destroy(threadpool thpool_p) {
    if (thpool_p == NULL) return;

    volatile int threads_total = thpool_p->num_threads_alive;

    /* End each thread 's infinite loop */
    threads_keepalive = 0;
    
    /* Give one second to kill idle threads */
    double TIMEOUT = 1.0;
    time_t start, end;
    double tpassed = 0.0;
    time(&start);
    while (tpassed < TIMEOUT && thpool_p->num_threads_alive) {
        bsem_post_all(thpool_p->jobqueue.has_jobs);
        time(&end);
        tpassed = difftime(end, start);
    }

    /* Poll remaining threads */
    while (thpool_p->num_threads_alive) {
        bsem_post_all(thpool_p->jobqueue.has_jobs);
        Sleep(1000);
    }

    /* Job queue cleanup */
    jobqueue_destroy(&thpool_p->jobqueue);

    /* Deallocs */
    for (int n = 0; n < threads_total; n++) {
        thread_destroy(thpool_p->threads[n]);
    }
    DeleteCriticalSection(&thpool_p->thcount_lock);
    CloseHandle(thpool_p->threads_all_idle);
    free(thpool_p->threads);
    free(thpool_p);
}

/* Pause all threads in threadpool */
void thpool_pause(threadpool thpool_p) {
    (void)thpool_p; // Silence the unused parameter warning
    threads_on_hold = 1;
}

/* Resume all threads in threadpool */
void thpool_resume(threadpool thpool_p) {
    // resuming a single threadpool hasn't been
    // implemented yet, meanwhile this supresses
    // the warnings    
    (void)thpool_p; // Silence the unused parameter warning
    threads_on_hold = 0;
}

int thpool_num_threads_working(threadpool thpool_p) {
    return thpool_p->num_threads_working;
}

/* ========================== JOBQUEUE ============================= */

static int jobqueue_init(jobqueue* jobqueue_p) {
    jobqueue_p->front = NULL;
    jobqueue_p->rear = NULL;
    jobqueue_p->len = 0;

    InitializeCriticalSection(&jobqueue_p->rwmutex);
    jobqueue_p->has_jobs = (bsem*)malloc(sizeof(bsem));
    if (jobqueue_p->has_jobs == NULL) {
        err("jobqueue_init(): Could not allocate memory for binary semaphore\n");
        return -1;
    }
    bsem_init(jobqueue_p->has_jobs, 0);

    return 0;
}

static void jobqueue_clear(jobqueue* jobqueue_p) {
    while (jobqueue_p->len) {
        job* job_p = jobqueue_pull(jobqueue_p);
        if (job_p) {
            free(job_p);
        }
    }
}

static void jobqueue_push(jobqueue* jobqueue_p, job* newjob_p) {
    EnterCriticalSection(&jobqueue_p->rwmutex);

    newjob_p->prev = NULL;
    if (jobqueue_p->len == 0) {
        jobqueue_p->front = newjob_p;
        jobqueue_p->rear = newjob_p;
    } else {
        jobqueue_p->rear->prev = newjob_p;
        jobqueue_p->rear = newjob_p;
    }
    jobqueue_p->len++;

    bsem_post(jobqueue_p->has_jobs);

    LeaveCriticalSection(&jobqueue_p->rwmutex);
}

static job* jobqueue_pull(jobqueue* jobqueue_p) {
    EnterCriticalSection(&jobqueue_p->rwmutex);

    job* job_p = jobqueue_p->front;
    if (jobqueue_p->len == 1) {
        jobqueue_p->front = NULL;
        jobqueue_p->rear = NULL;
    } else {
        jobqueue_p->front = jobqueue_p->front->prev;
    }
    jobqueue_p->len--;

    LeaveCriticalSection(&jobqueue_p->rwmutex);

    return job_p;
}

static void jobqueue_destroy(jobqueue* jobqueue_p) {
    jobqueue_clear(jobqueue_p);
    DeleteCriticalSection(&jobqueue_p->rwmutex);
    free(jobqueue_p->has_jobs);
}

/* ========================== BINARY SEMAPHORE ===================== */

static void bsem_init(bsem* bsem_p, int value) {
    InitializeCriticalSection(&bsem_p->cs);
    bsem_p->sem = CreateSemaphore(NULL, value, INT_MAX, NULL);
    bsem_p->v = value;
}

static void bsem_post(bsem* bsem_p) {
    EnterCriticalSection(&bsem_p->cs);
    ReleaseSemaphore(bsem_p->sem, 1, NULL);
    bsem_p->v++;
    LeaveCriticalSection(&bsem_p->cs);
}

static void bsem_post_all(bsem* bsem_p) {
    EnterCriticalSection(&bsem_p->cs);
    while (bsem_p->v > 0) {
        ReleaseSemaphore(bsem_p->sem, 1, NULL);
        bsem_p->v--;
    }
    LeaveCriticalSection(&bsem_p->cs);
}

static void bsem_wait(bsem* bsem_p) {
    WaitForSingleObject(bsem_p->sem, INFINITE);
    EnterCriticalSection(&bsem_p->cs);
    bsem_p->v--;
    LeaveCriticalSection(&bsem_p->cs);
}

/* ========================== THREAD MANAGEMENT ===================== */

static int thread_init(thpool_* thpool_p, thread** thread_p, int id) {
    *thread_p = (thread*)malloc(sizeof(thread));
    if (*thread_p == NULL) {
        err("thread_init(): Could not allocate memory for thread\n");
        return -1;
    }
    (*thread_p)->id = id;
    (*thread_p)->thpool_p = thpool_p;
    (*thread_p)->thread_event = CreateEvent(NULL, TRUE, FALSE, NULL);
    if ((*thread_p)->thread_event == NULL) {
        err("thread_init(): Could not create thread event\n");
        free(*thread_p);
        return -1;
    }

    (*thread_p)->thread_handle = (HANDLE)_beginthreadex(NULL, 0, thread_do, *thread_p, 0, NULL);
    if ((*thread_p)->thread_handle == NULL) {
        err("thread_init(): Could not create thread\n");
        CloseHandle((*thread_p)->thread_event);
        free(*thread_p);
        return -1;
    }

    return 0;
}

static void thread_destroy(thread* thread_p) {
    WaitForSingleObject(thread_p->thread_handle, INFINITE);
    CloseHandle(thread_p->thread_handle);
    CloseHandle(thread_p->thread_event);
    free(thread_p);
}

static unsigned __stdcall thread_do(void* param) {
    thread* thread_p = (thread*)param;

    wchar_t thread_name[128];
    mbstowcs(thread_name, "thread-pool-", sizeof(thread_name) / sizeof(thread_name[0]));
    swprintf(thread_name + wcslen(thread_name), sizeof(thread_name) / sizeof(thread_name[0]) - wcslen(thread_name), L"%d", thread_p->id);

    SetThreadDescription(GetCurrentThread(), thread_name);

    thpool_* thpool_p = thread_p->thpool_p;

    while (threads_keepalive) {
        bsem_wait(thpool_p->jobqueue.has_jobs);

        if (threads_keepalive) {
            EnterCriticalSection(&thpool_p->thcount_lock);
            thpool_p->num_threads_working++;
            LeaveCriticalSection(&thpool_p->thcount_lock);

            job* job_p = jobqueue_pull(&thpool_p->jobqueue);
            if (job_p) {
                job_p->function(job_p->arg);
                free(job_p);
            }

            EnterCriticalSection(&thpool_p->thcount_lock);
            thpool_p->num_threads_working--;
            if (!thpool_p->num_threads_working) {
                SetEvent(thpool_p->threads_all_idle);
            }
            LeaveCriticalSection(&thpool_p->thcount_lock);
        }
    }

    EnterCriticalSection(&thpool_p->thcount_lock);
    thpool_p->num_threads_alive--;
    LeaveCriticalSection(&thpool_p->thcount_lock);

    return 0;
}
