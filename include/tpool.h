/**
 * Implementation of Thread Pool using Cross-Platform Pthread wrapper found at:
 * https://nachtimwald.com/2019/04/12/thread-pool-in-c/
 */
#ifndef __TPOOL_H__
#define __TPOOL_H__

#include <stdbool.h>
#include <stddef.h>

struct tpool;
typedef struct tpool tpool_t;

typedef void (*thread_func_t)(void *arg);

tpool_t *tpool_create(size_t num);
void tpool_destroy(tpool_t *tm);

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg);
void tpool_wait(tpool_t *tm);

#endif /* __TPOOL_H__ */
