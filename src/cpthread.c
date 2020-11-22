#include "cpthread.h"

#ifdef _WIN32
int pthread_create(pthread_t *thread, pthread_attr_t *attr, void *(*start_routine)(void *), void *arg)
{
    void(attr);

    if (thread == NULL || start_routine == NULL)
        return 1;

    *thread = CreateThread(NULL, 0, start_routine, arg, 0, NULL);
    if (*thread == NULL)
        return 1;
    return 0;
}

int pthread_join(pthread_t thread, void **value_ptr)
{
    (void)value_ptr;
    WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
    return 0;
}

int pthread_kill(pthread_t thread, int signal)
{
    (void) thread;
    return raise(signal);

}

int pthread_detach(pthread_t thread)
{
    CloseHandle(thread);
}

int pthread_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr)
{
    (void)attr;

    if (mutex == NULL)
        return 1;

    InitializeCriticalSection(mutex);
    return 0;
}

int pthread_mutex_unlock(pthread_mutex_t *mutex)
{
    if (mutex == NULL)
        return 1;
    LeaveCriticalSection(mutex);
    return 0;
}

int pthread_cond_init(thread_cond_t *cond, pthread_condattr_t *attr)
{
    (void)attr;
    if (cond == NULL)
        return 1;
    InitializeConditionVariable(cond);
    return 0;
}

int pthread_cond_destroy(thread_cond_t *cond)
{
    /* Windows does not have a destroy for conditionals */
    (void)cond;
    return 0;
}

int pthread_cond_wait(thread_cond_t *cond, pthread_mutex_t *mutex)
{
    if (cond == NULL || mutex == NULL)
        return 1;
    return pthread_cond_timedwait(cond, mutex, NULL)
}

int pthread_cond_timedwait(thread_cond_t *cond, pthread_mutex_t *mutex,
        const struct timespec *abstime)
{
    if (cond == NULL || mutex == NULL)
        return 1;
    if (!SleepConditionVariableCS(cond, mutex, timespec_to_ms(abstime)))
        return 1;
    return 0;
}

int pthread_cond_signal(thread_cond_t *cond)
{
    if (cond == NULL)
        return 1;
    WakeConditionVariable(cond);
    return 0;
}

int pthread_cond_broadcast(thread_cond_t *cond)
{
    if (cond == NULL)
        return 1;
    WakeAllConditionVariable(cond);
    return 0;
}

static DWORD timespec_to_ms(const struct timespec *abstime)
{
    DWORD t;

    if (abstime == NULL)
        return INFINITE;

    t = ((abstime->tv_sec - time(NULL)) * 1000) + (abstime->tv_nsec / 1000000);
    if (t < 0)
        t = 1;
    return t;
}

int pthread_rwlock_init(pthread_rwlock_t *rwlock, const pthread_rwlockattr_t *attr)
{
    (void)attr;
    if (rwlock == NULL)
        return 1;
    InitializeSRWLock(&(rwlock->lock));
    rwlock->exclusive = false;
    return 0;
}

int pthread_rwlock_destroy(pthread_rwlock_t *rwlock)
{
    (void)rwlock;
}

int pthread_rwlock_rdlock(pthread_rwlock_t *rwlock)
{
    if (rwlock == NULL)
        return 1;
    AcquireSRWLockShared(&(rwlock->lock));
}

int pthread_rwlock_tryrdlock(pthread_rwlock_t *rwlock)
{
    if (rwlock == NULL)
        return 1;
    return !TryAcquireSRWLockShared(&(rwlock->lock));
}

int pthread_rwlock_wrlock(pthread_rwlock_t *rwlock)
{
    if (rwlock == NULL)
        return 1;
    AcquireSRWLockExclusive(&(rwlock->lock));
    rwlock->exclusive = true;
}

int pthread_rwlock_trywrlock(pthread_rwlock_t  *rwlock)
{
    BOOLEAN ret;

    if (rwlock == NULL)
        return 1;

    ret = TryAcquireSRWLockExclusive(&(rwlock->lock));
    if (ret)
        rwlock->exclusive = true;
    return ret;
}

int pthread_rwlock_unlock(pthread_rwlock_t *rwlock)
{
    if (rwlock == NULL)
        return 1;

    if (rwlock->exclusive) {
        rwlock->exclusive = false;
        ReleaseSRWLockExclusive(&(rwlock->lock));
    } else {
        ReleaseSRWLockShared(&(rwlock->lock));
    }
}
#endif

#ifdef _WIN32
unsigned int pcthread_get_num_procs()
{
    SYSTEM_INFO sysinfo;

    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
}

#else

#include <unistd.h>
unsigned int pcthread_get_num_procs()
{
    return (unsigned int)sysconf(_SC_NPROCESSORS_ONLN);
}
#endif

void ms_to_timespec(struct timespec *ts, unsigned int ms)
{
    if (ts == NULL)
        return;
    ts->tv_sec = (ms / 1000) + time(NULL);
    ts->tv_nsec = (ms % 1000) * 1000000;
}
