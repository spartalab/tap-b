#include "utils.h"

/*
waitForKey pauses while the user presses a key.
*/
void waitForKey() {
    getchar();
}

/*
SWAP is a generic swap function using void pointers to exchange arbitrary
memory.
*/
void SWAP(void* a, void* b, int size) {
  void* c = malloc(size);
  memcpy(c, a, size);
  memcpy(a, b, size);
  memcpy(b, c, size);
  free(c);
}

/*
Wrapper for fopen 
*/
FILE *openFile(const char *filename, const char *access) {
    FILE *handle = fopen(filename, access);
    if (handle == NULL) fatalError("File %s not found", filename);
    return handle;
}

void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    unsigned int result = fread(ptr, size, nmemb, stream);
    if (result != nmemb) fatalError("Error reading from file. Expecting %d bytes, read %d bytes\n", nmemb, result);
}

/*
updateElapsedTime helps with timing, adding an increment to elapsedTime.  Note
that this function does not reset startTime (in order to allow for calculations
to be selectively excluded from timing, e.g. when generating logs or debug
data).  You need to separately set startTime once you want to begin recording.
*/
double updateElapsedTime(clock_t startTime, double *elapsedTime) {
    *elapsedTime += ((double) (clock() - startTime)) / CLOCKS_PER_SEC;
    return *elapsedTime;
}

/*********************
 ** Status messages **
 *********************/

/*
displayMessage is a general-purpose printing function.  If the (global)
verbosity variable is high enough (exceeds the minVerbosity argument), then
prints the message indicated.  From least to greatest, verbosity levels are:
    NOTHING
    LOW_NOTIFICATIONS
    MEDIUM_NOTIFICATIONS
    FULL_NOTIFICATIONS
    DEBUG
    FULL_DEBUG
Output at the DEBUG and FULL_DEBUG levels is not printed to the screen; rather,
if DEBUG_MODE is enabled (by defining the appropriate preprocessor macro) these
messages will be written to the debug log (typically debug.txt)
*/
void displayMessage(int minVerbosity, const char *format, ...) {
    va_list message;
    if (verbosity < minVerbosity) return;
    if (minVerbosity < DEBUG) {
        va_start(message, format);
        vprintf(format, message);
        va_end(message);
        fflush(stdout);
    }
    #ifdef DEBUG_MODE
        va_start(message, format);
        vfprintf(debugFile, format, message);
        va_end(message);
        fflush(debugFile);
    #endif
}

/*
fatalError is a special version of displayMessage which further terminates the
program with the EXIT_FAILURE return code.
*/
void fatalError(const char *format, ...) {
    va_list message;
    va_start(message, format);
    printf("Fatal error: ");
    vprintf(format, message);
    va_end(message);
    printf("\n");
    fflush(stdout);
    #ifdef DEBUG_MODE
        va_start(message, format);
        fprintf(debugFile, "Fatal error: ");
        vfprintf(debugFile, format, message);
        va_end(message);
        fprintf(debugFile, "\n");
        fflush(debugFile);
    #endif
    if (PAUSE_ON_ERROR == TRUE) waitForKey();
    #ifdef DEBUG_MODE
        fclose(debugFile);
    #endif
    exit(EXIT_FAILURE);
}

/*
warning is a special version of displayMessage; by setting
the PAUSE_ON_WARNING macro this can suspend execution of the program
until the user intervenes.  If PAUSE_ON_WARNING is set to FALSE this
is identical to displayMessage except that the text "Warning: " is printed
first.
*/
void warning(int minVerbosity, const char *format, ...) {
    va_list message;
    if (verbosity < minVerbosity) return;
    if (verbosity < DEBUG) {
        va_start(message, format);
        printf("Warning: ");
        vprintf(format, message);
        va_end(message);
        fflush(stdout);
    }
    #ifdef DEBUG_MODE
        va_start(message, format);
        fprintf(debugFile, "Warning: ");
        vfprintf(debugFile, format, message);
        va_end(message);
        fflush(debugFile);
    #endif
    if (PAUSE_ON_WARNING == TRUE) waitForKey();
}

