#include "MexWrappers.h"

#ifdef MATLAB_MEX_FILE

// wrapper functions if compiling from Matlab
#include "mex.h"
inline void wrError(const char *errormsg) { mexErrMsgTxt(errormsg); }
inline void* wrCalloc(size_t num, size_t size) { return mxCalloc(num, size); }
inline void* wrMalloc(size_t size) { return mxMalloc(size); }
inline void wrFree(void * ptr) { mxFree(ptr); }

#else

// wrapper functions if compiling from C/C++
void wrError(const char *errormsg)
{
    UNREFERENCED_PARAMETER(errormsg);

    //throw errormsg; 
}

void mexErrMsgTxt(const char *errormsg)
{
    UNREFERENCED_PARAMETER(errormsg);
}

void* wrCalloc(size_t num, size_t size)
{
    return calloc(num, size);
}

void* wrMalloc(size_t size)
{
    return malloc(size);
}

void wrFree(void * ptr)
{
    free(ptr);
}

#endif

// platform independent aligned memory allocation (see also alFree)
void* alMalloc(size_t size, int alignment) {
    const size_t pSize = sizeof(void*), a = alignment - 1;
    void *raw = wrMalloc(size + a + pSize);
    void *aligned = (void*)(((size_t)raw + pSize + a) & ~a);
    *(void**)((size_t)aligned - pSize) = raw;
    return aligned;
}

// platform independent alignned memory de-allocation (see also alMalloc)
void alFree(void* aligned) {
    void* raw = *(void**)((char*)aligned - sizeof(void*));
    wrFree(raw);
}