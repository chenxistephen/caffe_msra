/*******************************************************************************
* Piotr's Image&Video Toolbox      Version 3.00
* Copyright 2012 Piotr Dollar.  [pdollar-at-caltech.edu]
* Please email me if you find bugs, or have suggestions or questions!
* Licensed under the Simplified BSD License [see external/bsd.txt]
*******************************************************************************/
#pragma once

#include <Windows.h>

#ifdef MATLAB_MEX_FILE

// wrapper functions if compiling from Matlab
#include "mex.h"
inline void wrError(const char *errormsg) { mexErrMsgTxt(errormsg); }
inline void* wrCalloc( size_t num, size_t size ) { return mxCalloc(num,size); }
inline void* wrMalloc( size_t size ) { return mxMalloc(size); }
inline void wrFree( void * ptr ) { mxFree(ptr); }

#else

// wrapper functions if compiling from C/C++
void wrError(const char *errormsg);

void mexErrMsgTxt(const char *errormsg);

void* wrCalloc(size_t num, size_t size);

void* wrMalloc(size_t size);

void wrFree(void * ptr);

#endif

// platform independent aligned memory allocation (see also alFree)
void* alMalloc(size_t size, int alignment);

// platform independent alignned memory de-allocation (see also alMalloc)
void alFree(void* aligned);
