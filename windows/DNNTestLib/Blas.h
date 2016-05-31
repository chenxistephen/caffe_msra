#pragma once

#define USE_BLAS 1

#if USE_BLAS

#define USE_MKL 0

#if USE_MKL

#include "..\\mkl\\include\\mkl.h"

#pragma comment(lib, "..\\mkl\\intel64\\mkl_intel_lp64.lib")
#pragma comment(lib, "..\\mkl\\intel64\\mkl_intel_thread.lib")
#pragma comment(lib, "..\\mkl\\intel64\\mkl_core.lib")
#pragma comment(lib, "..\\mkl\\compiler\\libiomp5md.lib")

#else

#include "..\\OpenBlas\\include\\cblas.h"
#pragma comment (lib, "..\\OpenBlas\\lib\\libopenblas.dll.a")

#endif


#endif

