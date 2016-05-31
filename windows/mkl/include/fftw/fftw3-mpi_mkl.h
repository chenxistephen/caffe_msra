/******************************************************************************
 * Copyright(C) 2010-2013 Intel Corporation. All Rights Reserved.
 *
 * The source code, information  and  material ("Material") contained herein is
 * owned  by Intel Corporation or its suppliers or licensors, and title to such
 * Material remains  with Intel Corporation  or its suppliers or licensors. The
 * Material  contains proprietary information  of  Intel or  its  suppliers and
 * licensors. The  Material is protected by worldwide copyright laws and treaty
 * provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
 * modified, published, uploaded, posted, transmitted, distributed or disclosed
 * in any way  without Intel's  prior  express written  permission. No  license
 * under  any patent, copyright  or  other intellectual property rights  in the
 * Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
 * implication, inducement,  estoppel or  otherwise.  Any  license  under  such
 * intellectual  property  rights must  be express  and  approved  by  Intel in
 * writing.
 *
 * *Third Party trademarks are the property of their respective owners.
 *
 * Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
 * this  notice or  any other notice embedded  in Materials by Intel or Intel's
 * suppliers or licensors in any way.
 *
 ******************************************************************************
 *
 * Definitions for MPI FFTW3 wrappers to MKL.
 *
 ******************************************************************************
 */

#ifndef FFTW3_MPI_MKL_H
#define FFTW3_MPI_MKL_H

#if defined(MKL_SINGLE)
#define REAL_T float
#define COMPLEX_T fftwf_complex
#define MKL_PREC DFTI_SINGLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_FLOAT(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_FLOAT(name)
#else
#define REAL_T double
#define COMPLEX_T fftw_complex
#define MKL_PREC DFTI_DOUBLE
#define FFTW_MPI_MANGLE(name) FFTW_MPI_MANGLE_DOUBLE(name)
#define FFTW_MANGLE(name) FFTW_MANGLE_DOUBLE(name)
#endif

#include "fftw3-mpi.h"
#include "fftw3_mkl.h"
#include "mkl_cdft.h"

typedef struct registered_plan_s registered_plan;
struct registered_plan_s
{
    fftw_mkl_plan plan;
    registered_plan *next;
    registered_plan *prev;
};

/* Global helper structure */
typedef struct fftw3_mpi_mkl_s fftw3_mpi_mkl_s;
struct fftw3_mpi_mkl_s
{
    registered_plan *registered_plans;
    void (*register_plan)(fftw_mkl_plan p);
    int (*unregister_plan)(fftw_mkl_plan p);
    void (*fftw_delete_plan)(fftw_mkl_plan p);
};
FFTW_EXTERN fftw3_mpi_mkl_s fftw3_mpi_mkl;

#if defined(USING_OPEN_MPI)
#define MKL_comm_c2m(x) MPI_Comm_c2f(x)
#else
#define MKL_comm_c2m(x) x
#endif

#endif /* FFTW3_MPI_MKL_H */

