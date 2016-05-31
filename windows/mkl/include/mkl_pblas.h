/*******************************************************************************
!   Copyright(C) 1999-2013 Intel Corporation. All Rights Reserved.
!
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!
!   *Third Party trademarks are the property of their respective owners.
!
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!  Content:
!      Intel(R) Math Kernel Library (MKL) interface for PBLAS routines
!******************************************************************************/

#ifndef _MKL_PBLAS_H_
#define _MKL_PBLAS_H_

#include "mkl_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// PBLAS Level 1 Routines

void    psamax( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdamax( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcamax( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzamax( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    psasum( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdasum( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    pscasum( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdzasum( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    psaxpy( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdaxpy( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcaxpy( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzaxpy( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    picopy( MKL_INT *n, MKL_INT *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, MKL_INT *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pscopy( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdcopy( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pccopy( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzcopy( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    psdot( MKL_INT *n, float *dot, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pddot( MKL_INT *n, double *dot, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pcdotc( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzdotc( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pcdotu( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzdotu( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    psnrm2( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdnrm2( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pscnrm2( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdznrm2( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    psscal( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdscal( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcscal( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzscal( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcsscal( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzdscal( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    psswap( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdswap( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcswap( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzswap( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );


// PBLAS Level 2 Routines

void    psgemv( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdgemv( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcgemv( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzgemv( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    psagemv( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdagemv( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcagemv( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzagemv( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    psger( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdger( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pcgerc( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzgerc( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pcgeru( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzgeru( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pchemv( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzhemv( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pcahemv( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzahemv( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pcher( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzher( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pcher2( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzher2( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pssymv( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdsymv( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    psasymv( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdasymv( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pssyr( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdsyr( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pssyr2( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdsyr2( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    pstrmv( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdtrmv( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pctrmv( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pztrmv( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    psatrmv( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdatrmv( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcatrmv( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzatrmv( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    pstrsv( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdtrsv( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pctrsv( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pztrsv( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );


// PBLAS Level 3 Routines

void    psgemm( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdgemm( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcgemm( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzgemm( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pchemm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzhemm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pcherk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzherk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pcher2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzher2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pssymm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdsymm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsymm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsymm( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pssyrk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdsyrk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsyrk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsyrk( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pssyr2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdsyr2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsyr2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsyr2k( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pstran( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdtran( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pctranu( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztranu( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pctranc( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztranc( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pstrmm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pdtrmm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pctrmm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pztrmm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );

void    pstrsm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pdtrsm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pctrsm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pztrsm( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );

void    psgeadd( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdgeadd( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcgeadd( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzgeadd( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    pstradd( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdtradd( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pctradd( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztradd( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );


/* OTHER NAMING CONVENSIONS FOLLOW */

// PBLAS Level 1 Routines

void    PSAMAX( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDAMAX( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCAMAX( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZAMAX( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSAMAX_( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDAMAX_( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCAMAX_( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZAMAX_( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    psamax_( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdamax_( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcamax_( MKL_INT *n, float *amax, MKL_INT *indx, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzamax_( MKL_INT *n, double *amax, MKL_INT *indx, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSASUM( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDASUM( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSASUM_( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDASUM_( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    psasum_( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdasum_( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSCASUM( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDZASUM( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSCASUM_( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDZASUM_( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pscasum_( MKL_INT *n, float *asum, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdzasum_( MKL_INT *n, double *asum, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSAXPY( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDAXPY( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCAXPY( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAXPY( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSAXPY_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDAXPY_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCAXPY_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAXPY_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psaxpy_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdaxpy_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcaxpy_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzaxpy_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PICOPY( MKL_INT *n, MKL_INT *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, MKL_INT *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSCOPY( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDCOPY( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCCOPY( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZCOPY( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PICOPY_( MKL_INT *n, MKL_INT *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, MKL_INT *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSCOPY_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDCOPY_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCCOPY_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZCOPY_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    picopy_( MKL_INT *n, MKL_INT *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, MKL_INT *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pscopy_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdcopy_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pccopy_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzcopy_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSDOT( MKL_INT *n, float *dot, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDDOT( MKL_INT *n, double *dot, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSDOT_( MKL_INT *n, float *dot, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDDOT_( MKL_INT *n, double *dot, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psdot_( MKL_INT *n, float *dot, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pddot_( MKL_INT *n, double *dot, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PCDOTC( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZDOTC( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCDOTC_( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZDOTC_( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcdotc_( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzdotc_( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PCDOTU( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZDOTU( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCDOTU_( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZDOTU_( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcdotu_( MKL_INT *n, float *dotu, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzdotu_( MKL_INT *n, double *dotu, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSNRM2( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDNRM2( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSCNRM2( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDZNRM2( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSNRM2_( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDNRM2_( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSCNRM2_( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDZNRM2_( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    psnrm2_( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdnrm2_( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pscnrm2_( MKL_INT *n, float *norm2, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdznrm2_( MKL_INT *n, double *norm2, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSSCAL( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDSCAL( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCSCAL( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZSCAL( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCSSCAL( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZDSCAL( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSSCAL_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDSCAL_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCSCAL_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZSCAL_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCSSCAL_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZDSCAL_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    psscal_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdscal_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcscal_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzscal_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pcsscal_( MKL_INT *n, float *a, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pzdscal_( MKL_INT *n, double *a, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSSWAP( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDSWAP( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCSWAP( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZSWAP( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSSWAP_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDSWAP_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCSWAP_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZSWAP_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psswap_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdswap_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcswap_( MKL_INT *n, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzswap_( MKL_INT *n, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );


// PBLAS Level 2 Routines

void    PSGEMV( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDGEMV( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCGEMV( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZGEMV( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSGEMV_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDGEMV_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCGEMV_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZGEMV_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psgemv_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdgemv_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcgemv_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzgemv_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSAGEMV( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDAGEMV( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCAGEMV( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAGEMV( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSAGEMV_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDAGEMV_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCAGEMV_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAGEMV_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psagemv_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdagemv_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcagemv_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzagemv_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSGER( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDGER( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PSGER_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDGER_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    psger_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdger_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PCGERC( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZGERC( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PCGERC_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZGERC_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pcgerc_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzgerc_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PCGERU( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZGERU( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PCGERU_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZGERU_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pcgeru_( MKL_INT *m, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzgeru_( MKL_INT *m, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PCHEMV( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZHEMV( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCHEMV_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZHEMV_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pchemv_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzhemv_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PCAHEMV( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAHEMV( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCAHEMV_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZAHEMV_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcahemv_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzahemv_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PCHER( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZHER( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PCHER_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZHER_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pcher_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzher_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PCHER2( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZHER2( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PCHER2_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PZHER2_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pcher2_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pzher2_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PSSYMV( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDSYMV( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSSYMV_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDSYMV_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pssymv_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdsymv_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSASYMV( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDASYMV( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSASYMV_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDASYMV_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psasymv_( char *uplo, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdasymv_( char *uplo, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSSYR( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDSYR( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PSSYR_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDSYR_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pssyr_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdsyr_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PSSYR2( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDSYR2( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PSSYR2_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    PDSYR2_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pssyr2_( char *uplo, MKL_INT *n, float *alpha, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );
void    pdsyr2_( char *uplo, MKL_INT *n, double *alpha, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca );

void    PSTRMV( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDTRMV( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCTRMV( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZTRMV( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSTRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDTRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCTRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZTRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pstrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdtrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pctrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pztrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );

void    PSATRMV( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDATRMV( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCATRMV( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZATRMV( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PSATRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PDATRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PCATRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    PZATRMV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    psatrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pdatrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pcatrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, float *beta, float *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );
void    pzatrmv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx, double *beta, double *y, MKL_INT *iy, MKL_INT *jy, MKL_INT *descy, MKL_INT *incy );

void    PSTRSV( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDTRSV( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCTRSV( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZTRSV( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PSTRSV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PDTRSV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PCTRSV_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    PZTRSV_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pstrsv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pdtrsv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pctrsv_( char *uplo, char *trans, char *diag, MKL_INT *n, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );
void    pztrsv_( char *uplo, char *trans, char *diag, MKL_INT *n, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *x, MKL_INT *ix, MKL_INT *jx, MKL_INT *descx, MKL_INT *incx );


// PBLAS Level 3 Routines

void    PSGEMM( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDGEMM( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCGEMM( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZGEMM( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSGEMM_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDGEMM_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCGEMM_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZGEMM_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    psgemm_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdgemm_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcgemm_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzgemm_( char *transa, char *transb, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PCHEMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHEMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCHEMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHEMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pchemm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzhemm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PCHERK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHERK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCHERK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHERK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcherk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzherk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PCHER2K( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHER2K( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCHER2K_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZHER2K_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcher2k_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzher2k_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSSYMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCSYMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZSYMM( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSSYMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCSYMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZSYMM_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pssymm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdsymm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsymm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsymm_( char *side, char *uplo, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSSYRK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYRK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCSYRK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZSYRK( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSSYRK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYRK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCSYRK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZSYRK_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pssyrk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdsyrk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsyrk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsyrk_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSSYR2K( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYR2K( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSSYR2K_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDSYR2K_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcsyr2k_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzsyr2k_( char *uplo, char *trans, MKL_INT *n, MKL_INT *k, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSTRAN( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDTRAN( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSTRAN_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDTRAN_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pstran_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdtran_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PCTRANU( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRANU( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCTRANU_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRANU_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pctranu_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztranu_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PCTRANC( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRANC( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCTRANC_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRANC_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pctranc_( MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztranc_( MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSTRMM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PDTRMM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PCTRMM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PZTRMM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PSTRMM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PDTRMM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PCTRMM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PZTRMM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pstrmm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pdtrmm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pctrmm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pztrmm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );

void    PSTRSM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PDTRSM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PCTRSM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PZTRSM( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PSTRSM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PDTRSM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PCTRSM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    PZTRSM_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pstrsm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pdtrsm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pctrsm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );
void    pztrsm_( char *side, char *uplo, char *transa, char *diag, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *b, MKL_INT *ib, MKL_INT *jb, MKL_INT *descb );

void    PSGEADD( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDGEADD( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCGEADD( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZGEADD( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSGEADD_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDGEADD_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCGEADD_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZGEADD_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    psgeadd_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdgeadd_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pcgeadd_( char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pzgeadd_( char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

void    PSTRADD( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDTRADD( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCTRADD( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRADD( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PSTRADD_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PDTRADD_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PCTRADD_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    PZTRADD_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pstradd_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pdtradd_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pctradd_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, float *alpha, float *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, float *beta, float *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );
void    pztradd_( char *uplo, char *trans, MKL_INT *m, MKL_INT *n, double *alpha, double *a, MKL_INT *ia, MKL_INT *ja, MKL_INT *desca, double *beta, double *c, MKL_INT *ic, MKL_INT *jc, MKL_INT *descc );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_PBLAS_H_ */
