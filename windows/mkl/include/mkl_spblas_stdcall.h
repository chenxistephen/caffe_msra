/*******************************************************************************
!   Copyright(C) 2005-2013 Intel Corporation. All Rights Reserved.
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
!    Intel(R) Math Kernel Library (MKL) stdcall interface for Sparse BLAS level 2,3
!    routines
!******************************************************************************/

#ifndef _MKL_SPBLAS_STDCALL_H_
#define _MKL_SPBLAS_STDCALL_H_

#include "mkl_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if defined(MKL_STDCALL)

/*Float*/
/* Sparse BLAS Level2 lower case */
void __stdcall mkl_scsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall mkl_scsrsv(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);
void __stdcall mkl_scsrgemv(char *transa, int transa_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_scsrgemv(char *transa, int transa_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_scsrsymv(char *uplo, int uplo_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_scsrsymv(char *uplo, int uplo_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_scsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_scsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);

void __stdcall mkl_scscmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall mkl_scscsv(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);

void __stdcall mkl_scoomv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x, float *beta, float *y);
void __stdcall mkl_scoosv(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x, float *y);
void __stdcall mkl_scoogemv(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall mkl_cspblas_scoogemv(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall mkl_scoosymv(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall mkl_cspblas_scoosymv(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall mkl_scootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall mkl_cspblas_scootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);

void __stdcall mkl_sdiamv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *x, float *beta, float *y);
void __stdcall mkl_sdiasv(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *x, float *y);
void __stdcall mkl_sdiagemv(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, float *x,  float *y);
void __stdcall mkl_sdiasymv(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, float *x,  float *y);
void __stdcall mkl_sdiatrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, float *x,  float *y);

void __stdcall mkl_sskymv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr, float *x, float *beta, float *y);
void __stdcall mkl_sskysv(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr,  float *x, float *y);

void __stdcall mkl_sbsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall mkl_sbsrsv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);
void __stdcall mkl_sbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_sbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_sbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_sbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_sbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall mkl_cspblas_sbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
/* Sparse BLAS Level3 lower case */

void __stdcall mkl_scsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall mkl_scsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall mkl_scscmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall mkl_scscsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall mkl_scoomm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall mkl_scoosm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall mkl_sdiamm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall mkl_sdiasm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *b, MKL_INT *ldb, float *c, MKL_INT *ldc);

void __stdcall mkl_sskysm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr,  float *b, MKL_INT *ldb, float *c, MKL_INT *ldc);
void __stdcall mkl_sskymm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);

void __stdcall mkl_sbsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall mkl_sbsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

/* Upper case declaration */
/* Sparse BLAS Level2 upper case */
void __stdcall MKL_SCSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall MKL_SCSRSV(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);
void __stdcall MKL_SCSRGEMV(char *transa, int transa_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCSRGEMV(char *transa, int transa_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_SCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_SCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);

void __stdcall MKL_SCSCMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall MKL_SCSCSV(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);

void __stdcall MKL_SCOOMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x, float *beta, float *y);
void __stdcall MKL_SCOOSV(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x, float *y);
void __stdcall MKL_SCOOGEMV(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCOOGEMV(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall MKL_SCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall MKL_SCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, float *x,  float *y);

void __stdcall MKL_SDIAMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *x, float *beta, float *y);
void __stdcall MKL_SDIASV(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *x, float *y);
void __stdcall MKL_SDIAGEMV(char *transa, int transa_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, float *x,  float *y);
void __stdcall MKL_SDIASYMV(char *uplo, int uplo_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, float *x,  float *y);
void __stdcall MKL_SDIATRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, float *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, float *x,  float *y);

void __stdcall MKL_SSKYMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr, float *x, float *beta, float *y);
void __stdcall MKL_SSKYSV(char *transa, int transa_len, MKL_INT *m, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr,  float *x, float *y);

void __stdcall MKL_SBSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *beta, float *y);
void __stdcall MKL_SBSRSV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *x, float *y);
void __stdcall MKL_SBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_SBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_SBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);
void __stdcall MKL_CSPBLAS_SBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, float *a, MKL_INT *ia,  MKL_INT *ja, float *x,  float *y);

/* Sparse BLAS Level3 upper case */

void __stdcall MKL_SCSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall MKL_SCSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall MKL_SCSCMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall MKL_SCSCSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall MKL_SCOOMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall MKL_SCOOSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

void __stdcall MKL_SDIAMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall MKL_SDIASM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, float *b, MKL_INT *ldb, float *c, MKL_INT *ldc);

void __stdcall MKL_SSKYSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr,  float *b, MKL_INT *ldb, float *c, MKL_INT *ldc);
void __stdcall MKL_SSKYMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *pntr, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);

void __stdcall MKL_SBSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb, float *beta, float *c, MKL_INT *ldc);
void __stdcall MKL_SBSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, float *alpha, char *matdescra, int matdescra_len, float  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, float *b, MKL_INT *ldb,  float *c, MKL_INT *ldc);

/*Double*/
/* Sparse BLAS Level2 lower case */
void __stdcall mkl_dcsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall mkl_dcsrsv(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);
void __stdcall mkl_dcsrgemv(char *transa, int transa_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dcsrgemv(char *transa, int transa_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_dcsrsymv(char *uplo, int uplo_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dcsrsymv(char *uplo, int uplo_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_dcsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dcsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);

void __stdcall mkl_dcscmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall mkl_dcscsv(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);

void __stdcall mkl_dcoomv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x, double *beta, double *y);
void __stdcall mkl_dcoosv(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x, double *y);
void __stdcall mkl_dcoogemv(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall mkl_cspblas_dcoogemv(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall mkl_dcoosymv(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall mkl_cspblas_dcoosymv(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall mkl_dcootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall mkl_cspblas_dcootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);

void __stdcall mkl_ddiamv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *x, double *beta, double *y);
void __stdcall mkl_ddiasv(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *x, double *y);
void __stdcall mkl_ddiagemv(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, double *x,  double *y);
void __stdcall mkl_ddiasymv(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, double *x,  double *y);
void __stdcall mkl_ddiatrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, double *x,  double *y);

void __stdcall mkl_dskymv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr, double *x, double *beta, double *y);
void __stdcall mkl_dskysv(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr,  double *x, double *y);

void __stdcall mkl_dbsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall mkl_dbsrsv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);
void __stdcall mkl_dbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_dbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_dbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall mkl_cspblas_dbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
/* Sparse BLAS Level3 lower case */

void __stdcall mkl_dcsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall mkl_dcsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall mkl_dcscmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall mkl_dcscsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall mkl_dcoomm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall mkl_dcoosm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall mkl_ddiamm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall mkl_ddiasm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *b, MKL_INT *ldb, double *c, MKL_INT *ldc);

void __stdcall mkl_dskysm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr,  double *b, MKL_INT *ldb, double *c, MKL_INT *ldc);
void __stdcall mkl_dskymm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);

void __stdcall mkl_dbsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall mkl_dbsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

/* Upper case declaration */
/* Sparse BLAS Level2 upper case */
void __stdcall MKL_DCSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall MKL_DCSRSV(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);
void __stdcall MKL_DCSRGEMV(char *transa, int transa_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCSRGEMV(char *transa, int transa_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_DCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_DCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);

void __stdcall MKL_DCSCMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall MKL_DCSCSV(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);

void __stdcall MKL_DCOOMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x, double *beta, double *y);
void __stdcall MKL_DCOOSV(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x, double *y);
void __stdcall MKL_DCOOGEMV(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCOOGEMV(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall MKL_DCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall MKL_DCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, double *x,  double *y);

void __stdcall MKL_DDIAMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *x, double *beta, double *y);
void __stdcall MKL_DDIASV(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *x, double *y);
void __stdcall MKL_DDIAGEMV(char *transa, int transa_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, double *x,  double *y);
void __stdcall MKL_DDIASYMV(char *uplo, int uplo_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, double *x,  double *y);
void __stdcall MKL_DDIATRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, double *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, double *x,  double *y);

void __stdcall MKL_DSKYMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr, double *x, double *beta, double *y);
void __stdcall MKL_DSKYSV(char *transa, int transa_len, MKL_INT *m, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr,  double *x, double *y);

void __stdcall MKL_DBSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *beta, double *y);
void __stdcall MKL_DBSRSV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *x, double *y);
void __stdcall MKL_DBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_DBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_DBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);
void __stdcall MKL_CSPBLAS_DBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, double *a, MKL_INT *ia,  MKL_INT *ja, double *x,  double *y);

/* Sparse BLAS Level3 upper case */

void __stdcall MKL_DCSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall MKL_DCSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall MKL_DCSCMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall MKL_DCSCSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall MKL_DCOOMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall MKL_DCOOSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

void __stdcall MKL_DDIAMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall MKL_DDIASM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, double *b, MKL_INT *ldb, double *c, MKL_INT *ldc);

void __stdcall MKL_DSKYSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr,  double *b, MKL_INT *ldb, double *c, MKL_INT *ldc);
void __stdcall MKL_DSKYMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *pntr, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);

void __stdcall MKL_DBSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
void __stdcall MKL_DBSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, double *alpha, char *matdescra, int matdescra_len, double  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb,  double *c, MKL_INT *ldc);

/*MKL_Complex8*/
/* Sparse BLAS Level2 lower case */
void __stdcall mkl_ccsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_ccsrsv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall mkl_ccsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_ccsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_ccsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall mkl_ccscmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_ccscsv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);

void __stdcall mkl_ccoomv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_ccoosv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall mkl_ccoogemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccoogemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_ccoosymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccoosymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_ccootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_ccootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall mkl_cdiamv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_cdiasv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall mkl_cdiagemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cdiasymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cdiatrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall mkl_cskymv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_cskysv(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr,  MKL_Complex8 *x, MKL_Complex8 *y);

void __stdcall mkl_cbsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall mkl_cbsrsv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall mkl_cbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_cbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_cbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall mkl_cspblas_cbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
/* Sparse BLAS Level3 lower case */

void __stdcall mkl_ccsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_ccsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall mkl_ccscmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_ccscsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall mkl_ccoomm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_ccoosm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall mkl_cdiamm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_cdiasm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall mkl_cskysm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr,  MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_cskymm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall mkl_cbsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_cbsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

/* Upper case declaration */
/* Sparse BLAS Level2 upper case */
void __stdcall MKL_CCSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CCSRSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall MKL_CCSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall MKL_CCSCMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CCSCSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);

void __stdcall MKL_CCOOMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CCOOSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall MKL_CCOOGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCOOGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall MKL_CDIAMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CDIASV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall MKL_CDIAGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CDIASYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CDIATRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex8 *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, MKL_Complex8 *x,  MKL_Complex8 *y);

void __stdcall MKL_CSKYMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CSKYSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr,  MKL_Complex8 *x, MKL_Complex8 *y);

void __stdcall MKL_CBSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *beta, MKL_Complex8 *y);
void __stdcall MKL_CBSRSV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *x, MKL_Complex8 *y);
void __stdcall MKL_CBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);
void __stdcall MKL_CSPBLAS_CBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex8 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex8 *x,  MKL_Complex8 *y);

/* Sparse BLAS Level3 upper case */

void __stdcall MKL_CCSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CCSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall MKL_CCSCMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CCSCSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall MKL_CCOOMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CCOOSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall MKL_CDIAMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CDIASM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall MKL_CSKYSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr,  MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CSKYMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *pntr, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);

void __stdcall MKL_CBSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb, MKL_Complex8 *beta, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CBSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, MKL_Complex8 *alpha, char *matdescra, int matdescra_len, MKL_Complex8  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex8 *b, MKL_INT *ldb,  MKL_Complex8 *c, MKL_INT *ldc);

/*Float*/
/* Sparse BLAS Level2 lower case */
void __stdcall mkl_zcsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zcsrsv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall mkl_zcsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zcsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zcsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall mkl_zcscmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zcscsv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);

void __stdcall mkl_zcoomv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zcoosv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall mkl_zcoogemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcoogemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zcoosymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcoosymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zcootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zcootrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall mkl_zdiamv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zdiasv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall mkl_zdiagemv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zdiasymv(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zdiatrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall mkl_zskymv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zskysv(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr,  MKL_Complex16 *x, MKL_Complex16 *y);

void __stdcall mkl_zbsrmv(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall mkl_zbsrsv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall mkl_zbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zbsrgemv(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zbsrsymv(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_zbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall mkl_cspblas_zbsrtrsv(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
/* Sparse BLAS Level3 lower case */

void __stdcall mkl_zcsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zcsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall mkl_zcscmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zcscsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall mkl_zcoomm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zcoosm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall mkl_zdiamm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zdiasm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall mkl_zskysm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr,  MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zskymm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall mkl_zbsrmm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zbsrsm(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

/* Upper case declaration */
/* Sparse BLAS Level2 upper case */
void __stdcall MKL_ZCSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZCSRSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall MKL_ZCSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall MKL_ZCSCMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZCSCSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);

void __stdcall MKL_ZCOOMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZCOOSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall MKL_ZCOOGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCOOGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCOOSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZCOOTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *rowind, MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall MKL_ZDIAMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZDIASV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall MKL_ZDIAGEMV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZDIASYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZDIATRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_Complex16 *val, MKL_INT *lval,  MKL_INT  *idiag, MKL_INT *ndiag, MKL_Complex16 *x,  MKL_Complex16 *y);

void __stdcall MKL_ZSKYMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZSKYSV(char *transa, int transa_len, MKL_INT *m, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr,  MKL_Complex16 *x, MKL_Complex16 *y);

void __stdcall MKL_ZBSRMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *k, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *beta, MKL_Complex16 *y);
void __stdcall MKL_ZBSRSV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *x, MKL_Complex16 *y);
void __stdcall MKL_ZBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZBSRGEMV(char *transa, int transa_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZBSRSYMV(char *uplo, int uplo_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_ZBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);
void __stdcall MKL_CSPBLAS_ZBSRTRSV(char *uplo, int uplo_len, char *transa, int transa_len, char *diag, int diag_len, MKL_INT *m, MKL_INT *lb, MKL_Complex16 *a, MKL_INT *ia,  MKL_INT *ja, MKL_Complex16 *x,  MKL_Complex16 *y);

/* Sparse BLAS Level3 upper case */

void __stdcall MKL_ZCSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZCSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall MKL_ZCSCMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZCSCSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall MKL_ZCOOMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZCOOSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *rowind,  MKL_INT *colind, MKL_INT *nnz, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall MKL_ZDIAMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZDIASM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *lval, MKL_INT *idiag,  MKL_INT *ndiag, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall MKL_ZSKYSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr,  MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZSKYMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *pntr, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);

void __stdcall MKL_ZBSRMM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb, MKL_Complex16 *beta, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZBSRSM(char *transa, int transa_len, MKL_INT *m, MKL_INT *n, MKL_INT *lb, MKL_Complex16 *alpha, char *matdescra, int matdescra_len, MKL_Complex16  *val, MKL_INT *indx,  MKL_INT *pntrb, MKL_INT *pntre, MKL_Complex16 *b, MKL_INT *ldb,  MKL_Complex16 *c, MKL_INT *ldc);

/* Converters lower case*/

void __stdcall mkl_dcsrbsr(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,double *Acsr,MKL_INT * AJ,MKL_INT * AI,double *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall mkl_dcsrcoo(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,double *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall mkl_ddnscsr(MKL_INT *job,MKL_INT *m,MKL_INT *n,double *Adns,MKL_INT *lda,double *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall mkl_dcsrcsc(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall mkl_dcsrdia(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,double *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall mkl_dcsrsky(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0, double *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall mkl_scsrbsr(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,float *Acsr,MKL_INT * AJ,MKL_INT * AI,float *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall mkl_scsrcoo(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,float *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall mkl_sdnscsr(MKL_INT *job,MKL_INT *m,MKL_INT *n,float *Adns,MKL_INT *lda,float *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall mkl_scsrcsc(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall mkl_scsrdia(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,float *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall mkl_scsrsky(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0, float *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall mkl_ccsrbsr(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,MKL_Complex8 *Acsr,MKL_INT * AJ,MKL_INT * AI,MKL_Complex8 *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall mkl_ccsrcoo(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,MKL_Complex8 *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall mkl_cdnscsr(MKL_INT *job,MKL_INT *m,MKL_INT *n,MKL_Complex8 *Adns,MKL_INT *lda,MKL_Complex8 *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall mkl_ccsrcsc(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex8 *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall mkl_ccsrdia(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex8 *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,MKL_Complex8 *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall mkl_ccsrsky(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0, MKL_Complex8 *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall mkl_zcsrbsr(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,MKL_Complex16 *Acsr,MKL_INT * AJ,MKL_INT * AI,MKL_Complex16 *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall mkl_zcsrcoo(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,MKL_Complex16 *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall mkl_zdnscsr(MKL_INT *job,MKL_INT *m,MKL_INT *n,MKL_Complex16 *Adns,MKL_INT *lda,MKL_Complex16 *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall mkl_zcsrcsc(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex16 *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall mkl_zcsrdia(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex16 *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,MKL_Complex16 *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall mkl_zcsrsky(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0, MKL_Complex16 *Asky,MKL_INT * pointers,MKL_INT * info);

/* Converters upper case*/

void __stdcall MKL_DCSRBSR(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,double *Acsr,MKL_INT * AJ,MKL_INT * AI,double *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall MKL_DCSRCOO(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,double *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall MKL_DDNSCSR(MKL_INT *job,MKL_INT *m,MKL_INT *n,double *Adns,MKL_INT *lda,double *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall MKL_DCSRCSC(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall MKL_DCSRDIA(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0,double *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,double *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall MKL_DCSRSKY(MKL_INT * job,MKL_INT * n,double *Acsr,MKL_INT * AJ0,MKL_INT * AI0, double *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall MKL_SCSRBSR(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,float *Acsr,MKL_INT * AJ,MKL_INT * AI,float *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall MKL_SCSRCOO(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,float *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall MKL_SDNSCSR(MKL_INT *job,MKL_INT *m,MKL_INT *n,float *Adns,MKL_INT *lda,float *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall MKL_SCSRCSC(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall MKL_SCSRDIA(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0,float *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,float *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall MKL_SCSRSKY(MKL_INT * job,MKL_INT * n,float *Acsr,MKL_INT * AJ0,MKL_INT * AI0, float *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall MKL_CCSRBSR(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,MKL_Complex8 *Acsr,MKL_INT * AJ,MKL_INT * AI,MKL_Complex8 *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall MKL_CCSRCOO(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,MKL_Complex8 *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall MKL_CDNSCSR(MKL_INT *job,MKL_INT *m,MKL_INT *n,MKL_Complex8 *Adns,MKL_INT *lda,MKL_Complex8 *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall MKL_CCSRCSC(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex8 *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall MKL_CCSRDIA(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex8 *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,MKL_Complex8 *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall MKL_CCSRSKY(MKL_INT * job,MKL_INT * n,MKL_Complex8 *Acsr,MKL_INT * AJ0,MKL_INT * AI0, MKL_Complex8 *Asky,MKL_INT * pointers,MKL_INT * info);

void __stdcall MKL_ZCSRBSR(MKL_INT * job,MKL_INT * m,MKL_INT * mblk,MKL_INT * ldAbsr,MKL_Complex16 *Acsr,MKL_INT * AJ,MKL_INT * AI,MKL_Complex16 *Absr, MKL_INT * AJB, MKL_INT * AIB, MKL_INT * info);
void __stdcall MKL_ZCSRCOO(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJR,MKL_INT * AIR,MKL_INT * nnz,MKL_Complex16 *Acoo, MKL_INT * ir, MKL_INT * jc, MKL_INT * info);
void __stdcall MKL_ZDNSCSR(MKL_INT *job,MKL_INT *m,MKL_INT *n,MKL_Complex16 *Adns,MKL_INT *lda,MKL_Complex16 *Acsr,MKL_INT *AJ,MKL_INT *AI,MKL_INT *info);
void __stdcall MKL_ZCSRCSC(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex16 *Acsc,MKL_INT * AJ1,MKL_INT * AI1,MKL_INT * info);
void __stdcall MKL_ZCSRDIA(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0,MKL_Complex16 *Adia,MKL_INT * ndiag,MKL_INT * distance,MKL_INT * idiag,MKL_Complex16 *Acsr_rem,MKL_INT * AJ0_rem,MKL_INT * AI0_rem,MKL_INT * info);
void __stdcall MKL_ZCSRSKY(MKL_INT * job,MKL_INT * n,MKL_Complex16 *Acsr,MKL_INT * AJ0,MKL_INT * AI0, MKL_Complex16 *Asky,MKL_INT * pointers,MKL_INT * info);


/* Sparse BLAS Level2 (CSR-CSR) lower case */
void __stdcall  mkl_dcsrmultcsr(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall mkl_dcsrmultd(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *ldc);
void __stdcall mkl_dcsradd(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ja, MKL_INT *ia, double *beta, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall mkl_scsrmultcsr(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall mkl_scsrmultd(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *ldc);
void __stdcall mkl_scsradd(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, float *a, MKL_INT *ja, MKL_INT *ia, float *beta, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall mkl_ccsrmultcsr(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall mkl_ccsrmultd(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall mkl_ccsradd(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *beta, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall mkl_zcsrmultcsr(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall mkl_zcsrmultd(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall mkl_zcsradd(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *beta, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);


/* Sparse BLAS Level2 (CSR-CSR) upper case */
void __stdcall  MKL_DCSRMULTCSR(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall MKL_DCSRMULTD(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, double *a, MKL_INT *ja, MKL_INT *ia, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *ldc);
void __stdcall MKL_DCSRADD(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, double *a, MKL_INT *ja, MKL_INT *ia, double *beta, double *b, MKL_INT *jb, MKL_INT *ib, double *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall MKL_SCSRMULTCSR(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall MKL_SCSRMULTD(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, float *a, MKL_INT *ja, MKL_INT *ia, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *ldc);
void __stdcall MKL_SCSRADD(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, float *a, MKL_INT *ja, MKL_INT *ia, float *beta, float *b, MKL_INT *jb, MKL_INT *ib, float *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall MKL_CCSRMULTCSR(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall MKL_CCSRMULTD(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *ldc);
void __stdcall MKL_CCSRADD(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_Complex8 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex8 *beta, MKL_Complex8 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex8 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

void __stdcall MKL_ZCSRMULTCSR(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);
void __stdcall MKL_ZCSRMULTD(char *transa, int transa_len,  MKL_INT *m, MKL_INT *n, MKL_INT *k, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *ldc);
void __stdcall MKL_ZCSRADD(char *transa, int transa_len, MKL_INT *job, MKL_INT *sort, MKL_INT *m, MKL_INT *n, MKL_Complex16 *a, MKL_INT *ja, MKL_INT *ia, MKL_Complex16 *beta, MKL_Complex16 *b, MKL_INT *jb, MKL_INT *ib, MKL_Complex16 *c, MKL_INT *jc, MKL_INT *ic, MKL_INT *nnzmax, MKL_INT *ierr);

#endif /* MKL_STDCALL */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_SPBLAS_STDCALL_H_ */
