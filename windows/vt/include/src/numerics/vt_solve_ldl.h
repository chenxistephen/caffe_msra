#pragma once

#include "vt_matrix.h"

namespace vt {

	//+-----------------------------------------------------------------------
	//
	//  Member:    VtInPlaceLDLSolve
	//
	//  Synopsis:  converts the matrix A in-place into LDL form
	//              ignores the upper triangular part of the source matrix
    //  Note:  2011-11-01-szeliski:  added option to not zero matrix on partial failure
	//
	//  Returns:   S_OK if all is well
	//
	//------------------------------------------------------------------------

	template <class T>
	HRESULT VtInPlaceLDLSolve(CMtx<T> &A, bool dontZeroOnError = false)
	{
		VT_HR_BEGIN();

		if(A.IsError())
			VT_HR_EXIT( A.GetError() );
		if(A.Rows() != A.Cols())
			VT_HR_EXIT( E_INVALIDARG );

		// Does the LDL^T decomposition of a matrix 
		int i, j, k;
		int n = A.Rows();

		// Do it row by row 
		for (j=0; j<n; j++)
		{
			// Compute v 
			for (i=0; i<=j-1; i++)
				A[i][j] = A[j][i] * A[i][i];

			// Now, the diagonal element 
			for (k=0; k<=j-1; k++)
				A[j][j] -= A[j][k]*A[k][j];

			// If this is zero, then we are rank deficient 
			if (A[j][j] == 0.0)  // An unlikely event, in floating point 
			{
				if (! dontZeroOnError)
                {
                    A = 0;
                }
				VT_HR_EXIT( E_INVALIDARG );
			}

			// Now, fix up the elements in the j-th column 
			for (i=j+1; i<n; i++)
			{
				T temp = A[i][j];
				for (k=0; k<j; k++)
					temp -= A[i][k]*A[k][j];
				A[i][j] = temp / A[j][j];
			}
		}

		VT_HR_END();
	}

	//+-----------------------------------------------------------------------
	//
	//  Member:    VtInPlaceLDLEquationSolve
	//
	//  Synopsis:  solves the equation Ax = b, where the LDL decomp of A is given.
	//
	//  Returns:   S_OK if all is well, x overwrites b.
	//
	//------------------------------------------------------------------------

	template <class T>
	HRESULT VtInPlaceLDLEquationSolve(const CMtx<T> &LDL, CVec<T> &b)
	{
		VT_HR_BEGIN();

		int n = LDL.Rows();

		if(b.IsError())
			VT_HR_EXIT( b.GetError() );
		if(LDL.IsError())
			VT_HR_EXIT( LDL.GetError() );

		if(n!=LDL.Cols() || b.Size()!=n)
			VT_HR_EXIT( E_INVALIDARG );

		int i, j;

		// First step of back substitution 
		for (i=1; i<n; i++)
			for (j=0; j<i; j++)
				b[i] -= LDL[i][j] * b[j];

		// Second set of back substitution 
		for (i=n-1; i>=0; i--)
		{
			for (j=i+1; j<n; j++)
				b[i] -= LDL[i][j] * b[j];
			b[i] /= LDL[i][i];
		}

		VT_HR_END();
	}

	template <class T>
	HRESULT VtSolveLeastSquares(const CMtx<T> &A, const CVec<T> &b, CVec<T> &x)
	{
		HRESULT hr = NOERROR;
		CMtx<T> mATA;   
		CVec<T> vATb;
		int n;

		VT_HR_EXIT( A.GetError() );
		VT_HR_EXIT( b.GetError() );
		n = A.Cols();
		if(b.Size()!=A.Rows())
			VT_HR_EXIT( E_INVALIDARG );

		VT_HR_EXIT( mATA.Create(n,n) );
		VT_HR_EXIT( vATb.Create(n) );

		int i,j;
		// zero only the lower triangular part
		// we do not use the upper triangular part
		for(i=0; i<n; i++)
			for(j=i; j<n; j++)
				mATA(j,i) = 0;
		vATb = 0;

		int k;
		// compute A'A and A'b
		// accumulate the lower triangular part of mATA
		for(k=0; k<A.Rows(); k++)
		{
			for(i=0; i<n; i++)
			{
				for(j=i; j<n; j++)
					mATA(j,i) += A(k,i) * A(k,j);
				vATb[i] += A(k,i) * b[k];
			}
		}

		// ignores the upper triangular part of the source matrix
		VT_HR_EXIT( VtInPlaceLDLSolve(mATA) );

		VT_HR_EXIT( VtInPlaceLDLEquationSolve(mATA, vATb) );

		x = vATb;
		VT_HR_EXIT( x.GetError() );

Exit:
		return hr;
	}

	//
	// LDL decomposition of a square symmetric (Hermitian) positive definite matrix
	// 
	//
	// A: n x n
	//
	// LDL decomposition produces the following factorization:
	// A = L D L*
	// where L is a lower triangular matrix with unity diagonal elements, D is a diagonal
	// matrix and L' is the Hermitian transpose of L. The result LDL, is a composite
	// matrix with lower triangle elements Lij from L, diagonal elements Dij from D,
	// and upper triangle elements Uij from L*.
	//
	// LDL factorization requires half the computation of Gaussian elimination
	// (LU decomposition), and is always stable. It is more efficient that Cholesky
	// factorization because it avoids computing the square roots of the diagonal elements.
	//
	// The algorithm requires that the input be square and symmtric positive definite.
	//
	// L: n x n
	// D: n x n
	// LDL: n x n
	//

	template <class T>
	class CSolveLDL : public CErrorBase
	{
	public:
		CSolveLDL() { m_bSym = false; }
		CSolveLDL(const CMtx<T> &mtx) { Solve(mtx); }

		HRESULT Solve(const CMtx<T> &A)
		{
			HRESULT hr = NOERROR;

			m_bSym = false;

			m_LDL = A;
			VT_HR_EXIT( m_LDL.GetError() );

			hr = VtInPlaceLDLSolve(m_LDL);

Exit:
			return SetError(hr);
		}

		void Free() { m_LDL.Free(); }

		CVec<T> EquationSolve(const CVec<T> &x) const
		{
			CVec<T> b = x;

			HRESULT hr = VtInPlaceLDLEquationSolve(m_LDL, b);
			if(FAILED(hr))
			{
				SetError(hr);
				b.SetError(hr);
			}
			return b;
		}

		CMtx<T> L() const
		{
			int n = m_LDL.Cols();
			CMtx<T> L(n, n);
			if(!IsError() && !L.IsError())
			{
				int i,j;
				for(i=0; i<n; i++)
					for(j=0; j<n; j++)
					{
						if(i==j)
							L(i,j) = 1;
						else if(j>i)
							L(i,j) = 0;
						else
							L(i,j) = m_LDL(i,j);
					}
			}
			return L;
		}

		CMtx<T> D() const
		{
			int n = m_LDL.Cols();
			CMtx<T> D(n, n);
			if(!IsError() && !D.IsError())
			{
				int i,j;
				for(i=0; i<n; i++)
					for(j=0; j<n; j++)
					{
						if(i==j)
							D(i,j) = m_LDL(i,j);
						else
							D(i,j) = 0;
					}
			}
			return D;
		}

		CMtx<T>& LDL()
		{
			if(!m_bSym)
			{
				// fix up the upper triangular symmetric part since VtInPlaceLDLSolve does
				// not give us this
				int i,j;
				int n = m_LDL.Rows();
				for(i=0; i<n-1; i++)
					for(j=i+1; j<n; j++)
						m_LDL[i][j] = m_LDL[j][i];
				m_bSym = true;
			}
			return m_LDL;
		}

	private:
		bool m_bSym;
		CMtx<T> m_LDL;
	};
	
	typedef CSolveLDL<float> CSolveLDLf;
	typedef CSolveLDL<double> CSolveLDLd;
}