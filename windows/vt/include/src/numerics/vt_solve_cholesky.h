#pragma once

#include "vt_matrix.h"

namespace vt {
			
	//
	// Cholesky decomposition of a square symmetric positive definite matrix.
	//
	// A: nxn
	//
	// The Cholesky decomposition returns an upper triangular matrix U such that
	// A = U'U
	//
	// The algorithm requires that the matrix be square and symmtric positive definite.
	//
	// U: nxn

	template <class T>
	class CSolveCholesky : public CErrorBase
	{
	public:
		CSolveCholesky() { }
		CSolveCholesky(const CMtx<T> &mtx) { Solve(mtx); }

		HRESULT Solve(const CMtx<T> &mtx)
		{
			HRESULT hr = NOERROR;

			int n = mtx.Rows();

			if(mtx.IsError())
				VT_HR_EXIT( mtx.GetError() );

			if(mtx.Rows() != mtx.Cols())
				VT_HR_EXIT( E_INVALIDARG );

			// Create the matrix U of the correct size (square)
			VT_HR_EXIT( m_U.Create(n, n) );

			int i, j, k;
			// Clear out the lower part of U
			for (i=1; i<n; i++)
				for (j=0; j<i; j++)
					m_U[i][j] = (T)0.0;

			for (i=0; i<n; i++)
				for (j=i; j<n; j++)
				{
					double sum = mtx[i][j];
					for (k=i-1; k>=0; k--)
						sum -= m_U[k][i] * m_U[k][j];
					if (i == j)
					{
						if (sum <= 0.0)
							VT_HR_EXIT( E_INVALIDARG );    // Not pos def 
						m_U[i][i] = (T)sqrt(sum);
					}
					else
						m_U[i][j] = (T)sum / m_U[i][i];
				}

Exit:
				return SetError(hr);
		}

		void Free() { m_U.Free(); }

		const CMtx<T> &U() const { return m_U; }

	private:
		CMtx<T> m_U;
	};
	
	typedef CSolveCholesky<float> CSolveCholeskyf;
	typedef CSolveCholesky<double> CSolveCholeskyd;
}