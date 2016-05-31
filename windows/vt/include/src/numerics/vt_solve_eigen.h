#pragma once

#include "vt_solve_qr.h"
#include "vt_solve_schur.h"

namespace vt {
	
	HRESULT EigHouseholderReduction(CMtxd &m, CVecd &vDiag, CVecd &vOffDiag);
	HRESULT EigTridiagonalQLImplicit(CVecd &vDiag, CVecd &vOffDiag, CMtxd &mRot);

	// This function computes the eigenvalue and eigenvectors of a real symmetric matrix A
	// This is faster and more convenient than CSolveEigen when doing PCA on covariance matrices
	// and works for larger dimensions. Eigenvalues and vectors are sorted into descending eigenvalue
	// size and eigenvectors are the column vectors of mV.
	HRESULT VtEigenDecomposition(const CMtxd &mA, CMtxd &mV, CVecd &vD);

	//
	// Solve real eigensystem
	// to give matrix of eigenvectors V and diagonal matrix of eigenvalues D
	//
	//  Synopsis:  Solve real valued eigensystem (square matrix)
	//              gives complex matrices D and V where D is the
	//              (square) diagonal matrix of eigenvalues and D
	//              is the corresponding eigenvector (column) matrix 
	//              see http://mathworld.wolfram.com/EigenDecomposition.html

	template <class T>
	class CSolveEigen : public CErrorBase
	{
	public:
		CSolveEigen() {}
		CSolveEigen(const CMtx<T> &mtx) { Solve(mtx); }

		HRESULT Solve(const CMtx<T> &mtx)
		{
			HRESULT hr = NOERROR;
			CSolveSchur<T> sh;
			CMtx< Complex<T> > Ucomplex;
			CVec< Complex<T> > v;
			int N, j;

			ClearError();
			VT_HR_EXIT(mtx.GetError());

			if(mtx.Cols()!=mtx.Rows())
				VT_HR_EXIT(E_INVALIDARG);

			VT_HR_EXIT(sh.Solve(mtx));

			Ucomplex = VtComplexify(sh.U());
			VT_HR_EXIT(Ucomplex.GetError());

			// We just need to compute the matrix of
			// eigenvectors of R. The matrix for A will
			// be U times that matrix.

			N = mtx.Rows();
			VT_HR_EXIT(m_V.Create(N,N));
			m_V = (T)0;
			VT_HR_EXIT(m_D.Create(N,N));
			m_D = (T)0;
			VT_HR_EXIT(v.Create(N));

			j = 0;
			while(j<N)
			{
				if((j == (N-1)) || (sh.R()(j+1,j) == 0))
				{
					// 1x1 block
					m_D(j,j)= sh.R()(j,j);

					Complex<T> lambda = sh.R()(j,j);
					v = Complex<T>(0);
					v[j] = (T)1;
					VT_HR_EXIT(BackSolve(sh.R(), lambda, v, j));
					m_V.SetCol(j, v / v.Magnitude());
					VT_HR_EXIT(m_V.GetError());
					j++;
				}
				else
				{
					// 2x2 block
					//   a b
					//   c a
					// 
					//    eigenvalues = a +/- sqrt(bc)
					//
					//    eigenvectors = [sqrt[b] ; +/- sqrt(c) ]
					//

					Complex<T> a = sh.R()(j,j);
					if(a != sh.R()(j+1,j+1))
					{
						VT_HR_EXIT(E_FAIL);
					}

					Complex<T> b = sh.R()(j, j+1);
					Complex<T> c = sh.R()(j+1, j);
					Complex<T> sqrt_bc = Complex<T>((T)0,(T)1) * VtSqrt(-b*c);
					Complex<T> sqrt_b = VtSqrt(b);
					Complex<T> sqrt_c = VtSqrt(c);

					m_D(j,j)     = a + sqrt_bc;
					m_D(j+1,j+1) = a - sqrt_bc;

					{
						Complex<T> lambda = m_D(j,j);
						v = Complex<T>(0);
						v[j] = sqrt_b;
						v[j+1] = sqrt_c;
						VT_HR_EXIT(BackSolve(sh.R(), lambda, v, j));
						m_V.SetCol(j, v / v.Magnitude());
						VT_HR_EXIT(m_V.GetError());
					}

					{
						Complex<T> lambda = m_D(j+1,j+1);
						v = Complex<T>(0);
						v[j]= sqrt_b;
						v[j+1]= -sqrt_c;
						VT_HR_EXIT(BackSolve(sh.R(), lambda, v, j));
						m_V.SetCol(j+1, v / v.Magnitude());
						VT_HR_EXIT(m_V.GetError());
					}      

					j+= 2;
				}
			}

			m_V = Ucomplex * m_V;

			VT_HR_EXIT(m_V.GetError());

Exit:
			return SetError(hr);
		}

		void Free() { m_V.Free(); m_D.Free(); }

		const CMtx< Complex<T> > &V() const { return m_V; }
		const CMtx< Complex<T> > &D() const { return m_D; }

	private:    
		HRESULT BackSolve(const CMtx<T> &R, Complex<T> lambda, CVec< Complex<T> > &v, int m)
		{
			HRESULT hr = NOERROR;

			int i,k;

			if(m<1)
				goto Exit;

			for(i = m-1; i >= 0; i--)
			{
				if((i>0) && (R(i,i-1)!=(T)0))
				{
					// 2x2 block: solve for v(i-1) and v(i)
					CVec2<Complex<T> > acc;
					acc = Complex<T>(0);

					for(k = i+1; k < R.Cols(); k++)
					{
						acc[0] += R(i-1,k) * v[k];
						acc[1] += R(i  ,k) * v[k];
					}

					CMtx2x2<Complex<T> > tmp;
					tmp(0,0) = R(i-1,i-1) - lambda;
					tmp(0,1) = R(i-1,i);
					tmp(1,0) = R(i,i-1);
					tmp(1,1) = R(i,i) - lambda;

					Complex<T> det_tmp = tmp.Det();
					v[i-1]= (- tmp(1, 1) * acc[0] + tmp(0, 1) * acc[1]) / det_tmp;
					v[i]  = (  tmp(1, 0) * acc[0] - tmp(0, 0) * acc[1]) / det_tmp;
				}
				else
				{
					Complex<T> acc = (T)0;
					for(k = i+1; k<R.Cols(); k++)
					{
						acc += R(i,k) * v[k];
					}

					Complex<T> denom = R(i,i) - lambda;
					if (denom != T(0))
						v[i] = -acc / denom;
					else
						v[i] = T(0); // OK?
				}
			}

Exit:
			return hr;
		}

		CMtx< Complex<T> > m_V, m_D;
	};

	typedef CSolveEigen<float> CSolveEigenf;
	typedef CSolveEigen<double> CSolveEigend;
}