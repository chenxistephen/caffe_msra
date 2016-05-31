#pragma once

#include "vt_solve_qr.h"

namespace vt {
	
	//
	// Performs the schur decomposition on a real square matrix A
	// giving matrices U and R such that
	//
	//          U' * A * U = R
	//
	//  where U is orthogonal and R is upper-right triangular
	//
	// http://mathworld.wolfram.com/SchurDecomposition.html
	//

	static const double SCHUR_ZERO_TOL = 1e-16;
	static const int SCHUR_MAX_ITER = 100;

	template <class T>
	class CSolveSchur : public CErrorBase
	{
	public:
		CSolveSchur() {}
		CSolveSchur(const CMtx<T> &mtx) { Solve(mtx); }

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::Solve
		//
		//  Synopsis:  Decomposes a square matrix into U' * A * U = R
		//             where U is orthogonal and R is upper-right triangular
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT Solve(const CMtx<T> &mtx)
		{
			VT_HR_BEGIN();

			ClearError();

			VT_HR_EXIT(mtx.GetError());

			if(mtx.Rows()!=mtx.Cols())
				VT_HR_EXIT(E_INVALIDARG);

			int N = mtx.Rows();
			if(N<1)
				VT_HR_EXIT(E_INVALIDARG);

			VT_HR_EXIT(m_U.Create(N, N));
			m_U.MakeI();

			m_R = mtx;
			VT_HR_EXIT(m_R.GetError());

			VT_HR_EXIT(Hessenberg(m_U, m_R));

			Truncate(m_R, 0, N-1);

			VT_HR_EXIT(QRIteration(m_U, m_R, 0, N-1));

			Truncate(m_R, 0, N-1);

			VT_HR_EXIT_LABEL();

			return SetError(hr);
		}

		void Free() { m_U.Free(); m_R.Free(); }

		const CMtx<T> &U() const { return m_U; }
		const CMtx<T> &R() const { return m_R; }

	private:
		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::Truncate
		//
		//  Synopsis:  
		//
		//  Returns:   none
		//
		//------------------------------------------------------------------------

		void Truncate(CMtx<T> &R, int m, int n) const
		{
			int i,j;
			for(i = m+2; i<=n; ++i)
				for(j = m; j<=i-2; ++j)
					R(i, j) = (T)(0);
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::QRIterationOne
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT QRIterationOne(CMtx<T> &U, CMtx<T> &R, int m, int n) const
		{
			VT_HR_BEGIN();

			CSolveQR<T> Rtempqr;
			CVec3<T> v;
			vector<int> vecrange;
			CMtx3x3<T> Q;
			CMtx<T> Rtemp;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(R.GetError());

			// p(t) = t^2 + at + b is the characteristic polynomial 
			//  of the trailing 2x2 matrix.
			Rtemp = R.Extract(n-1, n-1, 2, 2);
			VT_HR_EXIT(Rtemp.GetError());

			VT_HR_EXIT(Rtempqr.Solve(Rtemp));

			T a = -Rtemp.Tr();
			T b =  Rtempqr.Det();

			// first row of M = p(R)
			v[0] = R(m,m)*R(m,m) + R(m,m+1)*R(m+1,m) + a*R(m,m)+b;
			v[1] = R(m+1,m)*R(m,m) + R(m+1,m+1)*R(m+1,m) + a*R(m+1,m);
			v[2] = R(m+2,m+1)*R(m+1,m);

			T norm_v = v.Magnitude();

			if(v[0]>=(T)0)
				v[0] += norm_v;
			else
				v[0] -= norm_v;

			norm_v = v.Magnitude();

			if(norm_v > 0)
				v /= norm_v;
			else 
				v[0] = 1; // shoud this ever happen...?

			Q.MakeI();
			Q = Q - ((v ^ v) * (T)2); // subtract outer product

			VT_HR_EXIT(vecrange.resize(3));

			vecrange[0] = m;
			vecrange[1] = m+1;
			vecrange[2] = m+2;

			VT_HR_EXIT(Rotate(U, R, Q, vecrange));

			// now R has the form:
			//    * * * * * *
			//    * * * * * *
			//    x * * * * *
			//    x x * * * *
			//    0 0 0 * * *
			//    0 0 0 0 * *
			// so we have to chase away the unwanted non-zeros (marked x).
			VT_HR_EXIT(Hessenberg(U, R));

			VT_HR_END();	
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::Rotate
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT Rotate(CMtx<T> &U, CMtx<T> &R, const CMtx<T> &Q, const vector<int> &vecrange) const
		{
			HRESULT hr = NOERROR;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(R.GetError());
			VT_HR_EXIT(Q.GetError());

			R.SetRows(vecrange, Q.T() * R.GetRows(vecrange));
			VT_HR_EXIT(R.GetError());
			R.SetCols(vecrange, R.GetCols(vecrange) * Q);
			VT_HR_EXIT(R.GetError());
			U.SetCols(vecrange, U.GetCols(vecrange) * Q);
			VT_HR_EXIT(U.GetError());

Exit:
			return hr;
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::Hessenberg
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT Hessenberg(CMtx<T> &U, CMtx<T> &R) const
		{
			VT_HR_BEGIN();

			vector<int> vecrange;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(R.GetError());

			int N = R.Rows();
			int i,j;

			VT_HR_EXIT(vecrange.resize(2));

			for(j = 0; j < N-2; j++) 
			{
				// zap below R(j+1,j)
				for(i = N-1; i >= (j+2); i--) 
				{
					// zap below R(i,j)
					T c = R(j+1, j);
					T s = R(i, j);
					T r = VtHypot(c, s);

					if(r>0) 
					{
						c /= r;
						s /= r;

						CMtx2x2<T> Rot;
						Rot(0,0) = c; 
						Rot(0,1) = -s;
						Rot(1,0) = s; 
						Rot(1,1) =  c;

						vecrange[0] = j+1;
						vecrange[1] = i;

						VT_HR_EXIT(Rotate(U, R, Rot, vecrange));
					}
				}
			}

			VT_HR_END();
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSchur<T>::QRIteration
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT QRIteration(CMtx<T> &U, CMtx<T> &R, int m, int n) const
		{
			HRESULT hr = NOERROR;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(R.GetError());

			if(m==n)
			{
				// 1x1 matrix. nothing to do.
				goto Exit;
			}

			if(m+1==n) 
			{
				// 2x2 matrix.

				// make it equal to
				//   a b
				//   c a
				T theta = atan2(R(n,n) - R(n-1,n-1), R(n,n-1) + R(n-1,n)) / (T)(2);
				T c = cos(theta);
				T s = sin(theta);

				CMtx2x2<T> Rot;
				Rot(0,0) = c;
				Rot(0,1) = -s;
				Rot(1,0) = s;
				Rot(1,1) = c;

				vector<int> vecrange;
				VT_HR_EXIT(vecrange.resize(2));

				vecrange[0] = n-1;
				vecrange[1] = n;

				VT_HR_EXIT(Rotate(U, R, Rot, vecrange));

				{
					// make sure they are exactly equal.
					T mean = (R(n, n) + R(n-1, n-1)) / (T)(2);
					R(n-1, n-1) = mean;
					R(n, n) = mean;
				}

				if(R(n-1, n) * R(n, n-1) >= (T)(0)) 
				{
					// we can solve
					theta = atan2(sqrt(fabs(R(n, n-1))), sqrt(fabs(R(n-1, n))));
					c = cos(theta);
					s = sin(theta);

					Rot(0,0) = c;
					Rot(0,1) = -s;
					Rot(1,0) = s;
					Rot(1,1) = c;

					VT_HR_EXIT(Rotate(U, R, Rot, vecrange));

					// zap below the diagonal
					R(n, n-1) = (T)(0);
				}

				goto Exit;
			}

			int iteration;
			for (iteration = 0; iteration < SCHUR_MAX_ITER; ++iteration) 
			{
				int k;

				VT_HR_EXIT(QRIterationOne(U, R, m, n));

				Truncate(R, m, n);

				// zap negligible subdiagonal elements.
				for (k = m; k <= n-1; ++k) 
				{
					// look at R(k+1, k)
					//
					//  * *
					//  ? *
					T num = fabs(R(k+1, k));
					T den = fabs(R(k, k)) + fabs(R(k, k+1)) + fabs(R(k+1, k+1));
					if(den > (T)(0)) 
					{
						if (num / den < (T)(SCHUR_ZERO_TOL)) 
						{
							// zap it
							R(k+1, k) = (T)(0);
						}
					}
				}

				//    * * * * * *
				//    * * * * * *
				//    0 0 * * * *
				//    0 0 * * * *
				//    0 0 0 0 * *
				//    0 0 0 0 0 *

				// divide and conquer.
				for(k = m; k <= n-1; ++k)
				{
					if(R(k+1, k) == (T)(0))
					{
						VT_HR_EXIT(QRIteration(U, R, m, k));
						VT_HR_EXIT(QRIteration(U, R, k+1, n));
						goto Exit;
					}
				}
			}

			// failed to converge
			VT_HR_EXIT(E_TOOCOMPLEX);

Exit:
			return hr;
		}

	private:
		CMtx<T> m_U, m_R;
	};

	typedef CSolveSchur<float> CSolveSchurf;
	typedef CSolveSchur<double> CSolveSchurd;
}