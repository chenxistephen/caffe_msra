#pragma once

#include "vt_solve_qr.h"
#include <stdlib.h>

namespace vt {

	//
	// Singular value decomposition
	// http://mathworld.wolfram.com/SingularValueDecomposition.html
	//
	// function [U, D, V] = ms_svd(A)
	//
	// A: m x n
	//
	// Let k = min(m, n).
	//
	// U: m x k
	// D: k x k
	// V: n x k
	//

	const double SVD_ZERO_TOL = 2.5e-16;
	const int SVD_MAX_ITER = 1000; /* note this can be changed if AWF_DEBUG */

	template <class T>
	struct SortPair {
		int index;
		T value;
	};

	/// awf debug stuff -- please leave this in until say Oct 2011 in case we get new failure cases
#define AWF_DEBUG_ON 0
#define AWF_DEBUG if (!AWF_DEBUG_ON); else

#if AWF_DEBUG_ON
#undef SVD_MAX_ITER
#define SVD_MAX_ITER 100

	static inline char dbl2char(double d) 
	{
		if (d == 0) return '0';
		if (fabs(d) < 1e-12) return '.';
		if (fabs(d) < 1e-8) return '_';
		return '*';
	}

	template <class T>
	static inline void ShowMatrix(char const* tag, CMtx<T> const& D)
	{
		for(int i = 0; i < D.Rows(); ++i) {
			for(int j = 0; j < D.Cols(); ++j)
				DebuggerPrint("%c", dbl2char(D(i,j)));
			for(int j = 0; j < D.Cols(); ++j)
				DebuggerPrint("%10.5f", D(i,j));
			if (i==0) 
				DebuggerPrint(" \\__%s\n", tag);
			else
				DebuggerPrint("\n");
		}
	}

#define AWF_DEBUG_DUMP(M) ShowMatrix(#M,M)

	template <class T>
	static inline double fronorm(CMtx<T>& A)
	{
		double n = 0;
		for (int i = 0; i < A.Rows(); ++i)
			for (int j = 0; j < A.Cols(); ++j)
				n += double(A(i,j)*A(i,j));
		return n;
	}
#endif

	/// --- 

	//+-----------------------------------------------------------------------
	//
	//  Function:    SVDCompareFunction
	//
	//  Synopsis:  funciton for use with qsort to sort in descending order
	//
	//  Returns:   -1 or 1
	//
	//------------------------------------------------------------------------

	template <class T>
	int 
#if defined(MSRVT_WINDOWS_BUILD)
__cdecl
#endif
    SVDCompareFunction(const void *pc1, const void *pc2)
	{
		return ((SortPair<T> *)pc1)->value > ((SortPair<T> *)pc2)->value ? -1 : 1;
	}

	template <class T>
	class CSolveSVD : public CErrorBase
	{
	public:
		CSolveSVD() {}
		CSolveSVD(const CMtx<T> &mtx) { Solve(mtx); }

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::Solve
		//
		//  Synopsis:  Decomposes a matrix into USV where U V are orthogonal and
		//              D is diagonal
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT Solve(const CMtx<T> &mtx)
		{
			VT_HR_BEGIN();

			ClearError();

			VT_HR_EXIT(mtx.GetError());

			int m = mtx.Rows();
			int n = mtx.Cols();

			if(m==n) 
			{
				CMtx<T> Imm(m,m);
				VT_HR_EXIT(Imm.GetError());

				CMtx<T> Inn(n,n);
				VT_HR_EXIT(Inn.GetError());

				Imm.MakeI();
				Inn.MakeI();

				m_U = Imm;
				m_D = mtx;
				m_V = Inn;
			}
			else if(m > n) 
			{
				// portrait matrix
				CSolveQR<T> Mqr;
				VT_HR_EXIT(Mqr.Solve(mtx));

				// Q: m x n
				// R: n x n
				// p: 1 x n
				CMtx<T> Inn(n,n);
				VT_HR_EXIT(Inn.GetError());

				Inn.MakeI();

				m_U = Mqr.Q();
				m_D = Mqr.R();
				m_V = Inn.GetCols(Mqr.Perm());
			}
			else 
			{
				// landscape matrix
				CSolveQR<T> Mqr;
				VT_HR_EXIT(Mqr.Solve(mtx.T()));

				// Q: n x m
				// R: m x m
				// p: 1 x m
				CMtx<T> Imm(m,m);
				VT_HR_EXIT(Imm.GetError());

				Imm.MakeI();

				m_U = Imm.GetCols(Mqr.Perm());
				m_D = Mqr.R().T();
				m_V = Mqr.Q();
			}

			VT_HR_EXIT(m_U.GetError());
			VT_HR_EXIT(m_D.GetError());
			VT_HR_EXIT(m_V.GetError());

			hr = SquareSVD(m_U, m_D, m_V);

			VT_HR_EXIT_LABEL();

			return SetError(hr);
		}

		void Free() { m_U.Free(); m_D.Free(); m_V.Free(); }

		T D(int i) const { return m_D(i,i); }

		const CMtx<T> &U() const { return m_U; }
		const CMtx<T> &D() const { return m_D; }
		const CMtx<T> &V() const { return m_V; }



		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::PseudoInverse
		//
		//  Synopsis:  Calculate pseudo inverse. threshold determines how small
		//              the singular values must be to be set to zero in forming
		//              the inverse
		//
		//  Returns:   the inverse matrix
		//
		//------------------------------------------------------------------------

		CMtx<T> PseudoInverse(T threshold = 0) const
		{
			int k = m_D.Rows(); // m_D.Cols();
			CMtx<T> Sinv(k, k);
			Sinv.SetError(GetError());

			if(!Sinv.IsError())
			{
				Sinv = (T)0;
				int i;
				for(i=0; i<k; ++i)
					if(fabs(m_D(i, i)) > threshold)
						Sinv(i, i) = (T)1 / m_D(i, i);
					else
						Sinv(i, i) = (T)0;
			}

			return m_V * Sinv * m_U.T();
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::EquationSolve
		//
		//  Synopsis:  solve the equation Ax=b, given b and the svd of A
		//
		//  Returns:   x
		//
		//------------------------------------------------------------------------

		CVec<T> EquationSolve(const CVec<T> &b) const
		{
			int m = m_U.Rows(); // m x k
			int n = m_V.Rows(); // n x k
			int k = m_D.Rows(); // k x k

			// solve A x = b
			//       U D V* x = b
			//       D (V* x) = U* b
			//       V* x = D \ (U* b)
			//       x = V (D \ (U* b))

			int i,j;
			CVec<T> tmp(k);
			tmp.SetError(GetError());
			tmp.SetError(b.GetError());
			if(b.Size()!=m)
				tmp.SetError(E_INVALIDARG);

			if(tmp.IsError())
			{   
				// try to return a vector of the right size
				HRESULT hr = tmp.GetError();
				tmp.Create(n);
				tmp.ClearError();
				tmp.SetError(hr);
				return tmp;
			}

			//tmp = m_U.T() * b;

			for(j=0; j<k; ++j) 
			{
				T acc = (T)0;
				for(i=0; i<m; ++i)
					acc += m_U(i, j) * b(i);
				tmp(j) = acc;
			}

			for(j=0; j<k; ++j)
				if(m_D(j, j) != (T)0)
					tmp(j) /= m_D(j, j);
				else
					tmp(j) = (T)0;

			return m_V * tmp;
		}

		CVec<T> GetBestNullSpaceVector() const // returns v that minimizes ||Av||
		{
			T sv = m_D(0,0);
			int idx;
			int i;
			for(idx=0, i=1; i<m_V.Cols(); i++)
			{
				if( m_D(i,i) < sv )
				{
					idx = i;
					sv = m_D(i,i);
				}
			}
			return m_V.GetCol(idx);
		}


	private:

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::SquareSVD
		//
		//  Synopsis:  Calculate SVD from initialised matrices
		//              D must be an NxN. The matrices of left and right singular
		//              vectors of D are post-multiplied onto the given matrices
		//              U and V.
		//
		//  Returns:   S_OK if nothing went wrong
		//
		//------------------------------------------------------------------------

		HRESULT SquareSVD(CMtx<T> &U, CMtx<T> &D, CMtx<T> &V) const
		{
			VT_HR_BEGIN();

			CMtx<T> Utemp;
			CMtx<T> Vtemp;
			CMtx<T> Dtemp;
			vector<SortPair<T> > v;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(D.GetError());
			VT_HR_EXIT(V.GetError());

			int N = D.Rows();
			if(N != D.Cols() || N != U.Cols() || N!=V.Cols())
				VT_HR_EXIT(E_INVALIDARG); // really should be internal error?

			int i,j,k;

			// reduction to upper Hessenberg form
			for(k = 0; k < N; k++)
			{
				// clear below D(k,k)
				for(i = k+1; i < N; i++)
				{
					// rotate D(k, k) and D(i, k) to zap D(i, k).
					T c = D(k, k);
					T s = D(i, k);
					T r = VtHypot(c, s);

					if(r > (T)0)
					{
						c /= r;
						s /= r;

						// D([k i], :) = [c s; -s c] * D([k i], :);
						GivensRow(D, k, i, c, -s);

						// U(:, [k i]) = U(:, [k i]) * [c -s; s c];
						GivensCol(U, k, i, c, s);
					}
				}

				// clear to the right of D(k, k+1)
				for(j = k+2; j < N; j++) 
				{
					// rotate D(k, k+1) and D(k, j) to zap D(k, j).
					T c = D(k, k+1);
					T s = D(k, j);
					T r = VtHypot(c, s);

					if(r > (T)0) 
					{
						c /= r;
						s /= r;

						// D(:, [k+1 j]) = D(:, [k+1 j]) * [c -s; s c];
						GivensCol(D, k+1, j, c, s);

						// V(:, [k+1 j]) = V(:, [k+1 j]) * [c -s; s c];
						GivensCol(V, k+1, j, c, s);
					}
				}
			}

			Truncate(D, 0, N-1);

			VT_HR_EXIT(QRIteration(U, D, V, 0, N-1));

			// invert negative singular values.
			for(i = 0; i<N; i++) 
			{
				if(D(i, i) < (T)0) 
				{
					D(i, i) = -D(i, i);
					for(j = 0; j<U.Rows(); ++j)
						U(j, i) = -U(j, i);
				}
			}

			// find indexed sort of SVs
			for(i = 0; i < N; i++)
			{
				SortPair<T> sp;
				sp.index = i;
				sp.value = D(i,i);
				VT_HR_EXIT(v.push_back(sp));
			}

			// sort into descending order
			qsort((void *)&v[0], v.size(), sizeof(v[0]), &SVDCompareFunction<T>);

			Utemp.Create(U.Rows(), U.Cols());
			Vtemp.Create(V.Rows(), V.Cols());
			Dtemp.Create(D.Rows(), D.Cols());

			VT_HR_EXIT(Utemp.GetError());
			VT_HR_EXIT(Vtemp.GetError());
			VT_HR_EXIT(Dtemp.GetError());

			Dtemp = (T)0;

			for(i = 0; i < (int)v.size(); i++)
			{
				Utemp.SetCol(i, U.GetCol(v[i].index));
				Vtemp.SetCol(i, V.GetCol(v[i].index));
				Dtemp(i,i)= D(v[i].index, v[i].index);
			}

			U = Utemp;
			V = Vtemp;
			D = Dtemp;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(V.GetError());
			VT_HR_EXIT(D.GetError());

			VT_HR_END();
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::Truncate
		//
		//  Synopsis:  truncate non-zeros that should be zero
		//
		//  Returns:   none
		//
		//------------------------------------------------------------------------

		void Truncate(CMtx<T> &D, int m, int n) const
		{
			int i, j;

			for(i = m+1; i<=n; ++i)
				for(j = m; j<i; ++j)
					D(i, j) = (T)0;
			for(j = m+2; j<=n; ++j)
				for(i = m; i<=j-2; ++i)
					D(i, j) = (T)0;
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::QRIteration
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if all is well
		//
		//------------------------------------------------------------------------

		HRESULT QRIteration(CMtx<T> &U, CMtx<T> &D, CMtx<T> &V, int m, int n) const
		{
			VT_HR_BEGIN();

#if AWF_DEBUG_ON
			CMtx<T> A = U * D * V.T();
#endif

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(D.GetError());
			VT_HR_EXIT(V.GetError());

			int iteration = 0;
			int k;

			while(m<n)
			{
#if AWF_DEBUG_ON
				DebuggerPrint("iter %3d mn = %3d,%3d:\n",iteration, m,n);
				AWF_DEBUG_DUMP(D);
#endif

				if(++iteration >= SVD_MAX_ITER)
				{
					VT_HR_EXIT(E_TOOCOMPLEX);   
				}

				// do one QR iteration.
				VT_HR_EXIT(QRIterationOne(U, D, V, m, n));

				// zap negligible superdiagonal elements.
				for(k = m; k<n; k++) 
				{
					//    a b
					//    0 c
					// if |b| / (|a| + |c|) is small then we zap b.
					T num = fabs(D(k, k+1));
					T den = fabs(D(k, k)) + fabs(D(k+1, k+1));
					if(den > (T)0)
					{
						T rel = num / den;
						if(rel < (T)SVD_ZERO_TOL)
							// zap it.
							D(k, k+1) = (T)0;
					}
				}

				// zap negligible diagonal elements.
				for(k = m; k<=n; k++)
				{
					//  * a 0
					//    b c
					//      *
					// if |b| / (|a| + |c|) is small then we zap b.
					T num = fabs(D(k, k));
					T den;
					if(k == m)
						den = fabs(D(k, k+1));
					else if(k == n)
						den = fabs(D(k-1, k));
					else
						den = fabs(D(k-1, k)) + fabs(D(k, k+1));

					if(den > (T)0) 
					{
						T rel = num / den;
						if (rel < (T)SVD_ZERO_TOL)
							// zap it.
							D(k, k) = (T)0;
					}
				}

				// look for leading or trailing zeros.
				while(m<n && D(m,m+1)==(T)0)  // [> 0]
					m++;                      // [0 *]
				while(n>m && D(n-1,n)==(T)0)  // [* 0]
					n--;                      // [0 ^]

				// If a zero appears on the diagonal, zap its row
				for(k = m; k<n; k++)
					if (D(k,k) == T(0)) {
						// Premultiply pairs of rows (k,i) with Givens to transform:
						//  k-1: [ a b . . . ] 
						//  k:   [ . 0 d . . ]  // 0 is the zero we've just discovered
						//  k+1: [ . . e f . ]
						//  k+2: [ . . . g h ]
						//  k+3: [ . . . . j ]
						// into   
						//  k-1: [ a b . . . ] 
						//  k:   [ . . . D . ] <- these are the rows we rotated, D = d'
						//  k+1: [ . . E F . ] <-
						//  k+2: [ . . . g h ]
						//  k+3: [ . . . . j ]
						// into   
						//  k-1: [ a b . . . ] 
						//  k:   [ . . . . D ] <-
						//  k+1: [ . . E F . ]
						//  k+2: [ . . . G H ] <-
						//  k+3: [ . . . . j ]
						// into   
						//  k-1: [ a b . . . ] 
						//  k:   [ . . . . . ] <-
						//  k+1: [ . . E F . ]
						//  k+2: [ . . . G H ] 
						//  k+3: [ . . . . j ] <-
#if AWF_DEBUG_ON
						DebuggerPrint("ROWZAPPING\n");
						AWF_DEBUG_DUMP(D);
#endif
						for(int i = k+1; i <= n; ++i) {
#if AWF_DEBUG_ON
							DebuggerPrint("ROWZAPPING i=%d, norm = %g\n",i, fronorm(A - U*D*V.T()));
#endif
							// [ c s] [0 d 0] __\ [0 0 .. D]
							// [-s c] [0 e f]   / [0 E .. F]
							// So c d + s e = 0;
							T d = D(k,i);
							if (d == T(0)) {
								// If it's zero, then we're done, fall into the "blocked" case below.
#if AWF_DEBUG_ON
								DebuggerPrint("ROWZAPPING Early bailout\n");
#endif
								break;
							}
							T e = D(i,i);
							T c = e;
							T s = -d;
							T r = VtHypot(c, s); // Nonzero as d != 0
							c /= r;
							s /= r;
							// D([k i], :) = [c s; -s c] * D([k i], :);
							GivensRow(D, k, i, c, -s);
							// U(:, [k i]) = U(:, [k i]) * [c -s; s c];
							GivensCol(U, k, i, c,  s);
#if AWF_DEBUG_ON
							AWF_DEBUG_DUMP(D);
#endif
						}
#if AWF_DEBUG_ON
						DebuggerPrint("ROWZAPPING done. norm = %g\n", fronorm(A - U*D*V.T()));
#endif
					}

					// if a zero appears in the middle, recurse.
					for(k = m+1; k<n; k++)
					{
						// If diag is 
						//    a 0
						//    0 b
						// then matrix is blocked, hurrah.
						// Split into m:k and (k+1):n
						if(D(k,k+1)==(T)0) 
						{
							// recurse
							VT_HR_EXIT(QRIteration(U, D, V, m, k));
							VT_HR_EXIT(QRIteration(U, D, V, k+1, n));
							return NOERROR;
						}
					}
			}

			VT_HR_END();
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::QRIterationOne
		//
		//  Synopsis:  
		//
		//  Returns:   S_OK if all is well
		//
		//------------------------------------------------------------------------

		HRESULT QRIterationOne(CMtx<T> &U, CMtx<T> &D, CMtx<T> &V, int m, int n) const
		{
			VT_HR_BEGIN();

			CVec2<T> q;

			VT_HR_EXIT(U.GetError());
			VT_HR_EXIT(D.GetError());
			VT_HR_EXIT(V.GetError());

			if(n == m+1) 
			{
				// special cases for 2x2.
				//   a b
				//   0 c
				T a = D(m, m);
				T b = D(m, m+1);
				T c = D(m+1, m+1);

				if(b == (T)0)
					goto Exit; // nothing to do

				if(a == (T)0) 
				{
					// [0 b] = 1/r [b -c] [r 0] [0 1]
					// [0 c]       [c  b] [0 0] [1 0]
					T r = VtHypot(b, c); // non-zero since b is non-zero.

					D(m, m) = r;
					D(m, m+1) = (T)0;
					D(m+1, m+1) = (T)0;

					CMtx2x2<T> Mu;
					Mu(0, 0) = b / r; 
					Mu(0, 1) = -c / r;
					Mu(1, 0) = c / r; 
					Mu(1, 1) =  b / r;

					// U = U * [Mu]
					U.Update(0, m, U.Extract( 0, m, -1, 2) * Mu);
					VT_HR_EXIT(U.GetError());

					CMtx2x2<T> Mv;
					Mv(0, 0) = (T)0;
					Mv(0, 1) = (T)1;
					Mv(1, 0) = (T)1;
					Mv(1, 1) = (T)0;

					// V = V * [Mv];
					V.Update(0, m, V.Extract( 0, m, -1, 2) * Mv);
					VT_HR_EXIT(V.GetError());

					goto Exit;
				}

				if (c == (T)0) 
				{
					// [a b] = [1 0] [r 0] 1/r [a -b]^T
					// [0 0]   [0 1] [0 0]     [b  a]
					T r = VtHypot(a, b); // non-zero since a, b are non-zero.

					D(m, m) = r;
					D(m, m+1) = (T)0;
					D(m+1, m+1) = (T)0;

					CMtx2x2<T> Mv;
					Mv(0, 0) = a / r;
					Mv(0, 1) = - b / r;
					Mv(1, 0) = b / r;
					Mv(1, 1) =   a / r;

					// V = V * [v];
					V.Update(0, m, V.Extract( 0, m, -1, 2) * Mv);
					VT_HR_EXIT(V.GetError());

					goto Exit;
				}
			}

			// General case, larger than 2x2:

			// D has this pattern:
			//      * * 0 0 0
			//      0 * * 0 0
			//      0 0 * * 0
			//      0 0 0 * *
			//      0 0 0 0 *

			// so M = D' D has this pattern:
			//      * * 0 0 0
			//      * * * 0 0
			//      0 * * * 0
			//      0 0 * * *
			//      0 0 0 * *

			// calculate shift sigma = M(end, end) [single-shift method].
			T sigma = D(n,n)*D(n,n) + D(n-1,n)*D(n-1,n);

			// awf: we should instead compute the shift from Golub and van Loan (eq 8.3.3)
			{
				// The bottom right corner of D'*D = [a b ; 0 c]
				T a = D(n-1,n-1);
				T b = D(n-1,n);
				T c = D(n,n);
				T r = b*b + c*c;
				T aa = a*a;
				T ab = a*b;
				// Corner of D'D is [aa ab; ab r]
				T d = (r - aa)/2;
				T sign_d = (d > 0) ? T(+1) : T(-1);
				T mu = r + d - sign_d * sqrt(d*d + ab*ab);
				sigma = mu;
#if AWF_DEBUG_ON
				DebuggerPrint("sold = %g, snew = %g\n", oldsigma, sigma);
#endif
			}

			// We would compute the QR factorization of the matrix
			//   M - sigma I
			//
			// and then (QR iteration) form the matrix
			//   RQ = (DQ)' (DQ) - sigma I
			//
			// DQ has this pattern:
			//      * * * * *
			//      * * * * *
			//      0 * * * *
			//      0 0 * * *
			//      0 0 0 * *
			//
			// This has one sub-diagonal line and junk above. Use
			// Givens row operations to clean up below the diagonal.
			// This will automagically clean above the diagonal too.

			// But in the QR factorization, Q has this pattern:
			//      * * * * *
			//      * * * * *
			//      0 * * * *
			//      0 0 * * *
			//      0 0 0 * *
			// and we only need the first column q of Q:
			q[0]= D(m,m)*D(m,m) - sigma;
			q[1]= D(m,m)*D(m,m+1);

			// In fact we have the factorization
			//    [ * * * * * ] = [ * * 0 0 0 ] [ 1 0 0 0 0 ]
			//    [ * * * * * ]   [ * * 0 0 0 ] [ 0 * * * * ]
			//    [ 0 * * * * ]   [ 0 0 1 0 0 ] [ 0 * * * * ]
			//    [ 0 0 * * * ]   [ 0 0 0 1 0 ] [ 0 0 * * * ]
			//    [ 0 0 0 * * ]   [ 0 0 0 0 1 ] [ 0 0 0 * * ]
			// where the first factor on the right is formed from q.

			T nq = q.Magnitude();
			if (nq > (T)0) 
			{
				q /= nq;

				// replace D with D [q]

				// D(:, [m m+1]) = D(:, [m m+1]) * [q(1) -q(2); q(2) q(1)];
				GivensCol(D, m, m+1, q[0], q[1]);
				// V(:, [m m+1]) = V(:, [m m+1]) * [q(1) -q(2); q(2) q(1)];
				GivensCol(V, m, m+1, q[0], q[1]);

				// D now has this pattern:
				//      * * 0 0 0
				//      x * * 0 0
				//      0 0 * * 0
				//      0 0 0 * *
				//      0 0 0 0 *
				// so we just want to chase the x down the diagonal.

				VT_HR_EXIT(ChaseDown(U, D, V, m+1, m, m, n));
			}

			Truncate(D, m, n);

			VT_HR_END();
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::ChaseDown
		//
		//  Synopsis:  this function chases unwanted non-zero A(i, j) down the diagonal
		//
		//  Returns:   S_OK if all is well
		//
		//------------------------------------------------------------------------

		HRESULT ChaseDown(CMtx<T> &U, CMtx<T> &A, CMtx<T> &V, int i, int j, int /*m*/, int n) const
		{
			HRESULT hr = NOERROR;

			while((i<=n) && (j<=n))
			{
				if(i==(j+1))
				{
					//   * * 0 0 0 0
					//   0 * * 0 0 0
					// # 0 0 * * 0 0
					// # 0 0 x * * 0
					//   0 0 0 0 * *
					//   0 0 0 0 0 *

					// rotate A(i-1, j) and A(i, j).
					T c = A(i-1, j);
					T s = A(i, j);
					T r = VtHypot(c, s);

					if(r>(T)0) 
					{
						c/=r;
						s/=r;

						// A([i-1 i], :) = [c s; -s c] * A([i-1 i], :);
						GivensRow(A, i-1, i, c, -s);
						// U(:, [i-1 i]) = U(:, [i-1 i]) * [c -s; s c];
						GivensCol(U, i-1, i, c,  s);
					}

					i--;
					j+=2;
				}
				else if((i+1)==(j-1))
				{
					//       # #
					// * * 0 0 0 0 0
					// 0 * * 0 0 0 0
					// 0 0 * * x 0 0
					// 0 0 0 * * 0 0
					// 0 0 0 0 * * 0
					// 0 0 0 0 0 * *

					// rotate A(i, j-1) and A(i, j)
					T c = A(i, j-1);
					T s = A(i, j);
					T r = VtHypot(c, s);

					if(r>(T)0) 
					{
						c/=r;
						s/=r;

						// A(:, [j-1 j]) = A(:, [j-1 j]) * [c -s; s c];
						GivensCol(A, j-1, j, c, s);
						// V(:, [j-1 j]) = V(:, [j-1 j]) * [c -s; s c];
						GivensCol(V, j-1, j, c, s);
					}

					i+=2;
					j--;
				}
				else
					VT_HR_EXIT(E_FAIL);
			}

Exit:
			return hr;
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::GivensRow
		//
		//  Synopsis:  givens rotations
		//              pre-multiply rows i0, i1 by [c -s]
		//                                          [s  c]
		//
		//  Returns:   S_OK if all is well
		//
		//------------------------------------------------------------------------

		void GivensRow(CMtx<T> &X, int i0, int i1, T c, T s) const
		{
			int j;
			for (j=0; j<X.Cols(); ++j)
			{
				T tmp = c * X(i0, j) - s * X(i1, j);
				X(i1, j) = s * X(i0, j) + c * X(i1, j);
				X(i0, j) = tmp;
			}
		}

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveSVD<T>::GivensCol
		//
		//  Synopsis:  givens rotations
		//              post-multiply cols j0, j1 by [c -s]
		//                                           [s  c]
		//
		//  Returns:   S_OK if all is well
		//
		//------------------------------------------------------------------------

		void GivensCol(CMtx<T> &X, int j0, int j1, T c, T s) const
		{
			int i;
			for (i=0; i<X.Rows(); ++i) 
			{
				T tmp =  c * X(i, j0) + s * X(i, j1);
				X(i, j1) = -s * X(i, j0) + c * X(i, j1);
				X(i, j0) = tmp;
			}
		}

		CMtx<T> m_U, m_D, m_V;
	};

	typedef CSolveSVD<float> CSolveSVDf;
	typedef CSolveSVD<double> CSolveSVDd;

	// Eigenvalues and eigenvector of a 2x2 symmetric positive semi-definite matrix,
	//	A = U diag(e) U^T  or diag(e) = U^T A U
	//		(eigenvalues are all real and >= 0)
	template <class T>
	void VtSPSEigenValues(const CMtx2x2<T> &A, CMtx2x2<T> &U, CVec2<T>& e)
	{
		// Eigenvalues and eigenvector of a symmetric positive semi-definite matrix,
		//	A = U diag(e) U^T  or diag(e) = U^T A U
		//		(eigenvalues are all real and >= 0)

		// Form the characteristic equation, det|l I - A| = 0
		T b = A[0][0] + A[1][1];
		T c = A[0][0] * A[1][1]  - A[1][0] * A[0][1];
		VtNonNegativeRealQuadraticRoots(b, c, &e[0]);

		// Subtract the smaller eigenvalue from the diagonal
		T d0 = A[0][0] - e[1];
		T d1 = A[1][1] - e[1];

		// Find the first eigenvector from the larger row
		T row0 = abs(d0)+abs(A[0][1]);
		T row1 = abs(A[1][0])+abs(d1);
		CVec2<T> e0 = (row1 > row0) ? CVec2<T>(d1, -A[1][0]) :
			((row0 > 0) ? CVec2<T>(A[0][1], -d0) : CVec2<T>(1, 0));
		e0 = e0.Unit();

		// Use the eigenvector and its orthogonal complement
		U[0][1] =  e0[0], U[1][1] = e0[1];
		U[0][0] = -e0[1], U[1][0] = e0[0];
	}

	// SVD of 2x2 matrix
	//	A = U S V^T  or  S = U^T A V
	template <class T>
	inline void VtSVD2x2(const CMtx2x2<T>& A, CMtx2x2<T>& U, CVec2<T>& s, CMtx2x2<T>& VT)
	{
		// 2x2 singular value decomposition,
		//	A = U S V^T  or  S = U^T A V

		// Solve for the left singular vectors
		CVec2<T> e;
		VtSPSEigenValues(A * A.T(), U, e);
		s[0] = sqrt(e[0]);
		s[1] = sqrt(e[1]);

		// Solve for VT = S^-1 U^T A
		VT = U.T() * A;
		int i, j;
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
			{
				if (s[i] > 0)
					VT[i][j] /= s[i];
				else
					VT[i][j] = (i == 0) ? (j == 0) : // unit vector (0 eigenvalue)
					VT[0][1-j] * (j ? 1 : -1);
			}
	}

#if 0
	template <class Mtx, class Vec>
	double SVDTest(Mtx A, Vec s, int range = 2)
	{
		// Self-testing code:  try all matrices with integers in range [-range,range]
		int dim = A.Rows();
		int range2 = 2*range+1;
		int n_combinations = (int) pow(float(range2), float(dim*dim));
		int i, j, k;

		double max_dot = 0;		// maximum orthogonality error
		double max_err = 0;		// maximum reconstruction error
		for (i = 0; i < n_combinations; i++)
		{
			// Instantiate the matrix
			for (j = 0; j < dim; j++)
			{
				for (k = 0; k < dim; k++)
				{
					int l = j*dim + k;
					int pl = (int) pow(float(range2), float(l));
					int val = (i / pl) % range2 - range;
					A[j][k] = val;
				}
			}

			// Perform the SVD
			Mtx U, VT;
			VtSVD2x2(A, U, s, VT);

			// Check for orthogonality
			for (j = 0; j < dim; j++)
			{
				for (k = 0; k < dim; k++)
				{
					// Check U
					Vec& uj = *(Vec *) U[j];
					Vec& uk = *(Vec *) U[k];
					double udot = abs(uj * uk - (j == k));	// 0 or 1
					max_dot = __max(max_dot, udot); 

					// Check V
					Vec& vj = *(Vec *) VT[j];
					Vec& vk = *(Vec *) VT[k];
					double vdot = abs(vj * vk - (j == k));	// 0 or 1
					max_dot = __max(max_dot, vdot); 
				}
			}

			// Check the product
			Mtx S;  S = S.MakeDiag(s);
			Mtx B = U * S * VT;
			for (j = 0; j < dim; j++)
			{
				for (k = 0; k < dim; k++)
				{
					double err = abs(B[j][k] - A[j][k]);
					max_err = __max(max_err, err);
				}
			}
		}

		return __max(max_dot, max_err);
	}
#endif
}