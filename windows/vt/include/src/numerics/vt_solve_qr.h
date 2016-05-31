#pragma once

#include "vt_matrix.h"

namespace vt {

	// function [Q, R, p] = ms_qr(A)
	//
	// Economy size QR-decomposition with column pivoting.
	// http://mathworld.wolfram.com/QRDecomposition.html
	//
	// Input:
	//  A : m x n real matrix.
	//
	// Let k = min(m, n).
	//
	// Output:
	//  Q : m x k with orthonormal columns.
	//  R : k x n upper triangular.
	//  p : 1 x n permutation.
	//
	// where
	//  A(:, p) = Q R
	//
	// This is equivalent to [Q, R, p] = qr(A, 0).
	//

	template <class T>
	class CSolveQR : public CErrorBase
	{
	public:
		CSolveQR() {}
		CSolveQR(const CMtx<T> &mtx, bool permute = false) { Solve(mtx, permute); }

		HRESULT Solve(const CMtx<T> &mtx, bool permute = false)
		{
			VT_HR_BEGIN();

			CMtx<T> D, V;

			CVec<T> _v1, _v2, _v3;
			CVec<T> vtmp, vr, vd, v, vq;

			ClearError();

			VT_HR_EXIT(mtx.GetError());

			int i, j;
			int m = mtx.Rows();
			int n = mtx.Cols();
			int k = VtMin(m,n);

			VT_HR_EXIT(V.Create(m, k));
			V = (T)0;

			D = mtx;
			VT_HR_EXIT(D.GetError());

			// permutation vector
			VT_HR_EXIT(m_vP.resize(n));
			for(i=0; i<n; i++)
				m_vP[i] = i;

			VT_HR_EXIT(_v1.Create(m));
			VT_HR_EXIT(_v2.Create(m));
			VT_HR_EXIT(_v3.Create(m));

			// We are going to zap column j (= 3 below).
			//       #
			//   r r r r r r
			//   0 r r r r r
			// # 0 0 * * * *
			//   0 0 * * * *
			//   0 0 * * * *
			//   0 0 * * * *
			//   0 0 * * * *

			m_iDetQ = 1;
			m_iDetp = 1;

			for(j=0; j<k; j++)
			{
				// Permute if requested
				if( permute )
				{
					// find the column with the 
					//  largest norm
					int pivot = -1;
					T value = 0;
					for(i=j; i<n; i++)
					{
						//CVec<T> vtmp = D.GetCol(m_vP[i]).Extract(m-j, j);
						//VT_HR_EXIT(vtmp.GetError());

						vtmp.Wrap(_v1.Ptr(), m-j);
						D.GetColSlice(m_vP[i], j, m-j, vtmp);

						T thisvalue = vtmp.Magnitude();
						if(thisvalue > value)
						{
							pivot = i; 
							value = thisvalue;
						}
					}

					if((pivot != -1) && (pivot != j))
					{
						// swap it to position j
						int temp = m_vP[pivot];
						m_vP[pivot] = m_vP[j];
						m_vP[j] = temp;

						m_iDetp = -m_iDetp;
					}
				}

				// construct reflection vector to use.
				//CVec<T> vr = D.GetCol(m_vP[j]).Extract(m-j, j);
				//VT_HR_EXIT(vr.GetError());

				vr.Wrap(_v1.Ptr(), m-j);
				D.GetColSlice(m_vP[j], j, m-j, vr);
				VT_HR_EXIT(vr.GetError());

				if(vr[0] >= (T)0)
					vr[0] += vr.Magnitude();
				else
					vr[0] -= vr.Magnitude();

				if(vr.Size()>1) 
				{
					T norm_v = vr.Magnitude();
					if(norm_v > (T)0)
						vr /= norm_v;
					else
						vr[0] = (T)1;
				}
				else
					vr[0] = (T)1;

				//V.SetCol(j, V.GetCol(j).Update(vr,j));
				//VT_HR_EXIT(V.GetError());

				V.SetColSlice(j, j, m-j, vr);
				VT_HR_EXIT(V.GetError());

				// use reflection to eliminate stuff below position (j, j)
				for(i = j; i < n; i++)
				{
					//CVec<T> vd = D.GetCol(m_vP[i]).Extract(m-j, j);
					//VT_HR_EXIT(vd.GetError());

					vd.Wrap(_v2.Ptr(), m-j);
					D.GetColSlice(m_vP[i], j, m-j, vd);
					VT_HR_EXIT(vd.GetError());

					T fa = vr * vd;
					T fv = vr * vr;

					//vd -= vr * ((T)2.0 * fa/fv);
					//VT_HR_EXIT(vd.GetError());

					vtmp.Wrap(_v3.Ptr(), m-j);
					vtmp = vr;
					vtmp *= ((T)2.0 * fa/fv);
					vd -= vtmp;

					//D.SetCol(m_vP[i], D.GetCol(m_vP[i]).Update(vd, j));
					//VT_HR_EXIT(D.GetError());

					D.SetColSlice(m_vP[i], j, m-j, vd);
					VT_HR_EXIT(D.GetError());
				}

				m_iDetQ = -m_iDetQ;
			}

			// form Q.
			VT_HR_EXIT(m_Q.Create(m, k));
			m_Q.MakeI();

			// we have
			//    ... [v3] [v2] [v1] A = R
			// so
			//    Q = [v1] [v2] [v3] ...
			// and we can accumulate this in reverse without
			// knowing all m columns of Q.

			for(j = k-1; j >= 0; j--)
			{
				//CVec<T> v = V.GetCol(j).Extract(m-j, j);
				//VT_HR_EXIT(v.GetError());

				v.Wrap(_v1.Ptr(), m-j);
				V.GetColSlice(j, j, m-j, v);
				VT_HR_EXIT(v.GetError());

				T fvv = v * v;

				for(i=0; i < k; i++)
				{
					//CVec<T> vtmp = m_Q.GetCol(i).Extract(m-j, j);
					//VT_HR_EXIT(vtmp.GetError());

					vtmp.Wrap(_v2.Ptr(), m-j); 
					m_Q.GetColSlice(i, j, m-j, vtmp);
					VT_HR_EXIT(vtmp.GetError());

					T fa = v * vtmp;

					//CVec<T> vq = vtmp - v * (T)(((T)2.0) * fa/fvv);
					//VT_HR_EXIT(vq.GetError());

					vq.Wrap(_v3.Ptr(), m-j);
					vq = v;
					vq *= -(T)(((T)2.0) * fa/fvv);
					vq += vtmp;
					VT_HR_EXIT(vq.GetError());

					//m_Q.SetCol(i, m_Q.GetCol(i).Update(vq, j));
					//VT_HR_EXIT(m_Q.GetError());

					m_Q.SetColSlice(i, j, m-j, vq);
					VT_HR_EXIT(m_Q.GetError());
				}
			}

			// form R (just take upper triangle of D).
			VT_HR_EXIT(m_R.Create(k, n));
			m_R = (T)0;

			for(i=0; i < k; i++)
			{
				for(j=0; j < n; j++)
				{
					if(i > j)
						m_R(i,j) = (T)0;
					else
						m_R(i,j) = D(i, m_vP[j]);
				}
			}

			VT_HR_EXIT_LABEL();

			return SetError(hr);
		}

		void Free() { m_Q.Free(); m_R.Free(); m_vP.clear(); }

		const CMtx<T> &Q() const { return m_Q; }
		const CMtx<T> &R() const { return m_R; }
		const vector<int> &Perm() const { return m_vP; }

		//+-----------------------------------------------------------------------
		//
		//  Member:    CSolveQR<T>::Det
		//
		//  Synopsis:  Calculates the determinant from the QR decomp
		//
		//  Returns:   the determinant
		//
		//------------------------------------------------------------------------

		T Det() const
		{
			T detR = 1;
			int i;
			for (i=0; i<m_R.Rows() && i<m_R.Cols(); i++)
				detR *= m_R(i, i);
			return (T)m_iDetQ * detR * (T)m_iDetp;
		}

	private:
		CMtx<T> m_Q, m_R;
		vector<int> m_vP;
		int m_iDetQ, m_iDetp; // these are +/-1
	};

	typedef CSolveQR<float> CSolveQRf;
	typedef CSolveQR<double> CSolveQRd;
}