//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Quaternion class
//
//  History:
//      2003/11/12-swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_matrix2x2.h"
#include "vt_matrix3x3.h"
#include "vt_matrix4x4.h"

namespace vt {

	template <class T>
	class CQuaternion
	{
	public:
		CQuaternion() {}
		CQuaternion(const CVec3<T> &vec, const T &f) : v(vec), w(f) {}

		CVec4<T> GetVec4() const { return CVec4<T>(v.x, v.y, v.z, w); }
		CQuaternion<T> & MakeFromVec4(const CVec4<T> & v4) { v[0] = v4.x; v[1] = v4.y; v[2] = v4.z; w = v4.w; return *this; }

		T MagnitudeSq() const {return v.MagnitudeSq() + w * w; }
		T Magnitude() const { return sqrt(MagnitudeSq()); }
		CQuaternion<T> Unit() const { return (*this)/Magnitude(); }
		CQuaternion<T> Conj() const { return CQuaternion<T>(-v, w); }
		CQuaternion<T> Inv() const { return Conj()/MagnitudeSq(); }
		CVec3<T> RotateVec(const CVec3<T> &vec) const { return ((*this) * CQuaternion<T>(vec, 0) * Inv()).v; }
		CQuaternion<T> & MakeI() { v = CVec3<T>(0, 0, 0); w = T(1); return *this; }

		CQuaternion<T> operator*(const CQuaternion<T> &q) const
		{
			//  Cross(v, q.v) + w*q.v + q.w*v, w*q.w - Dot(v, q.v)
#if 0   // this formula does not return the exact value for *this or q = (0,0,0,1)
			T A = (w + v.x) * (q.w + q.v.x);
			T B = (v.z - v.y) * (q.v.y - q.v.z);
			T C = (v.x - w   ) * (q.v.y + q.v.z);
			T D = (v.y + v.z) * (q.v.x - q.w);
			T E = (v.x + v.z) * (q.v.x + q.v.y);
			T F = (v.x - v.z) * (q.v.x - q.v.y);
			T G = (w + v.y) * (q.w - q.v.z);
			T H = (w - v.y) * (q.w + q.v.z);
			CVec3<T> new_v(  A - 0.5f*( E + F + G + H),
				-C + 0.5f*( E - F + G - H),
				-D + 0.5f*( E - F - G + H));
			T new_w =        B + 0.5f*(-E - F + G + H);
#else   // added by Rick, does the right thing for (0,0,0,1) operands
			CVec3<T> new_v = v.Cross(q.v) + q.v*w + v*q.w;
			T        new_w = w*q.w - v * q.v;
#endif
			CQuaternion<T> new_q(new_v, new_w);
			return new_q;
		}

		CQuaternion<T>& 
			MakeFromEulerAngles(const T& pitch, const T& yaw, const T& roll,
			bool bSwapXY = false )
		{
			// Convert roll (z), pitch (x), yaw(y) into a 3D rotation
			// The particular order is chosen as the natural one to go with
			//  3D rotations rendered on a spherical surface
			// 22-Apr-04: for "vertical" panoramas, swing up by (-x) before right by y
			CQuaternion<T> qx, qy, qz;
			qx.MakeFromAxisAngle(CVec3<T>(1, 0, 0), pitch);
			qy.MakeFromAxisAngle(CVec3<T>(0, 1, 0), yaw);
			qz.MakeFromAxisAngle(CVec3<T>(0, 0, 1), roll);
			if (bSwapXY)
				*this = qx * (qy * qz);
			else
				*this = qy * (qx * qz);
			return *this;
		}

		CQuaternion<T> & MakeFromAxisAngle(const CVec3<T> &vAxis, 
			const T &fAngle)
		{
			T fSin, fCos;
			VtSinCos( T(0.5)*fAngle, &fSin, &fCos );

			if( vAxis.MagnitudeSq() == 0 )
			{
				v = vAxis;
			}
			else
			{
				v = vAxis.Unit() * fSin;
			}
			w = fCos;

			return *this;
		}

		void GetAxisAngle(CVec3<T> &vAxis, T &fAngle) const
		{
			CQuaternion<T> qu = Unit();
			fAngle = acos(qu.w);
			if( abs(fAngle) < 0.00001 )
			{
				vAxis = CVec3<T>(0, 0, 0); 
			}
			else
			{
				vAxis = qu.v / sin(fAngle);
				vAxis = vAxis.Unit();
			}
			fAngle *= T(2);
		}

		CMtx3x3<T> GetRotationMatrix() const
		{
			// This code must have originally been borrowed from D3D because
			//  it returned a rotation matrix for row vectors.  It's been changed
			//  (Rick, 29-Dec-03) to return correct column vector R, i.e., it's now
			//  consistent with Faugeras' book and the quaternion multiply code.
			CMtx3x3<T> m;

			T x2 = v.x + v.x;
			T y2 = v.y + v.y;
			T z2 = v.z + v.z;

			T wx2 = w * x2;
			T wy2 = w * y2;
			T wz2 = w * z2;
			T xx2 = v.x * x2;
			T xy2 = v.x * y2;
			T xz2 = v.x * z2;
			T yy2 = v.y * y2;
			T yz2 = v.y * z2;
			T zz2 = v.z * z2;

			m(0,0) = T(1) - yy2 - zz2;
			m(0,1) = xy2 - wz2;
			m(0,2) = xz2 + wy2;

			m(1,0) = xy2 + wz2;
			m(1,1) = T(1) - xx2 - zz2;
			m(1,2) = yz2 - wx2;

			m(2,0) = xz2 - wy2;
			m(2,1) = yz2 + wx2;
			m(2,2) = T(1) - xx2 - yy2;

			return m;
		}

		CMtx4x4<T> GetRotationMatrix4x4() const
		{
			CMtx3x3<T> m = GetRotationMatrix();
			return CMtx4x4<T>(m(0,0), m(0,1), m(0,2), 0,
				m(1,0), m(1,1), m(1,2), 0,
				m(2,0), m(2,1), m(2,2), 0,
				0,      0,      0,      1);
		}

		CQuaternion<T> & MakeFromInterpolation(const CQuaternion<T> &qS, 
			const CQuaternion<T> &qE,
			T f, bool *pErr = NULL)
		{
			T omega, cosOmega, sinOmega;

			CQuaternion<T> q1 = qS.Unit();
			CQuaternion<T> q2 = qE.Unit();

			//  choose the closest great-circle route
			if ((q1 - q2).MagnitudeSq() > (q1 + q2).MagnitudeSq())
				q2 = -q2;

			cosOmega = VtMin(T(1), q1.v * q2.v + q1.w * q2.w);
			sinOmega = sqrt(T(1) - cosOmega*cosOmega);

			//  if omega is nearly zero, just lerp
			//  if its nearly PI, then Slerp is ill defined

			if(pErr)
				*pErr = false;

			if ((cosOmega < 0) && (sinOmega < T(0.00001)))
			{
				if(pErr)
					*pErr = true;
				v = (T)0;
				w = (T)0;
				return *this;
			}

			if ((cosOmega > 0) && (sinOmega < T(0.00001)))
			{
				*this = q1*(T(1) - f) + q2*f;
				return *this;
			}

			omega = acos(cosOmega);
			*this = q1*sin((T(1) - f)*omega) + q2*sin(f*omega);
			*this /= sinOmega;

			return *this;
		}

		CQuaternion<T>& 
			MakeFromVectorRotation(const CVec3<T>& vPosOrig,
			const CVec3<T>& vPosDest)
		{
			// compute quaternion that rotates vPosOrig to vPosDest
			CVec3<T> vCross = vPosOrig.Cross(vPosDest);
			T fMag = vCross.Magnitude();
			if( fMag < T(0.00001) )
			{
				v = T(0);
				w = T(1);
			}
			else
			{
				T fCos2 = vPosOrig * vPosDest;
				T fSin  = sqrt((T(1)-fCos2)/T(2));
				T fCos  = sqrt((T(1)+fCos2)/T(2));

				v = vCross * fSin / fMag;
				w = fCos;
			}
			return *this;
		}


		CQuaternion<T> & MakeFromMatrix(const CMtx3x3<T>& M)
		{
			// Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
			// article "Quaternion Calculus and Fast Animation".  (Taken from GDMAG feb'98 p38)
			// Again, fixed by Rick (29-Dec-03) to be consistent with column vectors.
			T trace = M[0][0] + M[1][1] + M[2][2];
			T root;

			if ( trace > 0. )
			{
				// |w| > 1/2, may as well choose w > 1/2

				root = sqrt(trace + T(1));  // 2w
				w = T(0.5) * root;

				root = T(0.5) / root;  // 1/(4w)
				v.x = -(M[1][2] - M[2][1]) * root;
				v.y = -(M[2][0] - M[0][2]) * root;
				v.z = -(M[0][1] - M[1][0]) * root;
			}
			else
			{
				// |w| <= 1/2
				static const int next[3] = { 1, 2, 0 };

				int i;
				i = M[1][1] > M[0][0] ? 1 : 0;
				i = M[2][2] > M[i][i] ? 2 : i;

				int j = next[i];
				int k = next[j];

				root = sqrt(M[i][i] - M[j][j] - M[k][k] + T(1));
				v[i] = T(0.5) * root;

				if(0 != root)
					root = T(0.5) / root;

				w        = -(M[j][k] - M[k][j]) * root;
				v[j] = (M[i][j] + M[j][i]) * root;
				v[k] = (M[i][k] + M[k][i]) * root;
			}

			return *this;
		}

		CQuaternion<T> & MakeFromMatrix4x4(const CMtx4x4<T>& M)
		{
			CMtx3x3<T>  M3x3(M(0,0), M(0,1), M(0,2),
				M(1,0), M(1,1), M(1,2),
				M(2,0), M(2,1), M(2,2));
			return MakeFromMatrix(M3x3);
		}

		CVec3<T> GetEulerAngles(bool bSwapXY = false) const
		{
			// Compute the Euler angles corresponding to a "usual"
			//  mapping from a 3D rotation to a spherical map.
			// The formula for y = yaw(about y), r = roll(about z), p=pitch(about x)
			// sy = sin(yaw), cy = cos(yaw), etc
			//  [ cy*cr + sy*sr*sp  |-cy*sr + sy*cr*sp  | sy*cp ] 
			//  [ sr*cp             | cr*cp             | -sp   ]
			//  [ -sy*cr + cy*sr*sp | sy*sr+cy*cr*sp    | cy*cp ]
			// 22-Apr-04: for "bSwapXY" panoramas, we use this instead:

			CMtx3x3<T> R = GetRotationMatrix();
			CVec3<T> a;
			if (!bSwapXY)
			{
				double sp = R[1][2];
				a.x = (T) asin(-sp);
				double sy_cp = R[0][2];
				double cy_cp = R[2][2];
				double sr_cp = R[1][0];
				double cr_cp = R[1][1];
				if ((sy_cp != 0 || cy_cp != 0) &&
					(sr_cp != 0 || cr_cp != 0))
				{
					// Usual case, not looking at N/S pole
					a.y = (T) atan2(sy_cp, cy_cp);
					a.z = (T) atan2(sr_cp, cr_cp);
				}
				else
				{
					// Special case, looking at pole, set a.z = 0
					a.z = 0;
					double sy = -R[2][0];
					double cy =  R[0][0];
					a.y = (T) atan2(sy, cy);
				}
			}
			else
			{
				double sy = R[0][2];
				a.y = (T) asin(sy);
				double sx_cy = R[1][2];
				double cx_cy = R[2][2];
				double sz_cy =-R[0][1];
				double cz_cy = R[0][0];
				if ((sx_cy != 0 || cx_cy != 0) &&
					(sz_cy != 0 || cz_cy != 0))
				{
					// Usual case, not looking at pole
					a.x = (T) atan2(sx_cy, cx_cy);
					a.z = (T) atan2(sz_cy, cz_cy);
				}
				else
				{
					// Special case, looking at pole, set a.z = 0
					a.z = 0;
					double sx = -R[2][1];
					double cx =  R[1][1];
					a.x = (T) atan2(sx, cx);
				}
			}
			return a;
		}

		/*
		template <class T>
		void
		SphericalEulerSelfTest()
		{
		// Generate some random quaternions, test conversions
		static int n_iter = 100;
		CRand rand;
		double max_err = 0.0;
		for (int i = 0; i < n_iter; i++)
		{
		// Generate a random quaternion
		CQuaternion<T> q;
		for (int j = 0; j < 3; j++)
		q.v[j] = (T) rand.Gauss();
		q.w = (T) rand.Gauss();
		q = q.Unit();
		bool bSwapXY = (i & 1) != 0;

		// Map to angles and back again
		CVec3<T> a = q.GetSphericalEuler(bSwapXY);
		CQuaternion<T> q2;
		q2.MakeFromSphericalEuler(a, bSwapXY);

		// Measure the error
		CQuaternion<T> qd = q2 * q.Inv();
		double err = qd.v.Magnitude();
		max_err = (i == 0) ? err : __max(err, max_err);
		}
		VT_ASSERT(max_err < 1.0e-6);   // reasonable test, good breakpoint
		}
		*/

		T Dot(const CQuaternion &q) const { return v * q.v + w * q.w; }

		//NOTE: solid-body rotations are combined with *, not +
		CQuaternion<T> operator +() const { return *this; }
		CQuaternion<T> operator -() const { return CQuaternion<T>(-v, -w); }
		CQuaternion<T> operator +(const CQuaternion<T> &q) const { return CQuaternion<T>(v+q.v, w+q.w); }
		CQuaternion<T> operator -(const CQuaternion<T> &q) const { return CQuaternion<T>(v-q.v, w-q.w); }
		CQuaternion<T> operator *(const T &f) const { return CQuaternion<T>(v * f, w * f); }
		//    friend CQuaternion<T> operator *(const T &f, const CQuaternio<T> &q);
		CQuaternion<T> operator /(const T &f) const { return CQuaternion<T>(v / f, w / f); }

		CQuaternion<T> & operator *=(const CQuaternion<T> &q) { return (*this = *this * q); }
		CQuaternion<T> & operator +=(const CQuaternion<T> &q) { v += q.v; w += q.w; return *this; }
		CQuaternion<T> & operator -=(const CQuaternion<T> &q) { v -= q.v; w -= q.w; return *this; }
		CQuaternion<T> & operator *=(const T &f) { v *= f; w *= f; return *this; }
		CQuaternion<T> & operator /=(const T &f) { v /= f; w /= f; return *this; }

		bool operator ==(const CQuaternion<T> &q) const { return v==q.v && w==q.w; }
		bool operator !=(const CQuaternion<T> &q) const { return !operator==(q); }

	public:
		CVec3<T> v;
		T w;
	};

	// friend function
	template <class T>
	inline CQuaternion<T> operator*(const T &f, const CQuaternion<T> &q)
	{
		return CQuaternion<T>(q.v * f, q.w * f);
	}
	typedef CQuaternion<float> CQuaternionf;
	typedef CQuaternion<double> CQuaterniond;

};

