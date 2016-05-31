//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Numerical routines for solving matrix problems
//      These routines were developed by geoffrey cross MSR Cambridge and
//        adapted for vision tools by swinder.
//
//  History:
//      2003/11/16-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_complex.h"
#include "vt_matrix.h"
#include "vt_solve_svd.h"

namespace vt {

	//+-----------------------------------------------------------------------
	//
	//  Member:    VtBestNullSpaceVector
	//
	//  Synopsis:  returns the unit-norm vector v that minimizes ||Av||
	//
	//  Returns:   S_OK if all is well
	//
	//------------------------------------------------------------------------

	template <class T>
	inline HRESULT VtBestNullSpaceVector(const CMtx<T> &mtx, CVec<T> &vecrtn)
	{
		HRESULT hr = NOERROR;
		CSolveSVD<T> svd;
		if(mtx.Rows() < mtx.Cols())
			VT_HR_EXIT( E_INVALIDARG );
		VT_HR_EXIT( svd.Solve(mtx) );
		vecrtn = svd.GetBestNullSpaceVector();

Exit:
		return hr;
	}

	// solve a quadratic - returns the following values:
	// 0 - no roots e.g. when a and b are zero
	// 1 - one root e.g. when a=0
	// 2 - two real roots
	// 3 - two complex roots
	// roots are always returned largest first
	// coeffs should be provided in order a, b, c
	template <class T>
	int VtSolveQuadratic(const CVec3<T> &vCoeffs, CVec2<Complex<T> > &vRoots)
	{
		T a = vCoeffs[0];
		T b = vCoeffs[1];
		T c = vCoeffs[2];

		if(a==0 && b==0)
		{
			vRoots[0] = 0;
			vRoots[1] = 0;
			return 0; // and b are zero, and so there is no x and no roots
		}
		if(a==0)
		{
			vRoots[0] = Complex<T>(-c/b);
			vRoots[1] = 0;
			return 1;
		}
		T g = b*b - 4*a*c;
		if(g<0)
		{
			// complex roots
			T u = sqrt(-g)/(2*a);
			T v = -b/(2*a);
			vRoots[0] = Complex<T>(v, u);
			vRoots[1] = Complex<T>(v, -u);
			return 3;
		}

		// real roots
		T q;
		if(b>=0)
			q = -(T)0.5 * (b + sqrt(g));
		else
			q = -(T)0.5 * (b - sqrt(g));

		vRoots[0] = Complex<T>(q/a);
		vRoots[1] = Complex<T>(c/q);
		return 2;
	}


	// Compute the roots of the quadratic equation
	//	 x^2 - b x + c = 0
	//	Only non-negative real roots are allowed (for an symmetric matrix eigenvalues)
	template <class T>
	void VtNonNegativeRealQuadraticRoots(T b, T c, T x[2])
	{
		T avg = (T)0.5*b;
		T discr = avg*avg - c;
		discr = (discr > 0) ? sqrt(discr) : (T)0.0;
		x[0] = VtMin(avg + discr, b);	// larger eigenvalue
		x[1] = VtMax(avg - discr, T(0));	// smaller one
	}

	// solve a cubic ax^3 + bx^2 + cx + d = 0
	// if a=0 this calls the funtion above, return values are as follows:
	// 0 - no roots e.g. when a and b and c are zero
	// 1 - one root e.g. when a and b are zero
	// 2 - two real roots. when a is zero
	// 3 - two complex roots. when a is zero
	// 4 - three real roots
	// 5 - one real root and two complex roots (real root is first)
	// coeffs are in order a, b, c, d

	template <class T>
	int VtSolveCubic(const CVec4<T> &vCoeffs, CVec3<Complex<T> > &vRoots)
	{
		T s = vCoeffs[0];
		if(s==0)
		{
			CVec3<T> vC;
			vC[0] = vCoeffs[1];
			vC[1] = vCoeffs[2];
			vC[2] = vCoeffs[3];
			CVec2<Complex<T> > vRtn;
			int iStatus = VtSolveQuadratic(vC, vRtn);
			vRoots[0] = vRtn[0];
			vRoots[1] = vRtn[1];
			vRoots[2] = 0;
			return iStatus;
		}

		T a = vCoeffs[1];
		T b = vCoeffs[2];
		T c = vCoeffs[3];

		if(s!=1)
		{
			a/=s;
			b/=s;
			c/=s;
		}

		T aa = a*a;
		T q = (aa - 3*b)/9;
		T r = (2*a*aa - 9*a*b + 27*c)/54;
		T rr = r*r;
		T qqq = q*q*q;
		if(rr < qqq)
		{
			T th = acos(r/sqrt(qqq));
			T mag = -2*sqrt(q);
			T a3 = a/3;
			vRoots[0] = Complex<T>(mag*cos(th/3) - a3);
			vRoots[1] = Complex<T>(mag*cos((th + 2*VT_PI)/3) - a3);
			vRoots[2] = Complex<T>(mag*cos((th - 2*VT_PI)/3) - a3);
			return 4;
		}

		T w = pow(fabs(r) + sqrt(rr - qqq), 1/(T)3.0);
		if(w * r > 0)
			w *= -1;
		T f = 0;
		if(w!=0)
			f = q/w;
		T u = -(T)0.5*(w + f) - a/3;
		T v = (T)(VT_SQRT3_2) * (w - f);
		vRoots[0] = Complex<T>(w + f - a/(T)3.0);
		vRoots[1] = Complex<T>(u, v);
		vRoots[2] = Complex<T>(u, -v);
		return 5;
	}


	// Compute the roots of the cubic equation
	//	x^3 + a x^2 + b x + c = 0
	//	Only non-negative real roots are allowed (for an symmetric matrix eigenvalues)
	template <class T>
	void VtNonNegativeRealCubicRoots(T a, T b, T c,  T x[3])
	{
		// Compute the roots of the cubic equation
		//	x^3 + a x^2 + b x + c = 0
		//	Only non-negative real roots are allowed (for an symmetric matrix eigenvalues)
		x[0] = x[1] = x[2] = -a/3;		// default return
		T a2 = a*a;
		T a3 = a2*a;
		T q = (a2 - 3*b)/9;
		T r = (2*a3 - 9*a*b + 27*c)/54;
		if (q <= 0)
			return;
		T qr = sqrt(q);
		T ratio = VtMin(T(1.0), T(r / pow(qr, 3)));
		T t = acos(ratio);
		int i;
		for (i = 0; i < 3; i++)
			x[i] += 2*qr*cos((t + 2*i*(T)VT_PI)/3);
	}

};

