//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      functions to carry out Nelder-Mead method simplex minimization
// http://en.wikipedia.org/wiki/Nelder-Mead_method
//
//  History:
//      2006/10/22-swinder
//			Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

namespace vt {

#define SIMPLEX_TOL             0.00001f
#define SIMPLEX_MAX_MOVES       10000

	// use the simplex method to minimize the given function of N parameters
	// the number of parameters must be greater than one. this method does not require
	// the computation of derivatives. if you can easily compute derivatives, then
	// use a faster minimization technique, e.g. the marquardt class.
	// vStartEnd contains the vector of starting parameters and after the minimization
	// it gives the best ending parameters.
	// vDelta must be initialized to a reasonable stepsize along each of the dimensions
	// of the parameters. max and tol determine the ending conditions. the maximum number
	// of steps is given by max and the ratio of value change is given by tol

	template <class T>
	HRESULT SimplexFindBetter(int iHi, T fDelta, 
		HRESULT (*pFunc)(const vector<T> &vParams, T &fRtn, void *), 
		void *pUser, vector<vector<T> > &vvSet, vector<T> &vValue,
		vector<T> &vSum, T &fTry)
	{
		HRESULT hr = NOERROR;
		vector<T> vTry;

		int size = (int)vvSet[0].size();

		T fD1 = (1-fDelta)/size;
		T fD2 = fD1 - fDelta;
		int j;

		VT_HR_EXIT( vTry.resize(size) );

		for(j=0; j<size; j++)
			vTry[j] = vSum[j] * fD1 - vvSet[iHi][j] * fD2;

		VT_HR_EXIT( (*pFunc)(vTry, fTry, pUser) );

		if(fTry < vValue[iHi])
		{
			vValue[iHi] = fTry;
			for(j=0; j<size; j++)
			{
				vSum[j] += vTry[j] - vvSet[iHi][j];
				vvSet[iHi][j] = vTry[j];
			}
		}

Exit:
		return hr;
	}

	template <class T>
	HRESULT VtSimplexSearch(vector<T> &vStartEnd, const vector<T> &vDelta, 
		HRESULT (*pFunc)(const vector<T> &vParams, T &fValue, void *pUserData),
		void *pUser, int iMax, T fTol)
	{
		VT_HR_BEGIN();

		vector<vector<T> > vvSet;
		vector<T> vValue;
		vector<T> vSum;

		if(pFunc==NULL || vStartEnd.size()<1 || vStartEnd.size() != vDelta.size())
			VT_HR_EXIT( E_INVALIDARG );

		int size = (int)vStartEnd.size();
		int points = size + 1;

		VT_HR_EXIT( vvSet.resize(points) );
		VT_HR_EXIT( vValue.resize(points) );
		int i, j;
		for(i=0; i<points; i++)
		{
			VT_HR_EXIT( vvSet[i].resize(size) );
			for(j=0; j<size; j++)
				vvSet[i][j] = vStartEnd[j] + ((i==j) ? vDelta[j] : 0);
			VT_HR_EXIT( (*pFunc)(vvSet[i], vValue[i], pUser) );
		}

		VT_HR_EXIT( vSum.resize(size) );
		memset(&(vSum[0]), 0, size * sizeof(T));
		for(i = 0; i<points; i++)
			for(j = 0; j<size; j++)
				vSum[j] += vvSet[i][j];

		int iCount = 0;
		for (;;)
		{
			int iHi = 1;
			int iNHi = 0;
			int iLo = 0;
			if(vValue[0] > vValue[1])
			{
				iHi = 0; iNHi = 1;
			}
			for(i=0; i<points; i++)
			{
				if(vValue[i] <= vValue[iLo])
					iLo = i;
				if(vValue[i] > vValue[iHi])
				{
					iNHi = iHi;
					iHi = i;
				}
				else if(vValue[i] > vValue[iNHi] && i != iHi)
					iNHi = i;
			}

			T fRatio = 2 * fabs(vValue[iHi]- vValue[iLo]) / (fabs(vValue[iHi]) + fabs(vValue[iLo]));
			if(fRatio < fTol)
			{
				for(j=0; j<size; j++)
					vStartEnd[j] = vvSet[iLo][j];
				return NOERROR;
			}

			if(iCount > iMax)
			{
				for(j=0; j<size; j++)
					vStartEnd[j] = vvSet[iLo][j];
				VT_HR_EXIT( E_NOTFOUND );
			}

			iCount++;

			T fTry;
			VT_HR_EXIT( SimplexFindBetter(iHi, (T)-1.0, pFunc, pUser, vvSet, vValue, vSum, fTry) );
			if(fTry < vValue[iLo])
			{
				VT_HR_EXIT( SimplexFindBetter(iHi, (T)2.0, pFunc, pUser, vvSet, vValue, vSum, fTry) );
			}
			else if(fTry >= vValue[iNHi])
			{
				T fSave = vValue[iHi];
				VT_HR_EXIT( SimplexFindBetter(iHi, (T)0.5, pFunc, pUser, vvSet, vValue, vSum, fTry) );

				if(fTry >= fSave)
				{
					for(i=0; i<points; i++)
						if(i!=iLo)
						{
							for(j=0; j<size; j++)
								vvSet[i][j] = vSum[j] = (T)0.5 * (vvSet[i][j] + vvSet[iLo][j]);
							VT_HR_EXIT( (*pFunc)(vSum, vValue[i], pUser) );
						}

						memset(&(vSum[0]), 0, size * sizeof(T));
						for(i = 0; i<points; i++)
							for(j = 0; j<size; j++)
								vSum[j] += vvSet[i][j];
				}
			}
		}

		VT_HR_END();
	}

};

