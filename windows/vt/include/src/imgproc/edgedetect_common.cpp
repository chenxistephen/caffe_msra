//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      EdgeSegment class
//
//  History:
//      2011/11/10-kramnath
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_edgedetect_common.h"

#ifdef VT_GCC
#define _hypot hypot
#endif

using namespace vt;

void EdgeSegment::NormalizeNormalVector()
{
	length = (float)_hypot(n_x, n_y);
	if(length > 0)
	{
		n_x /= length;
		n_y /= length;
	}
	else
	{
		n_x = 1;
		n_y = 0;
	}
}

void EdgeSegment::ComputeThetaFromNormalVector()
{
	theta = (float)VT_180_PI * atan2(n_y, n_x);
	if (theta < -180)
		theta += 360;
	if (theta >= 180)
		theta -= 360;
}

void EdgeSegment::SetFromEndpoints(float x0, float y0, float x1, float y1)
{
	x = (x0 + x1) / 2;
	y = (y0 + y1) / 2;
	n_x = x1 - x0;
	n_y = y1 - y0;

	endpoint1.x = x0;
	endpoint1.y = y0;
	endpoint2.x = x1;
	endpoint2.y = y1;

	NormalizeNormalVector();

	ComputeThetaFromNormalVector();
}

void EdgeSegment::ComputeEndpoints() 
{
	float fltT = length / 2.0f;

	endpoint1.x = x - fltT * n_x;
	endpoint2.x = x + fltT * n_x;
	endpoint1.y = y - fltT * n_y;
	endpoint2.y = y + fltT * n_y;
}

// Transform a segment on the plane of projection and return the
// edgel corresponding to the projection of the transformed segment.
EdgeSegment operator*(const CMtx4x4d &transform, const EdgeSegment &edgel)
{
	float x0, y0, x1, y1;
	x0 = edgel.endpoint1.x;
	y0 = edgel.endpoint1.y;
    x1 = edgel.endpoint2.x;
	y1 = edgel.endpoint2.y;
	CVec4d v0(x0, y0, 1.0, 1.0);
	CVec4d v1(x1, y1, 1.0, 1.0);

	// Apply the transform
	CVec4d v0T = transform * v0;
	CVec4d v1T = transform * v1;

	// Project the results
	if ((v0T[2] != 0.0) && (v1T[2] != 0.0))
	{
		x0 = (float) (v0T[0] / v0T[2]);
		y0 = (float) (v0T[1] / v0T[2]);
		x1 = (float) (v1T[0] / v1T[2]);
		y1 = (float) (v1T[1] / v1T[2]);
	}
	else if (v0T[2] != 0.0)
	{
		x0 = (float) (v0T[0] / v0T[2]);
		y0 = (float) (v0T[1] / v0T[2]);
		x1 = x0;
		y1 = y0;
	}
	else if (v1T[2] != 0.0)
	{
		x1 = (float) (v1T[0] / v1T[2]);
		y1 = (float) (v1T[1] / v1T[2]);
		x0 = x1;
		y0 = y1;
	}
	else
	{
		// Not sure what to do in this case.
		x0 = y0 = x1 = y1 = 0.0;
	}

	EdgeSegment edgelRet;
	edgelRet.SetFromEndpoints(x0, y0, x1, y1);

	return edgelRet;
}


// Helper class
DisjointSets::DisjointSets()
{
}
DisjointSets::~DisjointSets()
{
}

HRESULT DisjointSets::Init(int n)
{
	VT_HR_BEGIN()
	VT_HR_EXIT(p.resize(n));
	VT_HR_EXIT(rank.resize(n));
	for(int i=0; i<n;i++)
	{
		rank[i] = 0;
		p[i]    = i;
	}
	VT_HR_END()
};

int DisjointSets::FindSet(int i)
{
	if (i!=p[i])
		p[i] = FindSet(p[i]);
	return p[i];
}

void DisjointSets::Union(int i, int j)
{
	Link(FindSet(i),FindSet(j));
}

void DisjointSets::Link(int i, int j)
{
	if (rank[i] > rank[j])
		p[j] = i; 
	else
	{
		p[i] = j;
		if (rank[i] == rank[j])
			rank[j]++;
	}
}

HRESULT DisjointSets::SetIds(vt::vector<int>& ids)
{
	VT_HR_BEGIN()
	VT_HR_EXIT(ids.resize(p.size()));
	for (int i=0;i<(int)p.size();i++)
		ids[i] = p[i];
	VT_HR_END()
}