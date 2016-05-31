//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Edge Segment class
//
//  History:
//      2011/11/10-kramnath
//          Created
//
//------------------------------------------------------------------------

#include "vtcommon.h"

namespace vt {

#pragma once

    /// \ingroup edgedetect
	/// <summary> Representation of an edge segment  </summary>
	struct EdgeSegment 
	{
		/// <summary> sub-pixel edge position </summary>
		float x, y;         
		/// <summary> orientation, as normal vector </summary>
		float n_x, n_y;     
		/// <summary>  orientation, as angle (degrees) </summary>
		float theta;        
		/// <summary> length of edgel within 1x1 rect </summary>
		float length;   
		/// <summary> strength of edgel (local gradient magnitude SQUARED) </summary>
		float strength;     
		/// <summary> in case of steerable filters this is the phase angle (degrees) </summary>
		float phi;			
		/// <summary> length of line (estimated from ellipsoid) </summary>
		float lineLength;  
		/// <summary> estimated std. dev. of edgel noise </summary>
		float sigma;        
		/// <summary> End point 1 of the edgel </summary>
		CVec2f endpoint1;
		/// <summary> End point 2 of the edgel </summary>
		CVec2f endpoint2;
		/// <summary> curvature of the line joining the current point and previous
		/// point on the curve </summary>
		float curvature;

		/// <summary> normalizes the normal vector (n_x, n_y) based on the length of the edgesegment </summary>
		void NormalizeNormalVector();

		/// <summary> gets the angle in degrees of the orientation of the edgesegment from normal vector (n_x, n_y) </summary>
		void ComputeThetaFromNormalVector();

		/// <summary> Sets x, y, n_x, n_y, endpoint1, endpoint2 and theta from the given endpoints. n_x and n_y will be normalized </summary>
		void SetFromEndpoints(float x0, float y0, float x1, float y1);

		/// <summary> Computes the end points using x, y n_x, n_y and length and sets them in ep1 and ep2 </summary>
		void ComputeEndpoints();
	};

    /// \ingroup edgedetect
	/// <summary> Transform a segment on the plane of projection and return the
	/// edgel corresponding to the projection of the transformed segment.</summary>
	EdgeSegment operator*(const CMtx4x4d &transform, const EdgeSegment &edgel);

	/// <summary> A line segment (used by the line detector)</summary>
	struct LineSegment {
		/// <summary> start point of the line segment</summary>
		CVec2d  start;
		/// <summary> endpoint of the line segment </summary>
		CVec2d  end;
		/// <summary> geometric length </summary>
		float length;      
	};

    /// \ingroup edgdetect
	/// <summary> Steerable filters can output two types of edge responses
	/// either step or bar. The enum allows to select whether to output only step
	/// or bar or both </summary>
	enum eSteerableEdgeType
	{
		/// <summary> output all edges</summary>
		eSteerableEdgeTypeAllEdges=0,
		/// <summary> output evenedges</summary>
		eSteerableEdgeTypeBarEdges=1,
		/// <summary> odd edges</summary>
		eSteerableEdgeTypesStepEdges=2
	};

    /// \ingroup edgedetect
	/// <summary> A vanishing point </summary>
	struct VanishingPoint 
	{
		/// <summary> score</summary>
		float score;
		/// <summary> 2D vanishing point in homogeneous coordinates</summary>
		CVec3d vpt;
		/// <summary> </summary>
		vector<int> idx;
	};


	// Helper class 
	struct DisjointSets
	{
		DisjointSets();
		~DisjointSets();
		HRESULT Init(int n);
		int     FindSet(int i);
		void    Union(int i, int j);
		void    Link(int i, int j);
		HRESULT SetIds(vector<int>& ids);
		vector<int> p;
		vector<int> rank;
	};

}