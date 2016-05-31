//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routine for extracting curves and lines
//
//  History:
//      2011/11/10-kramnath
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"
#include "vt_linedetector.h"

using namespace vt;

// local function declarations **** //

// Find connected components to find curves
HRESULT FindConnectedComponents(OUT vector< vector<int> >& curves, const vector<EdgeSegment>& edgelList, 
							 CIntImg &edgeIndexImg, const LineDetectorParams &ldParams, bool usedZeroCrossings);

// Get valid neighbors for a given edge segment
HRESULT GetNeighbors(OUT vector<CPoint>& nbrs, const vector<EdgeSegment> &edgelList, const CIntImg &edgeIndexImg, 
				CByteImg &status, int col, int row, const LineDetectorParams &ldParams, bool usedZeroCrossings);

// line simplification algorithm
HRESULT LineSimplify(vector<LineSegment>& lineCollection, const vector<int>& curveSeg, const vector<EdgeSegment>& edgelList, 
				const LineDetectorParams& ldParams);

// Merge segments (CURRENTLY UNUSED)
HRESULT MergeLineSegments(OUT vector<LineSegment>& mergedLines, const vector<LineSegment>& lines, int& nmerged, 
					const LineDetectorParams& ldParams);

// Support functions for merge segments ***
HRESULT ComputeDSets(vector<LineSegment>& data, vector<vector<int> >& sets, float distThreshold, float angleThreshold);
HRESULT LineSegmentsAreCollinearAndClose(bool& value, float distThreshold, float angleThreshold, LineSegment& a, LineSegment& b);
HRESULT ProjectOnLine(CVec2d& proj, CVecd& leq, CVec2d& p);

// Compute the curvature of the line joining two points. Curvature at a point is the curvature of the line 
// joining the current point and the previous point in the curve. The curvature changes sign on either side of the 
/// inflection point of the curve. This module records the curvature in the edgeSegment structure 
void ComputeCurvatureofCurves(vector< vector<int> >& curves, vector<EdgeSegment>& edgelList);

// *** support functions for merge segments END //

// **** local function declarations END //

HRESULT
	vt::VtFindCurveSegmentsFromEdges(OUT vector< vector<int> >& curves, vector<EdgeSegment>& edgelList, bool 
								roundToNearestQuadOrPixel, int width, int height, const LineDetectorParams& ldParams)
{
	VT_HR_BEGIN()

	// Initialize edgel image
	CIntImg edgeIndexImg;
	VT_HR_EXIT( edgeIndexImg.Create(width, height) );
	VT_HR_EXIT( edgeIndexImg.Fill((int)-1) );
	int x,y;
	for(int i = 0; i < (int) edgelList.size(); ++i)
	{
		const EdgeSegment& edgel = edgelList[i];
		if(roundToNearestQuadOrPixel)
		{	// if detected using zero crossings, round off to the nearest quad
			x = (int) (edgel.x);
			y = (int) (edgel.y);
		}
		else
		{	// if not detected using zero crossings, round off to the nearest pixel
			x = (int) (edgel.x+0.5f);
			y = (int) (edgel.y+0.5f);
		}
		if(x > 0 && x < width-1 && y > 0 && y < height-1)
			edgeIndexImg.Pix(x,y) = i;
	}

	// reserve
	VT_HR_EXIT( curves.resize(0) );
	VT_HR_EXIT( curves.reserve(width*height/4) );
	VT_HR_EXIT( FindConnectedComponents(curves, edgelList, edgeIndexImg, ldParams, roundToNearestQuadOrPixel) );

	// Compute curvature of curves and write it to the edgesegment list
	ComputeCurvatureofCurves(curves, edgelList);

	VT_HR_END()

}

HRESULT
	FindConnectedComponents(vector< vector<int> >& curves, const vector<EdgeSegment>& edgelList, CIntImg &edgeIndexImg, 
						const LineDetectorParams& ldParams, bool usedZeroCrossings)
{
	VT_HR_BEGIN()

	int width = edgeIndexImg.Width();
	int height = edgeIndexImg.Height();

	// status
	CLumaByteImg status;
	VT_HR_EXIT(status.Create(width, height));
	VT_HR_EXIT(status.Fill((Byte) 0));

	// Scan the image for edge segments
	for(int r = 1; r < height-1; ++r)
	{
		int *edgePix = edgeIndexImg.Ptr(0, r);
		Byte *statusPix = status.Ptr(0, r);
		for(int c = 1; c < width-1; ++c)
		{	
			// if not edgel or already visited, continue
			if(edgePix[c] < 0 || statusPix[c] > 0)
			{
				continue;
			}

			vector<int> curve;
			// Get neighbours for the root node
			vector<CPoint> nbrs;
			VT_HR_EXIT( GetNeighbors(nbrs, edgelList, edgeIndexImg, status, c, r, ldParams, usedZeroCrossings) );

			// indices
			int currCol = 0; 
			int currRow = 0;

			// Only one neighbor, grow one list
			if(nbrs.size() == 1)
			{
				// Push the root node
				VT_HR_EXIT( curve.push_back(edgeIndexImg.Pix(c, r)) );
				// set the status
				status.Pix(c,r) = 1;
				// Push the neighbor
				currCol = nbrs[0].x;
				currRow = nbrs[0].y;
				VT_HR_EXIT( curve.push_back(edgeIndexImg.Pix(currCol, currRow)) );
				// get neighbors and keep adding to the list
				for(;;)
				{
					status.Pix(currCol, currRow) = 1;
					vector<CPoint> chainNbrs;
					VT_HR_EXIT( GetNeighbors(chainNbrs, edgelList, edgeIndexImg, status, currCol, currRow, ldParams, usedZeroCrossings) );
					if(chainNbrs.size() >= 1)
					{
						// Greedily add the first neighbor that passes neighborhood test
						currCol = chainNbrs[0].x;
						currRow = chainNbrs[0].y;
						VT_HR_EXIT( curve.push_back(edgeIndexImg.Pix(currCol, currRow)) );
					}
					else
						break;
				}
			}
			else if(nbrs.size() >= 2) // two neighbors, grow on both sides
			{
				vector<int> l1, l2;
				// Push the root node
				VT_HR_EXIT( l1.push_back(edgeIndexImg.Pix(c, r)) );
				VT_HR_EXIT( l2.push_back(edgeIndexImg.Pix(c, r)) );
				// Mark the status of the root node
				status.Pix(c,r) = 1;
				// Add the first neighbor
				currCol = nbrs[0].x;
				currRow = nbrs[0].y;
				VT_HR_EXIT( l1.push_back(edgeIndexImg.Pix(currCol, currRow)) );

				// grow list 1
				for(;;)
				{
					// Set the status of the added pixel
					status.Pix(currCol, currRow) = 1;
					vector<CPoint> chainNbrs;
					VT_HR_EXIT( GetNeighbors(chainNbrs, edgelList, edgeIndexImg, status, currCol, currRow, ldParams, usedZeroCrossings) );
					if(chainNbrs.size() >= 1 )
					{	// Greedily add the first neighbor that passes neighborhood test					
						currCol = chainNbrs[0].x;
						currRow = chainNbrs[0].y;
						VT_HR_EXIT( l1.push_back(edgeIndexImg.Pix(currCol, currRow)) );
					}
					else
						break;
				}

				// Add the second neighbor
				currCol = nbrs[1].x;
				currRow = nbrs[1].y;
				VT_HR_EXIT( l2.push_back(edgeIndexImg.Pix(currCol, currRow)) );

				// grow list 2
				for(;;)
				{
					// Set the status of the added pixel
					status.Pix(currCol, currRow) = 1;
					vector<CPoint> chainNbrs;
					VT_HR_EXIT( GetNeighbors(chainNbrs, edgelList, edgeIndexImg, status, currCol, currRow, ldParams, usedZeroCrossings) );
					if(chainNbrs.size() >= 1 )
					{	// Greedily add the first neighbor that passes neighborhood test					
						currCol = chainNbrs[0].x;
						currRow = chainNbrs[0].y;
						VT_HR_EXIT( l2.push_back(edgeIndexImg.Pix(currCol, currRow)) );
					}
					else
						break;
				}

				// Copy the two lists, reversing one
				VT_HR_EXIT(curve.resize(l1.size() + l2.size() - 1));
				int m1 = (int) l1.size();
				int m2 = (int) l2.size();
				for (int j2 = 0; j2 < m2; j2++)
				{
					curve[j2] = l2[m2-1-j2];
				}
				for (int j1 = 1; j1 < m1; j1++)
				{
					curve[m2-1+j1] = l1[j1];
				}
			}

			// Add to curves if we pass the minimum edge segment count 
			if((int) curve.size() > ldParams.minEdgelCount)
			{
				VT_HR_EXIT( curves.push_back(curve) );
			}
		}
	}
	VT_HR_END()
}


HRESULT
	GetNeighbors(vector<CPoint>& nbrs, const vector<EdgeSegment>& edgelList, const CIntImg &edgeIndexImg, CByteImg &status, 
			int col, int row, const LineDetectorParams &ldParams, bool usedZeroCrossings)
{
	VT_HR_BEGIN()
	const EdgeSegment &curre = edgelList[edgeIndexImg.Pix(col, row)];
	int width = edgeIndexImg.Width();
	int height = edgeIndexImg.Height();

	// Using zero crossings - check for neighbouring end points
	if(usedZeroCrossings)
	{
		// endpoint1
		{
			float ep1x = curre.endpoint1.x;
			float ep1y = curre.endpoint1.y;

			// get neighbouring pixel position
			int x0 = (int) ep1x;
			int y0 = (int) ep1y;

			// Get the correct position of neighbor to look
			if(x0 == col && y0 == row)
			{
				if(floorf(ep1x) == ep1x)
					x0--;
				else
					y0--;
			}

			if(x0 > 0 && x0 < width && y0 > 0 && y0 < height)
			{
				// if pixel has an edge segment and status is not visited
				if(edgeIndexImg.Pix(x0, y0) >= 0 && status.Pix(x0, y0) == 0)
				{
					const EdgeSegment &nbre = edgelList[edgeIndexImg.Pix(x0, y0)];
					if(curre.endpoint1 == nbre.endpoint1 || curre.endpoint1 == nbre.endpoint2)
					{
						CPoint pt(x0, y0);
						VT_HR_EXIT( nbrs.push_back(pt) );
					}
				}
			}
		}

		// endpoint 2
		{
			float ep2x = curre.endpoint2.x;
			float ep2y = curre.endpoint2.y;

			// get neighbouring pixel position
			int x1 = (int) ep2x;
			int y1 = (int) ep2y;
			// Get the correct position of neighbor to look
			if(x1 == col && y1 == row)
			{
				if(floorf(ep2x) == ep2x)
					x1--;
				else
					y1--;
			}

			if(x1 > 0 && x1 < width && y1 > 0 && y1 < height)
			{
				// if pixel has an edge segment and status is not visited
				if(edgeIndexImg.Pix(x1, y1) >= 0 && status.Pix(x1, y1) == 0)
				{
					const EdgeSegment &nbre = edgelList[edgeIndexImg.Pix(x1, y1)];
					if(curre.endpoint2 == nbre.endpoint1 || curre.endpoint2 == nbre.endpoint2)
					{
						CPoint pt(x1, y1);
						VT_HR_EXIT( nbrs.push_back(pt) );
					}
				}
			}
		}
	}
	else // not using zero crossings - explore neighbourhood
	{
		// Explore the 3x3 neighborhood and push all potential (unassigned) edgels onto the queue
		static int nbs[8][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1},
		{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
		
		// Explore neighbors 
		for (int i = 0; i < 8; i++)
		{
			int c2 = col + nbs[i][0], r2 = row + nbs[i][1];
			if(c2 > 0 && c2 < width && r2 > 0 && r2 < height)
			{
				int edindex = edgeIndexImg.Pix(c2, r2);
				// if pixel has an edge segment and status is not visited
				if(edindex >= 0 && status.Pix(c2, r2) == 0)
				{
					const EdgeSegment &nbre = edgelList[edgeIndexImg.Pix(c2, r2)];
					CVec2f v(curre.x - nbre.x, curre.y - nbre.y);
					CVec2f u = v.Unit();

					float currd = fabs(u.x * curre.n_x + u.y * curre.n_y);
					float nbrd = fabs(u.x * nbre.n_x + u.y * nbre.n_y);

					if(currd > ldParams.edgeOrientationThreshold &&
						nbrd > ldParams.edgeOrientationThreshold)
					{
						CPoint pt(c2, r2);
						VT_HR_EXIT( nbrs.push_back(pt) );
					}
				}
			}
		}
	}
	VT_HR_END()
}
#define PI 3.14159265358979323846
void ComputeCurvatureofCurves(vector< vector<int> >& curves, vector<EdgeSegment>& edgelList)
{
	int length = (int) edgelList.size();
	vector<float> smooth_n_x;
	vector<float> smooth_n_y;
	smooth_n_x.resize(length);
	smooth_n_y.resize(length);

	for(int i = 0; i < length; ++i)
	{
		smooth_n_x[i] = edgelList[i].n_x;
		smooth_n_y[i] = edgelList[i].n_y;
	}

	// normals are quite noisy, smooth them out
	// smoothing iterations 
	int smoothingIterations = 4;

	// normals are quite noisy, smooth them out
	for(int i = 0; i < (int) curves.size(); ++i)
	{
		vector<int>& curve = curves[i];
		for(int b=0; b < smoothingIterations; ++b)
		{
			for(int j = 1; j < (int) curve.size(); ++j)
			{
				int jcm1 = curve[j-1];
				int jc = curve[j];
				smooth_n_x[jc] = (smooth_n_x[jcm1] + smooth_n_x[jc])/2;
				smooth_n_y[jc] = (smooth_n_y[jcm1] + smooth_n_y[jc])/2;
			}
		}
	}

	// for each curve
	for(int i = 0; i < (int) curves.size(); ++i)
	{
		vector<int>& curve = curves[i];
		float max = 0;
		// compute the curvature for each point
		for(int j = 1; j < (int) curve.size(); ++j)
		{
			int jcm1 = curve[j-1];
			int jc = curve[j];
			EdgeSegment& preve = edgelList[jcm1];
			EdgeSegment& curre = edgelList[jc];
			//EdgeSegment &actuale = edgelList[curve[j]];

			// line joining the two points
			CVec2f dVec(preve.x - curre.x, preve.y - curre.y);
			dVec = dVec.Unit();

			// dot product (in the code from zero crossings the normals are actually
			// rotated by 90 deg, so account for it)
			float theta1 = acosf(dVec.x * (smooth_n_y[jcm1]) - dVec.y * (smooth_n_x[jcm1]));
			float theta2 = acosf(dVec.x * (smooth_n_y[jc]) - dVec.y * (smooth_n_x[jc]));

			// Take the difference and compute curvature = 1/r = sin(theta)/||d|| but ||d|| = 1
			curre.curvature = sinf(theta1 - theta2) * 1000; // scale all the curvatures

			float absCur = fabs(curre.curvature);
			if(max < absCur)
				max = absCur;
		}

		// normalize the curvature
		for(int j = 1; j < (int) curve.size(); ++j)
		{
			EdgeSegment &actuale = edgelList[curve[j]];
			actuale.curvature /= max;
		}
	}
}

inline UINT fsign(float v) { return ((*(UINT*)&v) & 0x80000000) >> 31; }

HRESULT 
    vt::VtBreakCurveBasedOnCurvature(vector< vector<int> >& outputCurves, const vector<EdgeSegment>& edgelList, 
                                 const vector< vector<int> >& curves, float threshold)
{
	VT_HR_BEGIN()

        // dont break small curves
        int minPointsCount = 5;

	// go through all the curves
	for(int i = 0; i < (int) curves.size(); ++i)
	{
		const vector<int> &curve = curves[i];
        if((int) curve.size() < minPointsCount)
            continue;
		vector<int> outputCurve;

		EdgeSegment prev = edgelList[curve[1]];
		// first two values cannot be compared (no curvature for first one
		// and we need a previous to compare), just push them
		VT_HR_EXIT( outputCurve.push_back(curve[0]) );
		VT_HR_EXIT( outputCurve.push_back(curve[1]) );

        int currLength;

		// from the next point
		for(int j = 2; j < (int) curve.size(); ++j)
		{
			const EdgeSegment &curr = edgelList[curve[j]];
			// check for a sign change in the curvature
			UINT e1 = fsign(prev.curvature);
			UINT e2 = fsign(curr.curvature);
			int e3 = e1 ^ e2;

			currLength = (int) outputCurve.size();
			// if a long curve and the curvature values are greater than a threshold and a sign change is detected
			if(currLength > minPointsCount && fabs(prev.curvature) > threshold && fabs(curr.curvature) > threshold && e3)
			{
				// push the current point
				VT_HR_EXIT( outputCurve.push_back(curve[j]) );
				//push the current curve
				VT_HR_EXIT( outputCurves.push_back(outputCurve) );
				// reset the curve
				outputCurve.clear();
			}
			else // no need to create new curve, just add the point
			{	
				VT_HR_EXIT( outputCurve.push_back(curve[j]) );
			}

			// keep track of a previous curvature greater than a 
			// threshold
			if(fabs(curr.curvature) > threshold)
				prev = edgelList[curve[j]];
		}
        currLength = (int) outputCurve.size();
		// push back the curve
        if(currLength > minPointsCount)
		    VT_HR_EXIT( outputCurves.push_back(outputCurve) );
	}

	VT_HR_END()
}

HRESULT
	vt::VtFindLinesFromCurves(vector<LineSegment>& lineSegments, const vector<EdgeSegment>& edgelList, 
						const vector< vector<int> >& curves, const LineDetectorParams& ldParams)
{
	VT_HR_BEGIN()
	
	VT_HR_EXIT( lineSegments.resize(0) );
	// loop through each curve and line simplify the curves
	for(int i = 0; i < (int) curves.size(); ++i)
	{
		const vector<int>& curveSeg = curves[i];
		vector<LineSegment> lineSegmentsPerCurve;

		VT_HR_EXIT( LineSimplify(lineSegmentsPerCurve, curveSeg, edgelList, ldParams) );
		for(int p = 0; p < (int) lineSegmentsPerCurve.size(); ++p)
		{
			VT_HR_EXIT( lineSegments.push_back(lineSegmentsPerCurve[p]) );
		}
	}

	// NOT USING MERGE FOR NOW //
	/*int nmerged;
	VT_HR_EXIT( MergeLineSegments(LineSegments, tmp, nmerged, ldParams) );
	printf("merged segments - %d\n", LineSegments.size());*/

	VT_HR_END()
}

HRESULT
	LineSimplify(vector<LineSegment>& lineCollection, const vector<int>& curveSeg, const vector<EdgeSegment>& edgelList, const LineDetectorParams& ldParams)
{
	VT_HR_BEGIN()
		float maxDist = 0.f;
	int index = 0;
	int length = (int) curveSeg.size();

	// End points of the curve
	const EdgeSegment& edgel1 = edgelList[curveSeg[0]];
	float x1 = edgel1.x;
	float y1 = edgel1.y;

	const EdgeSegment& edgel2 = edgelList[curveSeg[length-1]];
	float x2 = edgel2.x;
	float y2 = edgel2.y;

	for(int k = 1; k < length; ++k)
	{
		// query Point - check to see if this point is the furthest from the line joining
		// the two end points (perpendicular distance from point to line)
		const EdgeSegment& edgelq = edgelList[curveSeg[k]];
		float qx =  edgelq.x;
		float qy =  edgelq.y;

		// slope of the line joining the two end points
		float slope = (y2-y1)/(x2-x1);
		// Perpendicular distance of a point from a line in terms of the points and the slope
		float distance = fabs((qy - slope*(qx - x1) - y1)) / sqrt(slope*slope + 1);

		// Get the index of the point and distance
		if(distance > maxDist)
		{
			index = k;
			maxDist = distance;
		}
	}

	// if the perpendicular distance is less than a certain threshold, then we do not have
	// to split the lines any further. If greater than threshold and min edgel count then
	// split the curve into two parts at the furthest point and recurse on the sub-parts
	if(maxDist >= ldParams.lineSimplificationThreshold && length > ldParams.minEdgelCount)
	{
		// create sub curve 1
		vector<int> curveSeg1;
		VT_HR_EXIT( curveSeg1.reserve(index) );
		for(int j = 0; j < index; ++j)
		{
			VT_HR_EXIT( curveSeg1.push_back(curveSeg[j]) );
		}

		// create sub curve 2
		vector<int> curveSeg2;
		curveSeg2.reserve(length - index);
		for(int j = index; j < length; ++j)
		{
			VT_HR_EXIT( curveSeg2.push_back(curveSeg[j]) );
		}

		// Recursively call line simplify on the two sub parts (divide)
		vector<LineSegment> lineSeg1;
		vector<LineSegment> lineSeg2;
		VT_HR_EXIT( LineSimplify(lineSeg1, curveSeg1, edgelList, ldParams) );
		VT_HR_EXIT( LineSimplify(lineSeg2, curveSeg2, edgelList, ldParams) );
		VT_HR_EXIT( lineCollection.resize(lineSeg1.size() + lineSeg2.size()) );

		// Finally merge the two sub parts together (conquer)
		int m1 = (int) lineSeg1.size();
		int m2 = (int) lineSeg2.size();

		for(int p = 0; p < m1; ++p)
		{
			lineCollection[p] = lineSeg1[p];
		}
		for(int p = 0; p < m2; ++p)
		{
			lineCollection[m1+p] = lineSeg2[p];
		}
	}
	else // we do not have to simplify any further, report results
	{
		// Only end points are stored
		const EdgeSegment &e1 = edgelList[curveSeg[0]];
		const EdgeSegment &e2 = edgelList[curveSeg[length-1]];

		CVec2f pt1;
		pt1.x = e1.x;
		pt1.y = e1.y;

		CVec2f pt2;
		pt2.x = e2.x;
		pt2.y = e2.y;

		// Create a line segment
		LineSegment lseg;

		lseg.start = pt1;
		lseg.end = pt2;
		lseg.length = (float) (pt1 - pt2).Magnitude();

		// store the line segment if length is greater than a min
		if(lseg.length >= ldParams.minLineSegmentLength)
			VT_HR_EXIT( lineCollection.push_back(lseg) );
	}
	VT_HR_END()
}



HRESULT
	MergeLineSegments(vector<LineSegment> &mergedLines, vt::vector<LineSegment>& lines, int& nmerged, LineDetectorParams &ldParams)
{
	VT_HR_BEGIN()

		int linesSize = (int) lines.size();

	// Compute Disjoint Sets of vector<int>-Segments
	vector<vector<int> > sets;
	VT_HR_EXIT( ComputeDSets(lines, sets, ldParams.lineSegmentMergingDistThreshold, ldParams.lineSegmentMergingAngleThreshold) );

	VT_HR_EXIT( mergedLines.reserve(linesSize) );
	int nsets = (int) sets.size();

	// Re-create the list of merged line segments
	for (int j=0;j<nsets;j++)
	{
		int ns = (int) sets[j].size();
		if (ns==1)
		{
			LineSegment lsc;
			int idx = sets[j][0];
			lsc.start   = lines[idx].start;
			lsc.end   = lines[idx].end;
			lsc.length = lines[idx].length;
			VT_HR_EXIT( mergedLines.push_back(lsc) );
		}
		else if (ns > 1)
		{
			LineSegment mergedSegment;
			CMtxd mat;     
			VT_HR_EXIT( mat.Create(2*ns, 3) );
			int   c = 0;
			for (int k=0;k<ns;k++)
			{
				LineSegment lc = lines[sets[j][k]];
				mat(c,0)   = lc.start.x; 
				mat(c,1)   = lc.start.y; 
				mat(c,2)   = 1.0;    
				c++;
				mat(c,0)   = lc.end.x; 
				mat(c,1)   = lc.end.y; 
				mat(c,2)   = 1.0;    
				c++;
			}
			CSolveSVDd mat_svd(mat);
			CVecd line_eq = mat_svd.GetBestNullSpaceVector();            
			vector<CVec2d> points; 
			VT_HR_EXIT( points.reserve(2*ns) );

			CVec2d center;
			center.x = 0;
			center.y = 0;

			for (int k=0;k<ns;k++)
			{
				CVec2d projp1, projp2;
				LineSegment lc = lines[sets[j][k]];
				ProjectOnLine(projp1, line_eq, lc.start);
				center = center + projp1;    
				VT_HR_EXIT( points.push_back(projp1) ) ;
				ProjectOnLine(projp2, line_eq, lc.end);    
				VT_HR_EXIT( points.push_back(projp2) );
				center = center + projp2;
			}
			center = (1.0/2*ns)*center;
			CVec2d ls, le;
			CVec2d lv = (points[0]-center);
			double maxdp = 0, mindp = 0.0;

			for (int k=0;k<(int)points.size();k++)
			{
				CVec2d vv   = points[k]-center;
				double dotp = vv * lv; 
				if (dotp > 0)
				{
					if (fabs(dotp) > maxdp)    {
						maxdp = fabs(dotp);
						le = points[k];
					}
				}
				else
				{
					if (fabs(dotp) > mindp)    {
						mindp = fabs(dotp);
						ls = points[k];
					}
				}
			}
			mergedSegment.start   = ls;
			mergedSegment.end   = le;
			mergedSegment.length = (float) (ls-le).Magnitude();
			VT_HR_EXIT( mergedLines.push_back(mergedSegment) );
		}
	}
	nmerged = (int) mergedLines.size();

	VT_HR_END()
}



HRESULT
	ComputeDSets(vector<LineSegment> &data, vector<vector<int> > &sets, float distThreshold, float angleThreshold)
{
	VT_HR_BEGIN()

		int  sz = (int) data.size();
	DisjointSets  ds;
	ds.Init(sz);
	vector<int>  setIds;
	for (int i=0;i<sz;i++)
	{
		for (int j=i+1;j<sz;j++)
		{
			bool value;
			VT_HR_EXIT( LineSegmentsAreCollinearAndClose(value, distThreshold, angleThreshold, data[i], data[j]) );
			if (value)
			{
				if (ds.FindSet(i)!=ds.FindSet(j))
					ds.Union(i,j);
			}
		}
	}
	ds.SetIds(setIds);
	int m=0;
	for (int i=0;i<(int)setIds.size();i++)
		m = VtMax(m,setIds[i]);
	vector<vector<int> > st;
	VT_HR_EXIT( st.resize(m+1) );
	VT_HR_EXIT( sets.reserve(sz) );

	for (int i=0;i<(int)setIds.size();i++)
	{
		int setId = setIds[i];
		VT_HR_EXIT( st[setId].push_back(i) );
	}
	for (int i=0;i<(int)st.size();i++)
	{
		if (st[i].size()>0)
			VT_HR_EXIT( sets.push_back(st[i]) );
	}
	VT_HR_END()
}



HRESULT 
	LineSegmentsAreCollinearAndClose(bool &value, float distThreshold, float angleThreshold, LineSegment &a, LineSegment &b)
{
	VT_HR_BEGIN()

	value = true;
	double len   = (double) (a.length+b.length);
	double lensq = len*len;
	double d1, d2, d3, d4, dmin, dmax;
	CVec2d v1, v2, v3, v4, va, vb, vsm;
	v1 = (a.start-b.start);
	v2 = (a.end-b.end);
	v3 = (a.start-b.end);
	v4 = (a.end-b.start);
	va = (a.start-a.end);
	vb = (b.start-b.end);
	double dotp = fabs((va*vb)/double(a.length * b.length));
	if (dotp > angleThreshold)
	{
		d1 = v1.MagnitudeSq();
		d2 = v2.MagnitudeSq();
		d3 = v3.MagnitudeSq();
		d4 = v4.MagnitudeSq();
		dmin = VtMin(d1,d2); dmin = VtMin(dmin,d3); dmin = VtMin(dmin,d4);
		dmax = VtMax(d1,d2); dmax = VtMax(dmax,d3); dmax = VtMax(dmax,d4);
		if (dmin > 0.04*lensq)
		{
			value = false;
			VT_HR_EXIT(S_OK);
		}
		if (dmin > a.length*a.length || dmin > b.length*b.length)
		{
			value = false;
			VT_HR_EXIT(S_OK);
		}
		if (dmax > lensq)
		{
			if (dmax==d1)         vsm = v1.Unit();
			else if (dmax==d2)    vsm = v2.Unit();
			else if (dmax==d3)    vsm = v3.Unit();
			else if (dmax==d4)    vsm = v4.Unit();

			if (fabs(vsm*va)/a.length > angleThreshold && fabs(vsm*vb)/b.length > angleThreshold )
			{
				CVec3d aeq, beq;
				CVec2d asp, aep, bsp, bep;
				aeq = (a.start.Hom()).Cross(a.end.Hom());
				beq = (b.start.Hom()).Cross(b.end.Hom());

				CVecd aeqv;
				aeqv.Create(3);
				CVecd beqv;
				beqv.Create(3);
				aeqv[0] = aeq.x;
				aeqv[1] = aeq.y;
				aeqv[2] = aeq.z;

				beqv[0] = beq.x;
				beqv[1] = beq.y;
				beqv[2] = beq.z;

				VT_HR_EXIT(ProjectOnLine(bsp, aeqv, b.start));
				VT_HR_EXIT(ProjectOnLine(bep, aeqv, b.end));
				VT_HR_EXIT(ProjectOnLine(asp, beqv, a.start));
				VT_HR_EXIT(ProjectOnLine(aep, beqv, a.end) );
				if ((a.start-asp).MagnitudeSq() < distThreshold && (a.end-aep).MagnitudeSq() < distThreshold &&
					(b.start-bsp).MagnitudeSq() < distThreshold && (b.end-bep).MagnitudeSq() < distThreshold)
				{
					value = true;
					VT_HR_EXIT(S_OK);
				}
				else
				{
					value = false;
					VT_HR_EXIT(S_OK);
				}
			}
			else
			{
				value = false;
				VT_HR_EXIT(S_OK);
			}
		}
		else
		{
			value = false;
			VT_HR_EXIT(S_OK);
		}
	}
	else
	{
		value = false;
		VT_HR_EXIT(S_OK);
	}

	VT_HR_END()
}

HRESULT
	ProjectOnLine(CVec2d& proj, CVecd& leq, CVec2d& p)
{
	VT_HR_BEGIN()
		if(leq.Size() < 3)
			VT_HR_EXIT(E_INVALIDARG);
	double a,b,c,d, dn;
	a = leq[0];
	b = leq[1];
	c = leq[2];
	d = a*p[1]-b*p[0];
	dn= a*a + b*b;
	proj[0] = (-b*d - a*c)/dn;
	proj[1] = (a*d - b*c)/dn;

	VT_HR_END()
}



