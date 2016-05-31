//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      EllipseSegment fitting from linked edge lists
//
//  History:
//      2011/07/07-szeliski
//          Created from earlier Sho/Python port of AWF's Matlab code
//          The comment lines with Python code (no terminating ;) are copied from fitellipse.py
//
//  AWF's original description/header:
//
//  FITELLIPSE  Least-squares fit of ellipse to 2D points.
//         A = FITELLIPSE(X,Y) returns the parameters of the best-fit
//         ellipse to 2D points (X,Y).
//         The returned vector A contains the center, radii, and orientation
//         of the ellipse, stored as (Cx, Cy, Rx, Ry, theta_radians)
//
//  Authors: Andrew Fitzgibbon, Maurizio Pilu, Bob Fisher
//  Reference: "Direct Least Squares Fitting of Ellipses", IEEE T-PAMI, 1999
//
//   @Article{Fitzgibbon99,
//    author = "Fitzgibbon, A.~W.and Pilu, M. and Fisher, R.~B.",
//    title = "Direct least-squares fitting of ellipses",
//    journal = pami,
//    year = 1999,
//    volume = 21,
//    number = 5,
//    month = may,
//    pages = "476--480"
//   }
//  
//  This is a more bulletproof version than that in the paper, incorporating
//  scaling to reduce roundoff error, correction of behaviour when the input 
//  data are on a perfect hyperbola, and returns the geometric parameters
//  of the ellipse, rather than the coefficients of the quadratic form.
//------------------------------------------------------------------------

#include "stdafx.h"
#include "vt_fitEllipse.h"

using namespace vt;

#define M_PI       3.14159265358979323846


//Local function declarations **** //

/// <summary> Fit the ellipse to the given edge points. </summary>
/// <param name="l"> List of edge indices </param>
/// <param name="e"> Returned ellipse parameters </param>
HRESULT Fit(const vector<int>& l, EllipseSegment& e, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams);

/// <summary> Filter out the edgels that aren't close enough to the fit. </summary>
/// <param name="l"> List of current edge indices </param>
/// <param name="e"> Current ellipse parameters (reset n_points)</param>
/// <param name="lf"> New filtered set of edge indices </param>
HRESULT FilterPoints(const vector<int>& l, EllipseSegment& e, vector<int>& lf, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams);

/// <summary> Iteratively re-fit the ellipse to the inlier edge points. </summary>
/// <param name="l"> List of edge indices </param>
/// <param name="e"> Final fitted ellipse parameters and inlier points </param>
HRESULT RobustFit(const vector<int>& l, EllipseSegment& e, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams);

/// <summary> Determine the min_t and max_t for ellipse segment and reverse order of points, if necessary. </summary>
/// <param name="e"> EllipseSegment parameters and inlier points </param>
void OrderEllipse(EllipseSegment &e, const vector<EdgeSegment> &edgelList);

/// <summary> Merge ellipses and keep the best ones </summary>
/// <param name="ellipses"> List of fitted ellipses </param>
HRESULT MergeEllipses(vector<EllipseSegment>& ellipses, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams);

/// <summary> Generate a synthetic noisy example and test the accuracy of the fit. </summary>
/// <param name="noise_level"> Standard deviation of noise added to ellipse </param>
/// <param name="extent"> Angular extent of test ellipse (full ellipse is 2 PI radians) </param>
/// <param name="n_points"> Number of generated points </param>
HRESULT UnitTest(float noiseLevel = 0.5f, float extent = 2.0f, int nPoints = 100);

// ****  local function declarations END//


// TODO: These are helper functions for the vector class, and should be moved there:
//  These versions assume that data types are elemental, i.e., they can be memcpy'd
template <typename T>
HRESULT FastAppend(vector<T>& dst, const vector<T>& src0, const vector<T>& src1)
{
	VT_HR_BEGIN()
    int n0 = (int) src0.size();
    int n1 = (&src1) ? (int) src1.size() : 0;
    VT_HR_EXIT( dst.resize(n0 + n1) );
    memcpy(&dst[0], &src0[0], n0*sizeof(T));
    if (&src1)
        memcpy(&dst[n0], &src1[0], n1*sizeof(T));
    VT_HR_END()
}

template <typename T>
HRESULT FastCopy(vector<T>& dst, const vector<T>& src0)
{
    return FastAppend(dst, src0, *(vector<T> *) 0);
}

HRESULT 
	Fit(const vector<int>& l, EllipseSegment& ellipse, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams)
{
    // Fit the ellipse to the given edge points
	VT_HR_BEGIN()
    memset(&ellipse, 0, sizeof(EllipseSegment));    // initialize to all 0 values (invalid)

	// Normalize the data (compute the extents min, max, mean)
	const EdgeSegment &e0 = edgelList[l[0]];
	float x_min = e0.x, x_max = e0.x, mx = 0;
	float y_min = e0.y, y_max = e0.y, my = 0;
	int n = (int) l.size();

	if (n < feParams.min_points)
		return S_OK;
	for (int i = 0; i < n; i++)
    {
        const EdgeSegment &e = edgelList[l[i]];
        x_min   = __min(x_min, e.x);
        x_max   = __max(x_max, e.x);
        mx += e.x;
        y_min   = __min(y_min, e.y);
        y_max   = __max(y_max, e.y);
        my += e.y;
    }
    mx /= n;
    my /= n;
    float sx = (x_max - x_min)/2;
    float sy = (y_max - y_min)/2;
    if (sx < feParams.min_extent || sy < feParams.min_extent)
        return S_OK;

    // Build design matrix
    //  D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ]
    // and the scatter matrix
    //  S = D.T * D
    double S[6][6];
    memset(S, 0, sizeof(S));
    double sxi = 1.0 / sx;
    double syi = 1.0 / sy;
    for (int i = 0; i < n; i++)
    {
        const EdgeSegment &e = edgelList[l[i]];
        double x  = sxi*(e.x - mx);
        double y  = syi*(e.y - my);
        double xx = x*x;
        double xy = x*y;
        double yy = y*y;
        double D[6] = {xx, xy, yy, x, y, 1};
        for (int j = 0; j < 6; j++)
        {
            for (int k = 0; k < 6; k++)
            {
                S[j][k] += D[j] * D[k];
            }
        }
    }
    //  print S

    // Build the 6x6 constraint matrix
    //  C = DoubleArray(6, 6)
    //  C[0,2], C[2,0], C[1,1] = -2, -2, 1
    //
    // Set up the eigensystem
    // Break into blocks
    //  tmpA = S[0:3,0:3]
    //  tmpB = S[0:3,3:6]
    //  tmpC = S[3:6,3:6]
    //  tmpD = C[0:3,0:3]
    CMtx3x3<double> tmpA, tmpB, tmpC;
    for (int j = 0; j < 3; j++)
    {
        tmpA.SetRow(j, *(const CVec3<double> *) &S[j+0][0]);
        tmpB.SetRow(j, *(const CVec3<double> *) &S[j+0][3]);
        tmpC.SetRow(j, *(const CVec3<double> *) &S[j+3][3]);
    }
    CMtx3x3<double> tmpD(0, 0, -2, 0, 1, 0, -2, 0, 0);
    CMtx3x3<double> tmpE = tmpC.Inv() * tmpB.T();
    CMtx3x3<double> tmpF = tmpB * tmpE;
    CMtx3x3<double> tmpG = 0.5 * (tmpF + tmpF.T());
    CMtx3x3<double> tmpH = tmpD.Inv() * (tmpA - tmpG);

    // Solve eigensystem
    //  eigv = Eigen(tmpH)
    //  print tmpA, tmpB, tmpC, tmpD, tmpE, tmpF, tmpG
    //  evcx = eigv.V.Real()
    //  evlx = eigv.D.Real()
    // print "eval_x, evec_x", evlx, evcx
    CMtxd mtx(tmpH);
    VT_HR_EXIT( mtx.GetError() );
    CSolveEigen<double> eigv;
    VT_HR_EXIT( eigv.Solve(mtx) );
    CMtxd evcx = VtRe(eigv.V());
    CMtxd evlx = VtRe(eigv.D());
    VT_HR_EXIT( evcx.GetError() );
    VT_HR_EXIT( evlx.GetError() );
    // Aargh!  This error checking is so ugly...wish all PGs would allow exceptions...

    //
    // Find the positive (as det(tmpD) < 0) eigenvalue
    //  sevs = evlx.Sort()
    //  minv = sevs[0,0]
    //  fmnv = list(Find(evlx <= minv))
    //  I    = fmnv[0].Row
    //  print evlx, I
    int I = (evlx[0][0] < evlx[1][1]) ?
        ((evlx[0][0] < evlx[2][2]) ? 0 : 2) :
        ((evlx[1][1] < evlx[2][2]) ? 1 : 2);
    //
    // Extract eigenvector corresponding to negative eigenvalue
    CVec3<double> A1(evcx.GetCol(I));
    //
    // Recover the bottom half...
    CVec3<double> A2 = -tmpE * A1;
    //  A = DoubleArray.VertStack(A1, A2)
    double A[6] = {A1[0], A1[1], A1[2], A2[0], A2[1], A2[2]};
    // print A
    //
    // Unnormalize
    double par[6] = {
        A[0]*sy*sy,
        A[1]*sx*sy,
        A[2]*sx*sx,
        -2*A[0]*sy*sy*mx - A[1]*sx*sy*my + A[3]*sx*sy*sy,
        -A[1]*sx*sy*mx - 2*A[2]*sx*sx*my + A[4]*sx*sx*sy,
        A[0]*sy*sy*mx*mx  + A[1]*sx*sy*mx*my + A[2]*sx*sx*my*my
        -A[3]*sx*sy*sy*mx - A[4]*sx*sx*sy*my + A[5]*sx*sx*sy*sy
    };
    //  print 'par =', par
    // New code (added by Rick 2011-07-21) to force par[0]-par[2] > 0
    //  so that thetarad is in [-45,45]
    if (par[0] < par[2])
    {
        for (int k = 0; k < 6; k++)
        {
            par[k] = -par[k];
        }
    }
    //
    // Convert to geometric radii, and centers
    //
    double thetarad = 0.5 * atan2(par[1], par[0] - par[2]);
    double cost, sint;
    VtSinCos(thetarad, &sint, &cost);
    double sin2 = sint * sint;
    double cos2 = cost * cost;
    double cssn = sint * cost;
    //
    double Ao  =  par[5];
    double Au  =  par[3] * cost + par[4] * sint;
    double Av  = -par[3] * sint + par[4] * cost;
    double Auu =  par[0] * cos2 + par[2] * sin2 + par[1] * cssn;
    double Avv =  par[0] * sin2 + par[2] * cos2 - par[1] * cssn;
    //
    // ROTATED = [Ao Au Av Auu Avv]
    // print 'Ao, Au, Av, Auu, Avv = ', Ao, Au, Av, Auu, Avv
    //
    double tuCentre = -Au/(2*Auu);
    double tvCentre = -Av/(2*Avv);
    double wCentre  =  Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;
    //
    double uCentre = tuCentre * cost - tvCentre * sint;
    double vCentre = tuCentre * sint + tvCentre * cost;
    //
    double Ru2 = -wCentre/Auu;
    double Rv2 = -wCentre/Avv;
    //  I can't think of why the sign actually matters, so drop it...
    //  Ru = sqrt(abs(Ru2))*sign(Ru2);
    //  Rv = sqrt(abs(Rv2))*sign(Rv2);
    double Ru = sqrt(abs(Ru2));
    double Rv = sqrt(abs(Rv2));
    //
    // Save the parameters
    ellipse.cx = (float) uCentre;
    ellipse.cy = (float) vCentre;
    ellipse.ax = (float) Ru;
    ellipse.ay = (float) Rv;
    ellipse.theta = (float) thetarad;
    ellipse.n_points = n;
    VT_HR_END()
}

// Generate sample points on the ellipse
HRESULT 
	EllipseSegment::GeneratePoints(int n_pts, vector<EdgeSegment>& edgels, vector<int>& l, float noiseLevel)
{
	VT_HR_BEGIN()
    // If no points specified, compute from axes and extent
    n_points = n_pts;
    if (n_points <= 0)
    {
        float fudge = 1.1f;
        n_points = (int) ceil(fudge * (max_t-min_t) * __max(ax, ay));
    }
    VT_HR_EXIT( edgels.resize(n_points) );
    VT_HR_EXIT( l.resize(n_points) );
    memset(&edgels[0], 0, n_points*sizeof(EdgeSegment));
    float dt = (max_t - min_t) / (n_points-1);
    for (int i = 0; i < n_points; i++)
    {
        EdgeSegment& e = edgels[i];
        float t = min_t + dt*i;
        float u = ax*cos(t);
        float v = ay*sin(t);
        e.x = u*cos(theta) - v*sin(theta) + cx;
        e.y = u*sin(theta) + v*cos(theta) + cy;
        // TODO: add the noise (need to figure out where Gaussian/normal noise is in VisionTools)
        // TODO:  fill in the normal information, etc.
		noiseLevel = 0;

        l[i] = i;
    }
    VT_HR_END()
}

// Generate a synthetic noisy example and test the accuracy of the fit
HRESULT 
	UnitTest(float noise_level, float extent, int n_points)
{
	VT_HR_BEGIN()
    // Create an ellipse
    EllipseSegment test_ellipse = {250.f, 150.f, 300.f, 200.f, 0.4f, 0, extent};
    vector<EdgeSegment> edgels;
    vector<int> l;
    VT_HR_EXIT( test_ellipse.GeneratePoints(n_points, edgels, l, noise_level) );

    // Try usingthe default fitting parameters
    FitEllipseParams par;

    // Now fit to the noisy samples
    EllipseSegment estimated_ellipse;
    VT_HR_EXIT( Fit(l, estimated_ellipse, edgels, par) );

    VT_HR_END()
}

//
// New code added to original fitellipse.m code in building up robust fitellipse.py code
//  2011-05-27-szeliski
//

//
// Estimate distance from ellipse and angle along arc of ellipse ('t'),
//  as well as the angle between the edgel direction and the ellipse tangent
//
void 
	EllipseSegment::td2aFromXY(const EdgeSegment& edgel, float& t, float& d2, float& a) const
{
    float x  = edgel.x;
    float y  = edgel.y;
    float nx = edgel.n_x;
    float ny = edgel.n_y;
    float dx = x - cx;              // distance from center
    float dy = y - cy;
    float cosr, sinr;       
    VtSinCos(theta, &sinr, &cosr);    // rotation of ellipse
    float rx =  dx*cosr + dy*sinr;  // projection onto major/minor axis
    float ry = -dx*sinr + dy*cosr;
    float mx =  nx*cosr + ny*sinr;  // rotated edge normal
    float my = -nx*sinr + ny*cosr;
    // After scaling by radii, distance is no longer truly Euclidean
    //  so the t estimate is not strictly the nearest point, but close
    float sx = rx / ax;             // scaled down by axis length
    float sy = ry / ay;
          t  = atan2(sy, sx);       // estimate of angle along ellipse
    float ex =  ax * cos(t);        // nearest ellipse point
    float ey =  ay * sin(t);
    // Approximate squared distance between point and ellipse:
          d2 = (rx - ex)*(rx - ex) + (ry - ey)*(ry - ey);
    float tx = -ax * sin(t);        // derivative of xy(t) function
    float ty =  ay * cos(t);
    float ts = sqrt(tx*tx + ty*ty); // tangential speed
    float vx = tx / ts;             // tangent direction
    float vy = ty / ts;
          a  = abs(mx*vy - my*vx);  // |sin<edgel dir, ellipse dir>|
    //  print 'x, y, dx, dy, cosr, sinr, rx, ry, sx, sy, t, ex, ey, d2'
    //  print  x, y, dx, dy, cosr, sinr, rx, ry, sx, sy, t, ex, ey, d2
    //  print 'nx, ny, rx, ry, mx, my, tx, ty, vx, vy, a'
    //  print  nx, ny, rx, ry, mx, my, tx, ty, vx, vy, a
}

// Filter out the large error points
HRESULT 
	FilterPoints(const vector<int>& l, EllipseSegment& e, vector<int>& lf, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams)
{
	VT_HR_BEGIN()
    
    float d2Thresh = feParams.max_dist * feParams.max_dist;
    lf.clear();     // initialize returned list to empty
    int n = (int) l.size();
    e.n_points = 0;     // in case we exit prematurely
    float sum_strength = 0;
    float sum_distance = 0;
    for (int i = 0; i < n; i++)
    {
        const EdgeSegment &el = edgelList[l[i]];
        float t, d2, a;     // parameter along ellipse, L2 distance, and angular distance
        e.td2aFromXY(el, t, d2, a);
        if (d2 <= d2Thresh && a <= feParams.max_angle)
        {
            VT_HR_EXIT( lf.push_back(l[i]) );
            if (lf.size() > 1)
            {
                const EdgeSegment &ep = edgelList[lf[(int)lf.size()-2]];
                float dist = (float) hypot(el.x-ep.x, el.y-ep.y);
                sum_distance += dist;
                sum_strength += sqrt(el.strength);
            }
        }
    }
    e.n_points = (int) lf.size();
    e.strength = sum_strength / __max(e.n_points, 1);
    e.density  = sum_distance / __max(e.n_points-1, 1);
    // TODO: The score below is just a heuristic, and could be fine-tuned
    e.score = e.n_points * e.strength / (e.density + 0.1f);
    VT_HR_END()
}

// Determine the min_t and max_t for ellipse segment and reverse order of points, if necessary
void 
	OrderEllipse(EllipseSegment &e, const vector<EdgeSegment> &edgelList)
{
    // Compute the distances and t values
    int n = (int) e.points.size();
    float t0, tp, tc, d2,  a;
	t0 = 0;
	tp = 0;
    int poss = 0;
    int negs = 0;
    for (int i = 0; i < n; i++)
    {
        e.td2aFromXY(edgelList[e.points[i]], tc, d2, a);
        if (i == 0)
        {
            t0 = tc;
        }
        else
        {
            // Count the number of positive and negative adjacent differences
            float delta = tc - tp;
            poss += (delta > 0);
            negs += (delta < 0);
        }
        tp = tc;
    }

    // Set the min and max extents
    if (poss > negs)
    {
        e.min_t = t0;
        e.max_t = tp;
    }
    else
    {
        e.min_t = tp;
        e.max_t = t0;
        // Reverse the order of the points
        for (int i = 0; i < n/2; i++)
        {
            int j = e.points[i];
            e.points[i] = e.points[n-1-i];
            e.points[n-1-i] = j;
        }
    }

    // Make sure that max_t > min_t
    if (e.max_t <= e.min_t)
    {
        if (e.min_t + e.max_t > 0)
        {
            e.min_t -= (float) (2*M_PI);
        }
        else
        {
            e.max_t += (float) (2*M_PI);
        }
    }
}

// Iteratively re-fit the ellipse to the inlier edge points
HRESULT 
	RobustFit(const vector<int>& l, EllipseSegment& e, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams)
{
	VT_HR_BEGIN()
    // Allocate the scratch line buffers and check for errors (always shrink?)
    vector<int> lt;
    VT_HR_EXIT( lt.reserve(l.size()) );
    vector<int>& lf = e.points;    // list of final filtered points
    //  print 'from', n0, 'points keep',
    lt = l;

    // Iteratively re-fit and find inliers
    for (int iter = 0; iter < feParams.n_iter && lt.size() > 0; iter++)
    {
        // Do the current fit
        VT_HR_EXIT( Fit(lt, e, edgelList, feParams) );
        if (e.n_points == 0)
            return hr;

        // Figure out which points are the inliers
        VT_HR_EXIT( FilterPoints(lt, e, lf, edgelList, feParams) );
        if (e.n_points == 0)
            return hr;
        //  print nf,
        lt = lf;
    }

    // Compute the distances and t values
    int n = (int) lf.size();
    vector<float> ts;
    VT_HR_EXIT ( ts.resize(n) );
    float d2,  a;
    for (int i = 0; i < n; i++)
    {
        e.td2aFromXY(edgelList[lf[i]], ts[i], d2, a);
    }

    // Compute the start and end t values and orient the ellipse CCW
    OrderEllipse(e, edgelList);
    // 2011-11-06-szeliski:  add in arc length diff_t as well
    // TODO:  score could be refined even further...
    e.score *= e.max_t-e.min_t;
    e.n_points = n;
    VT_HR_END()
}


// Compare two ellipses to see which one has a longer length/strength product
static int CompareEllipse(const void *p0, const void *p1)
{
    EllipseSegment e0 = *(EllipseSegment *) p0;
    EllipseSegment e1 = *(EllipseSegment *) p1;
    float s0 = e0.score;
    float s1 = e1.score;
    // We want longer curves FIRST, so choose the one with the LARGER score
    return (s0 > s1) ? -1 : (s0 < s1) ? 1 : 0;
}

// Fit ellipses to all the curves and sort by length/strength
HRESULT vt::VtFitEllipsesToCurves(OUT vector<EllipseSegment>& ellipses, 
        const vector< vector<int> > &curves, const vector<EdgeSegment> &edgelList, 
        const FitEllipseParams &feParams)
{
    // Fit each curve individually
    //  TODO:  break each curve up into singly curved regions first
	VT_HR_BEGIN()
    EllipseSegment ell;
    int n = (int) curves.size();
    for (int i = 0; i < n; i++)
    {
        VT_HR_EXIT( RobustFit(curves[i], ell, edgelList, feParams) );
        if (ell.n_points > 0)
        {
            VT_HR_EXIT( ellipses.push_back(ell) );
        }
    }

	// If no ellipses are found then just return pre-emptively
	if(ellipses.size() == 0)
		goto Exit;

    // Sort by length/strength
    qsort(ellipses.begin(), ellipses.size(), sizeof(EllipseSegment), CompareEllipse);

    // Merge ellipses whose endpoints are close
    if (feParams.max_merge_candidates > 0)
    {
        VT_HR_EXIT( MergeEllipses(ellipses, edgelList, feParams) );
    }
    VT_HR_END()
}

// Merge ellipses and keep the best ones
HRESULT 
	MergeEllipses(vector<EllipseSegment>& ellipses, const vector<EdgeSegment> &edgelList, const FitEllipseParams &feParams)
{
	VT_HR_BEGIN()
    // Create a new list
    vector<EllipseSegment> merged;

	// If number of ellipses is less than merge candidates then set the number to be the num of ellipses
	int mmc = ((int) ellipses.size() > feParams.max_merge_candidates) ? 
									feParams.max_merge_candidates : (int) ellipses.size();

    // Iterate over the original ellipse, trying to find match candidates
    for (int i = 0; i < mmc; i++)
    {
        EllipseSegment ei = ellipses[i];      // makes a copy, so points can be clobbered
        if (ei.n_points == 0)
            continue;

        // Test each subsequent ellipse;  if merged, set its count to 0
        for (int j = i+1; j < mmc; j++)
        {
            EllipseSegment& ej = ellipses[j];
            if (ej.n_points == 0)
                continue;

            // See if endpoints are close enough
            const EdgeSegment& ei0 = (edgelList)[ei.points.front()];
            const EdgeSegment& ein = (edgelList)[ei.points.back()];
            const EdgeSegment& ej0 = (edgelList)[ej.points.front()];
            const EdgeSegment& ejn = (edgelList)[ej.points.back()];
            float dist0 = (float) hypot(ei0.x-ejn.x, ei0.y-ejn.y);
            float dist1 = (float) hypot(ein.x-ej0.x, ein.y-ej0.y);
            if (dist0 > feParams.max_endpoint_dist &&
                dist1 > feParams.max_endpoint_dist)
                continue;

            // Merge the two lists
            vector<int>& l0 = (dist1 <  dist0) ? ei.points : ej.points;
            vector<int>& l1 = (dist1 >= dist0) ? ei.points : ej.points;
            vector<int>  merged_points;
            VT_HR_EXIT( FastAppend(merged_points, l0, l1) );

            // Fit the new merged points
            EllipseSegment ell;
            VT_HR_EXIT( RobustFit(merged_points, ell, edgelList, feParams) );

            // See if this is a better fit
            if (ell.n_points <= ei.n_points ||
                ell.n_points < feParams.min_merge_frac * (ei.n_points + ej.n_points))
                continue;

            // Copy the merged list to the first ellipse and "erase" the second one
            ei = ell;   // TODO: unsafe, since vector copy isn't error-checked, but I'm tired of this...
            ej.n_points = 0;
            ej.points.clear();
        }

        // Save the current (potentially merged) ellipse
        VT_HR_EXIT( merged.push_back(ei) );
    }

    // Sort by length/strength
    qsort(merged.begin(), merged.size(), sizeof(EllipseSegment), CompareEllipse);

    // Copy the merged list over the old one
    // TODO:  is this guaranteed to be safe?
    ellipses = merged;
    VT_HR_END()
}


//
//  Bi-tangent to two ellipses
//

void Conic(CMtx3x3d& C, CMtx3x3d& Ci, const EllipseSegment& el)
{
    // Compute the conic section matrix
    //  Note that I use double precision, since we are squaring floating point terms
    //   and there aren't a lot of ellipses being computed and fitted
    CMtx3x3d T, R, D, E;
    T.MakeTranslate(-el.cx, -el.cy);
    R.MakeRotation(-el.theta);
    CVec3d d(1.0f/el.ax, 1.0f/el.ay, 1.0f);
    D.MakeDiag(d);
    CVec3d e(1.0f, 1.0f, -1.0f);
    E.MakeDiag(e);
    CMtx3x3d Tf = D * R * T;
    C = Tf.T() * E * Tf;
    C.MakeSymmetric();

    // Compute its inverse at the same time
    CMtx3x3d Ti, Ri, Di;
    Ti.MakeTranslate(el.cx, el.cy);
    Ri.MakeRotation(el.theta);
    CVec3d di(el.ax, el.ay, 1.0f);
    Di.MakeDiag(di);
    CMtx3x3d Tb = Ti * Ri * Di;
    Ci = Tb * E * Tb.T();
    Ci.MakeSymmetric();
}

// Find 3 cubic roots of a non-degenerate cubic
HRESULT CubicRoots(CVec4d a, double roots[3])
{
	VT_HR_BEGIN();

    // Check that it's not a degenerate cubic
    if (a[3] == 0)
        VT_HR_EXIT( E_INVALIDARG );     // should be a real cubic, not degenerate

    // Find the cubic roots, using the formula from Press section 5.5
    double a1 = a[2]/a[3], a2 = a[1]/a[3], a3 = a[0]/a[3];
    double Q = (a1*a1 - 3*a2) / 9;
    double R = (2*a1*a1*a1 - 9*a1*a2 + 27*a3) / 54;
    if (Q*Q*Q < R*R)
        VT_HR_EXIT( E_INVALIDARG );     // should have 3 real roots
    double SQ = sqrt(Q);
    double theta = acos(R / (SQ*SQ*SQ));
    double s1 = -2 * SQ;
    double t3 = theta / 3;
    double tf = - a1 / 3;
    roots[0] = s1 * cos(t3               ) + tf;
    roots[1] = s1 * cos(t3 + 2 * M_PI / 3) + tf;
    roots[2] = s1 * cos(t3 + 4 * M_PI / 3) + tf;
    VT_HR_END();
}

void LDLTPermuted(CMtx3x3d& A, CVec3d& d, CMtx3x3d& L, bool useAbs, bool takeSqrt)
{
    // A = L D L^T = factorization of a symmetric matrix such that
    //  L is lower triangular with a unit diagonal
    // The actual L matrix returned is not lower triangular, since we use pivoting to put the smallest
    // (potentially 0) element at the bottom of the diagonal
    // If useAbs, sort by magnitude (for rank deficient); else, sort by signed magnitude (for conic)

    // Compute a permutation matrix
    CVec3d Ad(A[0][0], A[1][1], A[2][2]);
    CVec3d Av = (useAbs) ? CVec3d(fabs(Ad[0]), fabs(Ad[1]), fabs(Ad[2])) : Ad;
    int i0 = (Av[0] < Av[1]) ?
            ((Av[1] < Av[2]) ? 2 : 1) : 
            ((Av[0] < Av[2]) ? 2 : 0);
    int i1 = (i0+1)%3, i2 = (i0+2)%3;
    double D1 = A[i1][i1] - A[i0][i1]*A[i1][i0]/A[i0][i0];
    double D2 = A[i2][i2] - A[i0][i2]*A[i2][i0]/A[i0][i0];
    bool pick1 = (useAbs) ? fabs(D1) >= fabs(D2) : D1 >= D2;
    int p[3] = {i0, (pick1) ? i1 : i2, (pick1) ? i2 : i1};   // permutation vector
    CMtx3x3d P;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            P[i][j] = j == p[i];
        }
    }

    //  Call the VtInPlaceLDLSolve routine, but keep it from zeroing out return value on INVALIDARG
    CMtx<double> Am = P * A * P.T();
    VtInPlaceLDLSolve(Am, true);

    // Read out the diagonal and lower-triangular matrix
    CMtx3x3d Lp = (CMtx3x3d) Am;
    d[0] = Lp[0][0], d[1] = Lp[1][1], d[2] = Lp[2][2];
    if (d[0] < 0)
        d = -d;                                     // flip so leading number is positive
    Lp[0][0] = Lp[1][1] = Lp[2][2] = 1;             // set diagonal to 1
    Lp[0][1] = Lp[0][2] = Lp[1][2] = 0;             // set above diagonal to 0

    // Optionally take out the square root diagonal, as in Cholesky
    if (takeSqrt)
    {
        for (int i = 0; i < 3; i++)
        {
            double d2 = sqrt(fabs(d[i]));
            d[i] = (d[i] < 0) ? -1 : (d[i] > 0) ? 1 : 0;
            for (int j = i; j < 3 && d2 != 0; j++)
                Lp[j][i] *= d2;
        }
    }

    // Now apply the final permutation
    L = P.T() * Lp;
}

// Intersect two ellipses, assuming that there are exactly 4 intersections
HRESULT IntersectConics(CVec3d points[4] , CMtx3x3d& C0, CMtx3x3d& C1)
{
	VT_HR_BEGIN();

    // Find the scalar s s.t. | C0 + s*C1 | = 0
    //  This is a cubic equation, which we fit by sampling the determinant at 4 points
    CMtx4x4d A;
    CVec4d b;
    for (int i = 0; i < 4; i++)
    {
        int s = i-1;     // start @ -1
        CMtx3x3d Cb = C0 + C1 * (double) s;
        b[i] = Cb.Det();
        for (int j = 0, k = 1; j < 4; j++, k *= s)
        {
            A[i][j] = k;
        }
    }
    CVec4d a = A.Inv() * b;

    // Find the three solutions for s
    double roots[3];
    VT_HR_EXIT( CubicRoots(a, roots) );

    // For each of the roots, generate a pair of lines
    CVec3d lines[3][2];
    for (int ir = 0; ir < 3; ir++)
    {
        // Factor the rank-deficient matrix
        CMtx3x3d Cr = C0 + C1 * roots[ir];
        CVec3d d;
        CMtx3x3d L;
        LDLTPermuted(Cr, d, L, true, true);

        // Compute the line equation before L transform
        CVec3d m0(1,  1, 0);
        lines[ir][0] = L * m0;
        CVec3d m1(1, -1, 0);
        lines[ir][1] = L * m1;
    }

    // Compute the four intersection points from the first two pairs of lines
    for (int i0 = 0, k = 0; i0 < 2; i0++)
    {
        for (int i1 = 0; i1 < 2; i1++, k++)
        {
            points[k] = lines[0][i0].Cross(lines[1][i1]);
            double z = points[k][2];
            points[k] /= z;             // we don't expect ideal points (@ infinity), so easier for debugging

#ifdef _DEBUG
            // Check how well the points fit the conics
            double d0 = points[k] * (C0 * points[k]);
            double d1 = points[k] * (C1 * points[k]);
            double d2 = d0 + d1;
            VT_ASSERT(d2 == d2);       // remove compiler warning
#endif
        }
    }

    VT_HR_END();
}

// Find the *single* intersection point for this tangent line
CVec3d IntersectConicTangentLine(CMtx3x3d& C, CMtx3x3d& Li, CVec3d& l)
{
    // Transform the line by the "sqrt" matrix L, so that diagonal is (1,1,-1) (ellipse => unit circle)
    CVec3d m = Li * l;      // should be a scaled version of (n_x, n_y, 1)
    CVec3d y(m[0], m[1], -m[2]);
    CVec3d x = Li.T() * y;

    // Check that the point is on the line and conic
#ifdef _DEBUG
    double d1 = x * (C * x);
    double d2 = x * l;
    double d3 = hypot(d1, d2);
    VT_ASSERT(d3 >= 0);    // doesn't do anything, but just to remove compiler warnings
#else
    C[0][0] = C[0][0];  // doesn't do anything, but just to remove compiler warnings
#endif
    return x;
}

// Fit bi-tangents to the pair of dominant ellipses
HRESULT vt::VtFindBiTangentLines(OUT LineSegment bitangents[2], OUT CVec3f vanishingPoints[2],
    IN const vector<EllipseSegment>& ellipses)
{
	VT_HR_BEGIN();

	// If there are no ellipses, don't find bitangents
	if(ellipses.size() == 0)
		VT_HR_EXIT(E_INVALIDARG);	// need to have ellipses before computing bitangents

    // Compute the conic matrices for both ellipses
    CMtx3x3d C[2], Ci[2], L[2], Li[2];
    Conic(C[0], Ci[0], ellipses[0]);
    Conic(C[1], Ci[1], ellipses[1]);

    // Compute the intersections of the dual conics (tangent line bundles)
    CVec3d tangentLines[4];
    VT_HR_EXIT( IntersectConics(tangentLines, Ci[0], Ci[1]) );

    // Factor both conics and get ellipse centers
    CVec3d d;
    for (int i = 0; i < 2; i++)
    {
        LDLTPermuted(C[i], d, L[i], false, true);
        Li[i] = L[i].Inv();
    }
    CVec3d c0(ellipses[0].cx, ellipses[0].cy, 1);
    CVec3d c1(ellipses[1].cx, ellipses[1].cy, 1);

    // Test which tangent lines are outside the ellipses
    CVec3d p[2][2];
    for (int i = 0, k = 0; i < 4 && k < 2; i++)
    {
        CVec3d& l = tangentLines[i];
        double s0 = c0 * l;
        double s1 = c1 * l;
        if (s0 * s1 > 0)
        {
            // Compute the endpoints of each segment, i.e., ellipse/tangent intersections
            LineSegment& ls = bitangents[k];
            p[k][0] = IntersectConicTangentLine(C[0], Li[0], l);
            p[k][1] = IntersectConicTangentLine(C[1], Li[1], l);
            ls.start.x =  (float) p[k][0].Dehom().x;
            ls.start.y =  (float) p[k][0].Dehom().y;
            ls.end.x  =  (float) p[k][1].Dehom().x;
            ls.end.y  =  (float) p[k][1].Dehom().y;
            CVec2d v = ls.end - ls.start;
            ls.length = (float) v.Magnitude();
            k++;
        }
    }

    // Compute the two vanishing points
    for (int i = 0; i < 2; i++)
    {
        CVec3d l0 = (i) ? p[0][0].Cross(p[1][0]) : p[0][0].Cross(p[0][1]);
        CVec3d l1 = (i) ? p[0][1].Cross(p[1][1]) : p[1][0].Cross(p[1][1]);
        CVec3d vp = l0.Cross(l1);
        vp = (vp[2] < 0) ? -vp : vp;
        vanishingPoints[i] = CVec3f((float) vp[0], (float) vp[1], (float) vp[2]);
    }
    VT_HR_END();
}
