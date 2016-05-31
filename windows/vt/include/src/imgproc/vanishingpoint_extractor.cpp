//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routine for extracting curves and lines
//
//  History:
//      2011/11/25-kramnath
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"
#include "vt_vanishingpoint_extractor.h"

using namespace vt;

// Finds inlier lines given a vanishing point estimate and a cos threshold
HRESULT FindInliers(VanishingPoint& vanishingPoint, const vector<LineSegment>& lineSegments, 
    const vector<bool>& liveIndices, float cosThreshold);

// Finds inlier lines given a vanishing point estimate and a cos threshold
HRESULT FindInliers(VanishingPoint& vanishingPoint, const vector<LineSegment>& lineSegments, 
     float cosThreshold, float minLength);

// Estimates the vanishing point using least squares given a set of inlier indices
HRESULT ReEstimateVanishingPointUsingLeastSq(CVec3d& vp, const  vector<int>& inlierIndices, const vector<LineSegment>& lineSegments);

// checks to see if the lines are self-intersecting (T junctions)
bool SelfIntersecting(const CVec3d& p, const LineSegment& la);

//// checks to see if the lines are self-intersecting (T junctions)
bool SelfIntersecting(const CVec3d &vp, const LineSegment &la, const LineSegment &lb);

// used by find inliers and eval hypothesis
inline double FindSqCosAngleBetLines(const CVec2d& s, const CVec2d& e, const CVec3d& vp);

// computes the score (total line length) given vanishing point estimate and threshold
float EvalHypothesis( const CVec3d& vp, const vector<LineSegment>& lineSegments, 
    const vector<bool>& liveIndices, float cosThreshold);

// initializes the live indices list given a region of interest
HRESULT IntializelLiveIndicesUsingROI(vector<bool>& liveIndices, const vector<CVec3d>& lineEqs,
    const CRect* regionOfInterest, const vector<LineSegment>& lineSegments, float minLength);

// calls from all the main functions use this for performing the ransac estimation
HRESULT ExtractVanishingPointInternal(VanishingPoint& vanishingPoint, bool& success, vector<int>& validPositions, 
    vector<bool>& liveIndices, const vector<LineSegment>& lineSegments, const vector<CVec3d>& lineEquations, 
    const VPExtractorParams& vpParams, const CVec3d& horizonLine = CVec3d(0.0));

HRESULT 
    vt::VtExtractVanishingPoints(vector<VanishingPoint>& vanishingPoints, const vector<LineSegment>& lineSegments, const VPExtractorParams& vpParams)
{
    VT_HR_BEGIN()

    int nLines = (int) lineSegments.size();

    vector<int> validPos;
    vector<bool> liveIndices;
    vector<CVec3d> lineEqs;

    //VT_HR_EXIT( vanishingPoints.resize(0) );
    VT_HR_EXIT( vanishingPoints.resize(vpParams.numVPoints) );
    for(int i = 0; i < vpParams.numVPoints; ++i)
    {
        // initialze everything
        vanishingPoints[i].score = 0;
        vanishingPoints[i].vpt = 0;
        vanishingPoints[i].idx.clear();
    }

    VT_HR_EXIT( validPos.reserve(nLines) );
    VT_HR_EXIT( liveIndices.resize(nLines) );
    VT_HR_EXIT( lineEqs.resize(nLines) );

    // Pre-compute line equations for all line segments
    for (int i = 0; i < nLines; i++)
    {
        if(lineSegments[i].length > vpParams.minLength)
        {
            VT_HR_EXIT( validPos.push_back(i) );
            liveIndices[i] = true;
        }
        else
            liveIndices[i] = false;
        lineEqs[i] =  (lineSegments[i].start.Hom()).Cross(lineSegments[i].end.Hom());
    }

    // Extract vanishing points by running RANSAC
    for(int iter = 0; iter < vpParams.numVPoints; ++iter)
    {
        bool success;
        ExtractVanishingPointInternal(vanishingPoints[iter], success, validPos, liveIndices, lineSegments,
                                lineEqs, vpParams);
        if(!success)
            break;
    } // for

    // reset live lines
    for (int i = 0; i < nLines; i++)
        if(lineSegments[i].length > vpParams.minLength)
            liveIndices[i] = true;
        else
            liveIndices[i] = false;

    // re-compute inliers now with the entire set of lines
    for (int i = 0; i < vpParams.numVPoints; i++)
    {    
        vanishingPoints[i].idx.clear();
        if (vanishingPoints[i].score > 0)
        {
            VT_HR_EXIT( FindInliers( vanishingPoints[i], lineSegments, liveIndices, vpParams.fitLinesToVPCosAngleThreshold) );
        }
    } // for 

    for (int i = 0; i < vpParams.numVPoints; i++) 
    {
        if ((int) vanishingPoints[i].idx.size() < vpParams.inlierThreshold)
            vanishingPoints[i].score = 0;
    }

    VT_HR_END()
}


HRESULT
    vt::VtExtractVanishingPointUsingPrior(
    VanishingPoint& vanishingPoint, const CVec3d& priorVanishingPoint, const vector<LineSegment>& lineSegments, 
    const VPExtractorParams& vpParams, const CRect* roiRect)
{
    VT_HR_BEGIN()
    int nLines = (int) lineSegments.size();

    vector<int> validPos;
    vector<bool> liveIndices;
    vector<CVec3d> lineEqs;

    vanishingPoint.score = 0;
    vanishingPoint.vpt = priorVanishingPoint;
    vanishingPoint.idx.clear();

    VT_HR_EXIT( liveIndices.resize(nLines) );
    VT_HR_EXIT( lineEqs.resize(nLines) );

    for(int i = 0; i < nLines; ++i)
    {
         lineEqs[i] =  (lineSegments[i].start.Hom()).Cross(lineSegments[i].end.Hom());
    }

    if(roiRect != NULL)
    {
        VT_HR_EXIT( IntializelLiveIndicesUsingROI(liveIndices, lineEqs, roiRect, lineSegments, vpParams.minLength) );

        // Get inliers corresponding to the initial estimate
        VT_HR_EXIT( FindInliers(vanishingPoint, lineSegments, liveIndices, vpParams.fitLinesToPriorVPCosineAngleThreshold) );
    }
    else
    {
        // Get inliers corresponding to the initial estimate
        VT_HR_EXIT( FindInliers(vanishingPoint, lineSegments, vpParams.fitLinesToPriorVPCosineAngleThreshold, vpParams.minLength) );
    }

    // reset live indices (ZeroMemory?)
    ZeroMemory(&liveIndices[0], nLines*sizeof(bool));

    int numInl = (int) vanishingPoint.idx.size();
    VT_HR_EXIT( validPos.reserve(numInl) );
    for (int i = 0; i < (int) numInl; i++)
    {
        if(lineSegments[vanishingPoint.idx[i]].length > vpParams.minLength)
        {   
            liveIndices[vanishingPoint.idx[i]] = true;
            VT_HR_EXIT( validPos.push_back(vanishingPoint.idx[i]) );
        }
    }

    bool success;
    VT_HR_EXIT( ExtractVanishingPointInternal(vanishingPoint, success, validPos, liveIndices, lineSegments,
        lineEqs, vpParams) );

    if(success)
    {
        if ((int) vanishingPoint.idx.size() < vpParams.inlierThreshold)
            vanishingPoint.score = 0;
    }
    VT_HR_END()
}

HRESULT
    vt::VtExtractVanishingPointOnHorizon(
    VanishingPoint& vanishingPoint, const CVec3d& horizonLine, const vector<LineSegment>& lineSegments, 
    const VPExtractorParams& vpParams, const CRect* roiRect)
{
    VT_HR_BEGIN()

    int nLines = (int) lineSegments.size();

    vector<int> validPos;
    vector<bool> liveIndices;
    vector<CVec3d> lineEqs;

    // initialze everything
    vanishingPoint.score = 0;
    vanishingPoint.vpt = 0;
    vanishingPoint.idx.clear();

    VT_HR_EXIT( validPos.reserve(nLines) );
    VT_HR_EXIT( liveIndices.resize(nLines) );
    VT_HR_EXIT( lineEqs.resize(nLines) );

    if(roiRect != NULL)
    {
        for (int i = 0; i < nLines; i++)
        {
            lineEqs[i] =  (lineSegments[i].start.Hom()).Cross(lineSegments[i].end.Hom());
        }
        VT_HR_EXIT( IntializelLiveIndicesUsingROI(liveIndices, lineEqs, roiRect, lineSegments, vpParams.minLength) );
        for (int i = 0; i < nLines; i++)
        {
            if(liveIndices[i])
                VT_HR_EXIT( validPos.push_back(i) );
        }
    }
    else
    {
        // Pre-compute line equations for all line segments
        for (int i = 0; i < nLines; i++)
        {
            if(lineSegments[i].length > vpParams.minLength)
            {
                VT_HR_EXIT( validPos.push_back(i) );
                liveIndices[i] = true;
            }
            else
                liveIndices[i] = false;

            lineEqs[i] =  (lineSegments[i].start.Hom()).Cross(lineSegments[i].end.Hom());
        }
    }

    bool success;
    ExtractVanishingPointInternal(vanishingPoint, success, validPos, liveIndices, lineSegments,
        lineEqs, vpParams, horizonLine);

    if(success)
    {
        if ((int) vanishingPoint.idx.size() < vpParams.inlierThreshold)
            vanishingPoint.score = 0;
    }

    VT_HR_END()

}


HRESULT
    IntializelLiveIndicesUsingROI(vector<bool>& liveIndices, const vector<CVec3d>& lineEqs,
    const CRect* regionOfInterest, const vector<LineSegment>& lineSegments, float minLength)
{
     VT_HR_BEGIN()

    int nLines = (int) lineSegments.size();
    if(liveIndices.size() != lineSegments.size())
        VT_HR_EXIT( liveIndices.resize(nLines) );

    vector<CVec3d> roiLineEqs;
    VT_HR_EXIT( roiLineEqs.resize(4) );

    CVec2d topLeft = CVec2d(regionOfInterest->TopLeft().x, regionOfInterest->TopLeft().y);
    CVec2d topRight  = CVec2d(regionOfInterest->TopRight().x, regionOfInterest->TopRight().y);
    CVec2d bottomRight = CVec2d(regionOfInterest->BottomRight().x, regionOfInterest->BottomRight().y);
    CVec2d bottomLeft = CVec2d(regionOfInterest->BottomLeft().x, regionOfInterest->BottomLeft().y);

    roiLineEqs[0] = (topLeft.Hom()).Cross(topRight.Hom());
    roiLineEqs[1] = (topRight.Hom()).Cross(bottomRight.Hom());
    roiLineEqs[2] = (bottomRight.Hom()).Cross(bottomLeft.Hom());
    roiLineEqs[3] = (bottomLeft.Hom()).Cross(topLeft.Hom());

    // Pre-compute line equations for all line segments
    vector<CVec3d> intPts;
    intPts.resize(4);
    for (int i = 0; i < nLines; i++)
    {
        if(lineSegments[i].length > minLength)
        {
            //**** make this or 
            if(lineSegments[i].start.x > topLeft.x && lineSegments[i].start.y > topLeft.y && 
                lineSegments[i].end.x < bottomRight.x && lineSegments[i].end.y < bottomRight.y)
            {
                liveIndices[i] = true;
            }
            else
            {
                for(int k = 0; k < 4; ++k)
                {
                    intPts[k] = lineEqs[i].Cross(roiLineEqs[k]);
                }
                int numInl = 0;
                CVec2d pt = intPts[0].Dehom();
                if((lineSegments[i].start - pt) * (lineSegments[i].end - pt) < 0 && (topLeft - pt) * (topRight - pt) < 0)
                    numInl++;
                pt = intPts[1].Dehom();
                if((lineSegments[i].start - pt) * (lineSegments[i].end - pt) < 0 && (topRight - pt) * (bottomRight - pt) < 0)
                    numInl++;
                pt = intPts[2].Dehom();
                if((lineSegments[i].start - pt) * (lineSegments[i].end - pt) < 0 && (bottomRight - pt) * (bottomLeft - pt) < 0)
                    numInl++;
                pt = intPts[3].Dehom();
                if((lineSegments[i].start - pt) * (lineSegments[i].end - pt) < 0 && (bottomLeft - pt) * (topLeft - pt) < 0)
                    numInl++;

                if(numInl > 0 && numInl <= 2)
                    liveIndices[i] = true;
                else
                    liveIndices[i] = false;
            }
        }
        else
            liveIndices[i] = false;
    }
    VT_HR_END()
}


HRESULT 
    ExtractVanishingPointInternal(VanishingPoint& vanishingPoint, bool& success, vector<int>& validPositions, vector<bool>& liveIndices,
     const vector<LineSegment>& lineSegments, const vector<CVec3d>& lineEquations, const VPExtractorParams& vpParams, const CVec3d& horizonLine)
{
    VT_HR_BEGIN()

    success = false;
    int nValidPos = (int) validPositions.size();

	if(nValidPos == 0)
		VT_HR_EXIT(E_INVALIDARG);


    // cumulative sums
    vector<int> cumSum;
    VT_HR_EXIT( cumSum.resize(nValidPos) );
    // length of each segment
    vector<int> lineLength;
    VT_HR_EXIT( lineLength.resize(nValidPos) );
    // cdf that has the indices of each line
    // repeated length of line number of times
    vector<int> cdf;
    VT_HR_EXIT( cdf.reserve(nValidPos*5) ); // avg guess
    
    int numTries = VtMax(vpParams.minRansacTries, nValidPos);

    int totalSum = 0;
    for(int k = 0; k < nValidPos; ++k)
    {
        cumSum[k] = totalSum;
        float length = lineSegments[validPositions[k]].length;
        lineLength[k] = (int) (length+0.5f);
        totalSum += lineLength[k];

        for(int q = 0; q < lineLength[k]; ++q)
        {
            VT_HR_EXIT( cdf.push_back(k) );
        }
    }

    // if horizon line is not given
    if(horizonLine.x == 0.0 && horizonLine.y == 0.0 && horizonLine.z == 0.0)
    {
        for(int p = 0; p < numTries; ++p)
        {
            // Get a random index from the set
            int sample = (int) ( (rand()/(float) (RAND_MAX + 1)) * totalSum);

            // choose that sample
            int s1 = cdf[sample];

            // We remove the length of the currently selected line
            // so that later on we do not include indices that belong to a line
            // that is already selected. Becomes apparent in the later steps
            int reducedSamplingSpace = totalSum - lineLength[s1];

            // Select a sample from the reduced set of indices. If the index corresponds
            // to a line already selected then we increment it to select a random value
            sample = (int) ( (rand()/(float) (RAND_MAX + 1)) * reducedSamplingSpace);

            int s2 = s1;
            // if the sample selected is less than the 
            // index of the selected line, then we can just use the
            // selected random value as the one
            if(sample < cumSum[s1])
            {
                s2 = cdf[sample];
            }
            else
            {
                // if not, add to the sample, the length of the current line
                // this will give a random sampling from the next line not the 
                // current line
                if(sample + lineLength[s1] < (int) cdf.size())
                    s2 = cdf[sample + lineLength[s1]];
                
            }
            // Got the valid positions
            s1 = validPositions[s1];
            s2 = validPositions[s2];

            // They should not be the same.
            if(s1 == s2)
            {
                continue;
            }
            // two random lines are available at this point s1 and s2
            CVec3d vpHyp = lineEquations[s1].Cross(lineEquations[s2]);

            // If the two lines are not self intersecting then the hyp is good
            if(! SelfIntersecting(vpHyp, lineSegments[s1], lineSegments[s2]))	
            {
                float score = EvalHypothesis(vpHyp, lineSegments, liveIndices, vpParams.ransacAndLeastSqCosAngleThreshold);
                if (score > vanishingPoint.score) 
                {	// update best hypothesis
                    vanishingPoint.score = score;
                    vanishingPoint.vpt   = vpHyp;
                }
            }

        }// over all ransac iterations
    }
    else
    {
         for(int p = 0; p < numTries; ++p)
        {
            // Get a random index from the set
            int sample = (int) ( (rand()/(float) (RAND_MAX + 1)) * totalSum);

            // choose that sample
            int s = cdf[sample];
            // Got the valid positions
            s = validPositions[s];

            // two random lines are available at this point s1 and s2
            CVec3d vpHyp = lineEquations[s].Cross(horizonLine);
            float score = EvalHypothesis(vpHyp, lineSegments, liveIndices, vpParams.ransacAndLeastSqCosAngleThreshold);
            if (score > vanishingPoint.score) 
            {	// update best hypothesis
                vanishingPoint.score = score;
                vanishingPoint.vpt   = vpHyp;
            }

        }// over all ransac iterations
    }

    if(vanishingPoint.score > 0)
    {
        VT_HR_EXIT( FindInliers(vanishingPoint, lineSegments, liveIndices,  vpParams.ransacAndLeastSqCosAngleThreshold) );

        // least sq 
        VT_HR_EXIT( ReEstimateVanishingPointUsingLeastSq(vanishingPoint.vpt, vanishingPoint.idx,  lineSegments) );

        // Mark lines that are within a threshold of the previously found VP estimate
        // as used
        VT_HR_EXIT( FindInliers(vanishingPoint, lineSegments, liveIndices,  vpParams.fitLinesToVPCosAngleThreshold) );

        for (int i = 0; i < (int)vanishingPoint.idx.size(); i++)
            liveIndices[vanishingPoint.idx[i]] = false;

        int totLines = (int) lineSegments.size();
        validPositions.clear();
        VT_HR_EXIT( validPositions.reserve(totLines) );
        for (int i = 0; i < totLines; i++)
            if (liveIndices[i])
                VT_HR_EXIT( validPositions.push_back(i) );

        success = true;
    }

    VT_HR_END()
}


bool SelfIntersecting(const CVec3d &vp, const LineSegment &la, const LineSegment &lb)
{
    if (vp[2]==0)
        return false;

    CVec2d u, v, p2 = vp.Dehom();
    u = la.start - p2;
    v = la.end - p2;
    if (u*v < 0)
        return true;

    u = lb.start - p2;
    v = lb.end - p2;
    if (u*v < 0)
        return true;

    return false;
}

bool SelfIntersecting(const CVec3d& p, const LineSegment& la)
{
    if (p[2]==0)
        return false;
    CVec2d u, v, p2 = p.Dehom();
    u = la.start - p2;
    v = la.end - p2;
    if (u*v < 0)
        return true;
    return false;
}


HRESULT ReEstimateVanishingPointUsingLeastSq(CVec3d& vp, const  vector<int>& inlierIndices, const vector<LineSegment>& lineSegments)
{
    VT_HR_BEGIN()
        int numInliers = (int)inlierIndices.size();

    CMtxd M;
    VT_HR_EXIT( M.Create(numInliers,3) );

    for (int i = 0; i < numInliers; i++) 
    {
        const LineSegment& ls = lineSegments[inlierIndices[i]];
        CVec3d    le = ls.start.Hom().Cross(ls.end.Hom());
        double    dn = sqrt(le[0]*le[0] +le[1]*le[1]);
        double    wt = ls.length/dn;
        M(i,0) =  wt * le[0];
        M(i,1) =  wt * le[1];
        M(i,2) =  wt * le[2];
    } // for i

    if (M.Size() != 0)
    {
        CSolveSVDd  svdM(M);
        vp = svdM.GetBestNullSpaceVector();
    }

    VT_HR_END();
}

HRESULT FindInliers(VanishingPoint& vanishingPoint, const vector<LineSegment>& lineSegments, 
    const vector<bool>& liveIndices, float cosThreshold)
{
    VT_HR_BEGIN()
    vector<int>& indices = vanishingPoint.idx;
    CVec3d& vp = vanishingPoint.vpt;
    int numLines = (int) liveIndices.size();

    double threshold = (double)cosThreshold*(double)cosThreshold;

    VT_HR_EXIT( indices.reserve(numLines) );
    for (int i=0; i<numLines; i++) 
    {
        if (!liveIndices[i])
            continue;

        double cosSq = FindSqCosAngleBetLines(lineSegments[i].start, lineSegments[i].end, vp);

        if (cosSq > threshold)
        {
            if (!SelfIntersecting(vp, lineSegments[i]))
            {
                VT_HR_EXIT( indices.push_back(i) );
            }
        }
    } // for i

    VT_HR_END()
}

HRESULT FindInliers(VanishingPoint& vanishingPoint, const vector<LineSegment>& lineSegments, 
     float cosThreshold, float minLength)
{
    VT_HR_BEGIN()
    vector<int>& indices = vanishingPoint.idx;
    CVec3d& vp = vanishingPoint.vpt;
    int numLines = (int) lineSegments.size();

    double threshold = (double)cosThreshold*(double)cosThreshold;

    VT_HR_EXIT( indices.reserve(numLines) );
    for (int i=0; i<numLines; i++) 
    {
        if (lineSegments[i].length < minLength)
            continue;

        double cosSq = FindSqCosAngleBetLines(lineSegments[i].start, lineSegments[i].end, vp);

        if (cosSq > threshold)
        {
            if (!SelfIntersecting(vp, lineSegments[i]))
            {
                VT_HR_EXIT( indices.push_back(i) );
            }
        }
    } // for i

    VT_HR_END()
}

float EvalHypothesis( const CVec3d& vp, const vector<LineSegment>& lineSegments, 
    const vector<bool>& liveIndices, float cosThreshold)
{
    int    numLines = (int) liveIndices.size();
    double threshold = (double)cosThreshold*(double)cosThreshold;
    float score = 0.0;
    for (int i=0; i<numLines; i++) 
    {
        if (!liveIndices[i])
            continue;

        double cosSq = FindSqCosAngleBetLines(lineSegments[i].start, lineSegments[i].end, vp);

        if (cosSq > threshold)
        {
            if (!SelfIntersecting(vp, lineSegments[i]))
            {
                score+= lineSegments[i].length;
            }
        }
    } // for i
    return score;
}


double FindSqCosAngleBetLines(const CVec2d& s, const CVec2d& e, const CVec3d& vp)
{
    double lu_lu, lv_lv, lu_lv;
    CVec3d le;
    CVec2d midpt, lu, lv;
    midpt = 0.5*(s+e);
    le    = midpt.Hom().Cross(vp);
    lv[0] = -le[1];
    lv[1] =  le[0];
    lu    = s - e;
    lu_lv = lu*lv;
    lu_lu = lu*lu;
    lv_lv = lv*lv;
    // cos(theta)sq = (a.b)^2/(a.a * b.b)
    return (lu_lv*lu_lv)/(lu_lu * lv_lv); 
}