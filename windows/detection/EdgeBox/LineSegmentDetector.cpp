#include <memory>
//#include <gdiplus.h>

#define PI  3.14159265358979323846
#define PI_HALF 1.570796326794897

#include "LineSegmentDetector.h"

namespace WhiteboardCleanup
{
static const int c_maxImageSizeInPixels = INT_MAX;  // 512 1024
static const int c_angleBuffer = 45;            //60; //
static const float c_maxLineMergingDistancePerpendicular = 0.11f;
static const float c_maxLineMergingDistanceColinear = 0.3f;

HRESULT CLineSegmentDetector::ProjectOnLine(vt::CVec2d& proj, vt::CVecd& leq, vt::CVec2d& p)
{
    VT_HR_BEGIN()
        if (leq.Size() < 3)
            VT_HR_EXIT(E_INVALIDARG);
    double a, b, c, d, dn;
    a = leq[0];
    b = leq[1];
    c = leq[2];
    d = a*p[1] - b*p[0];
    dn = a*a + b*b;
    proj[0] = (-b*d - a*c) / dn;
    proj[1] = (a*d - b*c) / dn;

    VT_HR_END()
}

HRESULT CLineSegmentDetector::ProjectOnLine(vt::CVec2d& proj,
    vt::CVec3d& leq,
    vt::CVec2d& p)
{
    double a, b, c, d, dn;
    a = leq.x;
    b = leq.y;
    c = leq.z;
    d = a*p[1] - b*p[0];
    dn = a*a + b*b;
    proj[0] = (-b*d - a*c) / dn;
    proj[1] = (a*d - b*c) / dn;

    return S_OK;
}

int CLineSegmentDetector::PointBetweenLineSegments(vt::CVec2d& a,
    vt::CVec2d& lineStart,
    vt::CVec2d& lineEnd)
{
    double y1 = lineEnd.y - a.y;
    double y2 = lineStart.y - a.y;
    double x1 = lineEnd.x - a.x;
    double x2 = lineStart.x - a.x;

    int dy1 = fabs(y1) < 1.0 ? 0 : ((y1 > 0) ? 1 : -1);
    int dy2 = fabs(y2) < 1.0 ? 0 : ((y2 > 0) ? 1 : -1);
    int dx1 = fabs(x1) < 1.0 ? 0 : ((x1 > 0) ? 1 : -1);
    int dx2 = fabs(x2) < 1.0 ? 0 : ((x2 > 0) ? 1 : -1);

    int dy = dy1 * dy2;
    int dx = dx1 * dx2;

    if ((dx > 0 && dy >= 0) || (dx >= 0 && dy > 0))
        return false;	// -- outside
    else
        return true;	// -- inside or overlap
}

CLineSegmentDetector::CLineSegmentDetector(void)
{
	m_params.gradientThreshold = 20;
	m_params.anchorThreshold = 2;
	m_params.minLineLen = 15;
	m_params.lineFitErrThreshold = 1.8f;
    m_params.maxImageDims = c_maxImageSizeInPixels;
    m_params.lineSegmentMergingAngleThreshold = static_cast<float>(cos(m_params.maxImageDims * PI / 180.0));
	m_params.lineSegmentMergingDistThreshold  = 4.f;		// in pixels
}

void OffsetLineSegment(const vt::CVec2d& offset,
    LineSegmentEx& ls)
{
    ls.start += offset;
    ls.end += offset;
    ls.center += offset;

    // update line param
    ls.para = (ls.start.Hom()).Cross(ls.end.Hom());    
    //ls.para[2] = ls.start.x * ls.end.y - ls.end.x * ls.start.y;
}

HRESULT CLineSegmentDetector::MergeLineSegments(vt::vector<LineSegmentEx> &mergedLines, 
    vt::vector<LineSegmentEx> &lines, 
    int &nmerged, 
    const LineParams &ldParams)
{
    VT_HR_BEGIN()

        int linesSize = (int)lines.size();
    VT_HR_EXIT(mergedLines.reserve(linesSize));

    // Compute Disjoint Sets of vector<int>-Segments
    vt::vector<vt::vector<int>> sets;
    VT_HR_EXIT(ComputeDSets(lines, sets, ldParams.lineSegmentMergingAngleThreshold, ldParams.lineSegmentMergingDistThreshold));
    int nsets = (int)sets.size();

    // Re-create the list of merged line segments
    for (int j = 0; j < nsets; j++)
    {
        int ns = (int)sets[j].size();
        if (ns == 1)
        {
            LineSegmentEx lsc;
            int idx = sets[j][0];
            lsc.start = lines[idx].start;
            lsc.end = lines[idx].end;
            lsc.length = lines[idx].length;
            VT_HR_EXIT(mergedLines.push_back(lsc));
        }
        else if (ns > 1)
        {
            LineSegmentEx mergedSegment;
            vt::vector<double> thetas;
            double min = 2 * PI;

            for (int k = 0; k < ns; k++)
            {
                LineSegmentEx lc = lines[sets[j][k]];
                double theta = std::atan2(lines[sets[j][k]].end.y - lines[sets[j][k]].start.y, lines[sets[j][k]].end.x - lines[sets[j][k]].start.x);

                if (theta < 0.0f)
                {
                    theta += PI;
                }

                if (theta < min)
                {
                    min = theta;
                }

                thetas.push_back(theta);
            }

            if (min < PI / 2)
            {
                for (int k = 0; k < ns; k++)
                {
                    if (thetas[k] > PI * 3 / 4)
                    {
                        thetas[k] -= PI;
                    }
                }
            }

            double weightSum = 0;
            double angleSum = 0;

            for (int k = 0; k < ns; k++)
            {
                double lineWeight = lines[sets[j][k]].length;
                weightSum += lineWeight;
                double angle = thetas[k];

                angleSum += angle * lineWeight;
            }

            double averageAngle = angleSum / weightSum;

            double slope = std::tan(averageAngle);
            vt::CVec2d geometricCenter;
            geometricCenter.x = 0;
            geometricCenter.y = 0;

            for (int k = 0; k < ns; k++)
            {
                geometricCenter += lines[sets[j][k]].start + lines[sets[j][k]].end;
            }

            geometricCenter = geometricCenter / (ns * 2);

            vt::CVecd line_eq(3);
            line_eq[0] = slope;
            line_eq[1] = -1;
            line_eq[2] = geometricCenter.y - slope * geometricCenter.x;


            vt::vector<vt::CVec2d> points;
            VT_HR_EXIT(points.reserve(2 * ns));

            vt::CVec2d center;
            center.x = 0;
            center.y = 0;

            for (int k = 0; k < ns; k++)
            {
                vt::CVec2d projp1, projp2;
                LineSegmentEx lc = lines[sets[j][k]];

                ProjectOnLine(projp1, line_eq, lc.start);
                VT_HR_EXIT(points.push_back(projp1));
                center += projp1;

                ProjectOnLine(projp2, line_eq, lc.end);
                VT_HR_EXIT(points.push_back(projp2));
                center += projp2;
            }

            center *= 1.0 / (2.0*ns);
            vt::CVec2d ls, le;
            vt::CVec2d lv = (points[0] - center);
            double maxdp = 0.0, mindp = 0.0;

            for (int k = 0; k < (int)points.size(); k++)
            {
                vt::CVec2d vv = points[k] - center;
                double dotp = vv*lv;
                if (dotp > 0)
                {
                    if (fabs(dotp) > maxdp)
                    {
                        maxdp = fabs(dotp);
                        le = points[k];
                    }
                }
                else
                {
                    if (fabs(dotp) > mindp)
                    {
                        mindp = fabs(dotp);
                        ls = points[k];
                    }
                }
            }

            mergedSegment.start = ls;
            mergedSegment.end = le;
            mergedSegment.length = (float)(ls - le).Magnitude();
            VT_HR_EXIT(mergedLines.push_back(mergedSegment));
        }
    }

    nmerged = (int)mergedLines.size();

    VT_HR_END()
}

HRESULT CLineSegmentDetector::ComputeDSets(vt::vector<LineSegmentEx> &data, 
    vt::vector<vt::vector<int>> &sets, 
    float angleThreshold, 
    float distThreshold)
{
    VT_HR_BEGIN()

        int sz = (int)data.size();

    // Calculate the angle of each line segment, so we can quickly compare the angle
    // of lines against each other.

    vt::vector<double> segmentThetas;
    VT_HR_EXIT(segmentThetas.reserve(sz));

    for (int i = 0; i < sz; ++i)
    {
        double theta = std::atan2(data[i].end.y - data[i].start.y, data[i].end.x - data[i].start.x);
        if (theta < 0.0f)
        {
            theta += PI;
        }

        VT_HR_EXIT(segmentThetas.push_back(theta));
    }

    // Calculate the threshold, expressed as the raw angle.
    double angleThresholdTheta = std::acos(angleThreshold);

    // Compare each LineSegment against every other one, to separate them into disjoint sets of
    // close collinear lines. To conserve memory, ensure each comparison is done only once and
    // store the results in a way that is easy to pack eight booleans into a single byte. This
    // is accomplished by looking at the next size/2 lines for each line. This generates a grid
    // like the following, which can be stored nicely in a sz * sz/2 array. Then, since the
    // rows are uniform, we can further reduce the array to sz * sz/2/8.

    /*   0 1 2 3 4 5 6 7 8
    * 0   X X X X
    * 1     X X X X
    * 2       X X X X
    * 3         X X X X
    * 4           X X X X
    * 5 X           X X X
    * 6 X X           X X
    * 7 X X X           X
    * 8 X X X X
    */

    // The number of comparisons needed is sz*(sz-1)/2. For simplicity, the outer loop is sz.
    // For odd numbers (sz-1)/2 doesn't have a fraction, but for even numbers it does. To keep
    // the array uniform, we'll use sz/2 for even numbers and skip the last comparison for the
    // latter half of the rows. For odd numbers, (sz-1)/2 = sz/2 due to integer arithmetic, so
    // we can use sz/2 for both. To ensure we always have bytes for numbers that are not
    // divisible by 8, we add 7 before the integer division.

    int rowSize = (sz / 2 + 7) / 8;

    vt::vector<unsigned char> areCollinearAndClose;
    if (sz * rowSize > 0)
    {
        VT_HR_EXIT(areCollinearAndClose.resize(sz * rowSize));
        memset(&areCollinearAndClose[0], 0, sz * rowSize * sizeof(unsigned char));
    }

#ifdef USE_PPL
    Concurrency::parallel_for(0, sz, [&](int i)
#else
#pragma omp parallel for
    for (int i = 0; i < sz; i++)
#endif
    {
        // As mentioned above, with even numbers the last comparison can be skipped for the
        // latter half of the rows.

        int skip = (sz % 2 == 0 && i >= sz / 2) ? 1 : 0;

        for (int k = 0; k < sz / 2 - skip; k++)
        {
            // Compare to the kth element, starting with the number after i. Wrap around to
            // the beginning of the row when the end is reached.

            int j = (k + i + 1) % sz;

            // Compare line angles to make sure they're close.

            if (fabs(segmentThetas[i] - segmentThetas[j]) < angleThresholdTheta
                || fabs(fabs(segmentThetas[i] - segmentThetas[j]) - PI) < angleThresholdTheta)
            {
                bool value;

                // We can't use VT_HR_EXIT here if we're running in parallel, so just do nothing on
                // failure. It just means we might not merge as many lines as we could.

                if (SUCCEEDED(LineSegmentsAreCollinearAndClose(value, angleThreshold, distThreshold, data[i], data[j])) && value)
                    //if (SUCCEEDED(LineSegmentsAreCollinearAndClose(value, angleThreshold, data[i], data[j])) && value)
                {
                    areCollinearAndClose[i*rowSize + k / 8] |= (1 << (k % 8));
                }
            }
        }
    }
#ifdef USE_PPL
    );
#endif

    // vt::DisjointSets isn't thread-safe, so we saved the results of
    // LineSegmentsAreCollinearAndClose into areCollinearAndClose. Use
    // the results here to build up the disjoint sets.

    DisjointSets ds;
    ds.Init(sz);
    for (int i = 0; i < sz; i++)
    {
        int skip = (sz % 2 == 0 && i >= sz / 2) ? 1 : 0;

        for (int k = 0; k < sz / 2 - skip; k++)
        {
            int j = (k + i + 1) % sz;

            if (areCollinearAndClose[i*rowSize + k / 8] & (1 << (k % 8)))
            {
                if (ds.FindSet(i) != ds.FindSet(j))
                {
                    ds.Union(i, j);
                }
            }
        }
    }

    vt::vector<int> setIds;
    ds.SetIds(setIds);
    int m = 0;
    for (int i = 0; i < (int)setIds.size(); i++)
    {
        m = vt::VtMax(m, setIds[i]);
    }

    vt::vector<vt::vector<int>> st;
    VT_HR_EXIT(st.resize(m + 1));
    VT_HR_EXIT(sets.reserve(sz));

    for (int i = 0; i < (int)setIds.size(); i++)
    {
        int setId = setIds[i];
        VT_HR_EXIT(st[setId].push_back(i));
    }

    for (int i = 0; i < (int)st.size(); i++)
    {
        if (st[i].size() > 0)
        {
            VT_HR_EXIT(sets.push_back(st[i]));
        }
    }

    VT_HR_END()
}

HRESULT CLineSegmentDetector::DetectEdges(const vt::CByteImg& imgSrc,
    bool smoothInput,
    std::vector<POINT>& edgePts)
{
    VT_HR_BEGIN()

    //// scaling
    //const int iWidth = imgSrc.Width();
    //const int iHeight = imgSrc.Height();
    //m_flLineScalar = (float)__max(iWidth, iHeight) / m_params.maxImageDims;

    //vt::CImg imgScaled;

    //if (m_flLineScalar < 1.f)
    //{
    //    m_flLineScalar = 1.f;

    //    imgSrc.Share(imgScaled);
    //}
    //else
    //{
    //    const int iDstWid = (int)(iWidth / m_flLineScalar + 0.5f);
    //    const int iDstHei = (int)(iHeight / m_flLineScalar + 0.5f);

    //    VT_HR_EXIT(VtResizeImage(imgScaled, vt::CRect(0, 0, iDstWid, iDstHei), imgSrc, vt::eSamplerKernel::eSamplerKernelBilinear));
    //}

    //PrepareInternalScaledImage(imgScaled);

    //vt::vector<vt::CByteImg> imgSrcChannels;
    //VT_HR_EXIT(SeperateRGBtoSingleChannel(imgSrcChannels, m_imgScaled));

    //DetectEdgesImpl(imgSrcChannels[0], smoothInput, edgePts);

    VT_HR_EXIT(DetectEdgesImpl(imgSrc, smoothInput, edgePts));

    VT_HR_END()
}

HRESULT CLineSegmentDetector::DetectLines(
    const vt::CImg& imgSrc,
    bool smoothInput,
    vt::vector<LineSegmentEx>& vecAllLines)
{
    // only supporting byte pixel type
    if (!imgSrc.IsValid()
        || EL_FORMAT(imgSrc.GetType()) != EL_FORMAT_BYTE)
    {
        return E_INVALIDARG;
    }

    VT_HR_BEGIN()

    // scaling
    const int iWidth = imgSrc.Width();
    const int iHeight = imgSrc.Height();
    m_flLineScalar = (float)__max(iWidth, iHeight) / m_params.maxImageDims;

    vt::CImg imgScaled;

    if (m_flLineScalar < 1.f)
    {
        m_flLineScalar = 1.f;

        imgSrc.Share(imgScaled);
    }
    else
    {
        const int iDstWid = (int)(iWidth / m_flLineScalar + 0.5f);
        const int iDstHei = (int)(iHeight / m_flLineScalar + 0.5f);

        VT_HR_EXIT(VtResizeImage(imgScaled, vt::CRect(0, 0, iDstWid, iDstHei), imgSrc, vt::eSamplerKernel::eSamplerKernelBilinear));
    }

    PrepareInternalScaledImage(imgScaled);

    // detection
    vecAllLines.clear();

    vt::vector<vt::CByteImg> imgSrcChannels;
    VT_HR_EXIT(SeperateRGBtoSingleChannel(imgSrcChannels, m_imgScaled));

    // edge detection in entire image
    if (imgSrc.Bands() >= 3)
    {
        vt::vector<vt::vector<LineSegmentEx>> rgbLines;
        rgbLines.resize(3);    // r, g, b

#ifdef USE_PPL
        Concurrency::parallel_for(0, 3, [&](int i)
#else
        for (int i = 0; i < 3; ++i)
#endif    
        {
            CLineSegmentDetector lsd;

            lsd.LineDetection(imgSrcChannels[i], smoothInput, rgbLines[i]);
        }
#ifdef USE_PPL
        );
#endif

        vecAllLines.reserve(rgbLines[0].size() + rgbLines[1].size() + rgbLines[2].size());
        for (int i = 0; i < 3; ++i)
        {
            for (size_t j = 0; j < rgbLines[i].size(); ++j)
            {
                vecAllLines.push_back(rgbLines[i][j]);
            }
        }
    }
    else
    {
        VT_ASSERT(imgSrc.Bands() == 1);

        LineDetection(imgSrcChannels[0], smoothInput, vecAllLines);
    }

	VT_HR_END()
}

HRESULT CLineSegmentDetector::PrepareInternalScaledImage(const vt::CImg& scaledImg)
{
    VT_ASSERT(EL_FORMAT(scaledImg.GetType()) == EL_FORMAT_BYTE);

    HRESULT hr = S_OK;

    if (scaledImg.GetType() == m_imgScaled.GetType())
    {
        hr = scaledImg.Share(m_imgScaled);
        if (FAILED(hr)) return hr;
    }
    else
    {
        VT_ASSERT(scaledImg.Bands() == 3 || scaledImg.Bands() == 1);

        const int width = scaledImg.Width();
        const int height = scaledImg.Height();

        hr = m_imgScaled.Create(width, height);
        if (FAILED(hr)) return hr;

        if (scaledImg.Bands() == 3)
        {
            for (int y = 0; y < height; ++y)
            {
                vt::Byte* ptrDst = m_imgScaled.BytePtr(y);
                const vt::Byte* ptrSrc = scaledImg.BytePtr(y);

                for (int x = 0; x < width; ++x)
                {
                    memcpy(ptrDst, ptrSrc, 3 * sizeof(vt::Byte));

                    ptrDst += 4;
                    ptrSrc += 3;
                }   // x
            }   // y
        }
        else
        {   
            for (int y = 0; y < height; ++y)
            {
                vt::Byte* ptrDst = m_imgScaled.BytePtr(y);
                const vt::Byte* ptrSrc = scaledImg.BytePtr(y);

                for (int x = 0; x < width; ++x)
                {
                    memset(ptrDst, *ptrSrc, 3);

                    ptrDst += 4;
                    ptrSrc++;
                }   // x
            }   // y
        }
    }

    return S_OK;
}

HRESULT CLineSegmentDetector::DetectEdgesImpl(const vt::CByteImg& imgSrc,
    bool smooth,
    std::vector<POINT>& edgePts)
{
    VT_HR_BEGIN()

    //step 1: smooth image
    vt::CByteImg imgSrcFilt;
    if (smooth)
        vt::VtGaussianSmooth(imgSrcFilt, imgSrc, 2.f);   //1.5f //1.0f); // 0.667f
    else
        imgSrc.Share(imgSrcFilt);

    vt::CIntImg imgDx, imgDy;

#ifdef USE_PPL
    Concurrency::parallel_for(0, 2, [&](int i)
    {
        if (i == 0)
        {
            SobelDetector(imgSrcFilt, imgDx, true);
        }
        else
        {
            SobelDetector(imgSrcFilt, imgDy, false);
        }
    });
#else
    VT_HR_EXIT(SobelDetector(imgSrcFilt, imgDx, true));
    VT_HR_EXIT(SobelDetector(imgSrcFilt, imgDy, false));
#endif

    vt::CIntImg imgOrgGrad, imgGrad, imgDirect;
    VT_HR_EXIT(ComputeGradient_Direction(imgDx, imgDy, imgOrgGrad, imgGrad, imgDirect));

    vt::vector<CVec2i> vecAnchor;
    VT_HR_EXIT(DetectAnchors(imgGrad, imgDirect, vecAnchor));

    EdgeChains edges;
    VT_HR_EXIT(LinkEdgesBetweenAnchor(imgGrad, imgDirect, vecAnchor, edges));

    //printf("%d\n", edges.pts.size());

    edgePts.resize(edges.pts.size());

    for (int i = 0; i < edges.pts.size(); ++i)
    {
        edgePts[i].x = edges.pts[i].x;
        edgePts[i].y = edges.pts[i].y;
    }

    VT_HR_END()
}

HRESULT CLineSegmentDetector::LineDetection(
    const vt::CByteImg& imgSrc,
	bool smooth,
    vt::vector<LineSegmentEx>& vecAllLines
	)
{
	VT_HR_BEGIN()

	//step 1: smooth image
    vt::CByteImg imgSrcFilt;
	if(smooth)
		vt::VtGaussianSmooth(imgSrcFilt, imgSrc, 0.667f);//1.0f); // 0.667f
	else
		imgSrc.Share(imgSrcFilt);

    vt::CIntImg imgDx, imgDy;

#ifdef USE_PPL
    Concurrency::parallel_for(0, 2, [&](int i)
    {
        if (i == 0)
        {
            SobelDetector(imgSrcFilt, imgDx, true);
        }
        else
        {
            SobelDetector(imgSrcFilt, imgDy, false);
        }
    });
#else
	VT_HR_EXIT( SobelDetector(imgSrcFilt, imgDx, true) );
	VT_HR_EXIT( SobelDetector(imgSrcFilt, imgDy, false) );
#endif

    vt::CIntImg imgOrgGrad, imgGrad, imgDirect;
	VT_HR_EXIT( ComputeGradient_Direction(imgDx, imgDy, imgOrgGrad, imgGrad, imgDirect) );

    vt::vector<CVec2i> vecAnchor;
	VT_HR_EXIT( DetectAnchors(imgGrad, imgDirect, vecAnchor) );

	EdgeChains edges;
	VT_HR_EXIT( LinkEdgesBetweenAnchor(imgGrad, imgDirect, vecAnchor, edges) );

    printf("%d\n", edges.pts.size());

    if (edges.numOfEdges > 0 && !edges.pts.empty())
    {
        //printf("edges: %d\n", edges.numOfEdges);
        //printf("pixels: %d\n", edges.pts.size());

        LineChains lines;
        VT_HR_EXIT(DectectLinesFromEdges(imgDirect, imgDx, imgDy, imgOrgGrad, edges, lines));

        //printf("line detected: %d\n", lines.numOfLines);

        const int iMaxDimension = __max(imgSrc.Width(), imgSrc.Height());
        const float fLineMinLengthRatio = 0.05f;
        const float lenThres = iMaxDimension * fLineMinLengthRatio;
        RecognizeKeyLines(lines, lenThres, vecAllLines);
    }

	VT_HR_END()
}

HRESULT CLineSegmentDetector::MergeCoLines(
    vt::vector<LineSegmentEx>& vecAllLines,
    vt::vector<LineSegmentEx>& vecMergedLines
	)
{
	VT_HR_BEGIN()

	// Merge lines so if something looks like ----- ----
	// we can merge into                      ----------
	int nmergedOld = (int)vecAllLines.size();
	int nmerged;

	float flMaxImageDim = (float) vt::VtMax(m_imgScaled.Height(), m_imgScaled.Width());
	m_params.lineSegmentMergingDistThreshold = flMaxImageDim * 0.006f;

	VT_HR_EXIT( MergeLineSegments(vecMergedLines, vecAllLines, nmerged, m_params) );
	
	while (nmerged != nmergedOld)
	{
		nmergedOld = nmerged;
		vt::vector<LineSegmentEx> mergedLinesTemp(vecMergedLines);
		vecMergedLines.clear();
		VT_HR_EXIT( MergeLineSegments(vecMergedLines, mergedLinesTemp, nmerged, m_params) );
	}

	const double fWidthLimited = (double) m_imgScaled.Width() - 1.0;
	const double fHeightLimited = (double) m_imgScaled.Height() - 1.0;
	for(size_t i = 0; i < vecMergedLines.size(); i ++)
	{
		LineSegmentEx& a = vecMergedLines[i];
		a.para = (a.start.Hom()).Cross(a.end.Hom());
		a.start.x = vt::VtClamp(a.start.x, 0.0, fWidthLimited);
		a.start.y = vt::VtClamp(a.start.y, 0.0, fHeightLimited);
		a.end.x = vt::VtClamp(a.end.x, 0.0, fWidthLimited);
		a.end.y = vt::VtClamp(a.end.y, 0.0, fHeightLimited);
		a.center.x = (a.start.x + a.end.x) * 0.5;
		a.center.y = (a.start.y + a.end.y) * 0.5;
		a.length = (float)(a.start - a.end).Magnitude();
	}

	VT_HR_END()
}

HRESULT CLineSegmentDetector::GetColorInformationOfAllLines(
    vt::vector<LineSegmentEx>& leftLines,
    vt::vector<LineSegmentEx>& rightLines,
    vt::vector<LineSegmentEx>& topLines,
    vt::vector<LineSegmentEx>& bottomLines
	)
{
	VT_HR_BEGIN()

	const int xLims = m_imgScaled.Width() - 1;
	const int yLims = m_imgScaled.Height() - 1;

	for(size_t i = 0; i < leftLines.size(); i ++)
	{
		ColorInformationOfVerticalLines(leftLines[i], xLims, yLims, true);
	}

	for(size_t i = 0; i < rightLines.size(); i ++)
	{
		ColorInformationOfVerticalLines(rightLines[i], xLims, yLims, false);
	}

	for(size_t i = 0; i < topLines.size(); i ++)
	{
		ColorInformationOfHorizontalLines(topLines[i], xLims, yLims, true);
	}

	for(size_t i = 0; i < bottomLines.size(); i ++)
	{
		ColorInformationOfHorizontalLines(bottomLines[i], xLims, yLims, false);
	}

	VT_HR_END()
}

HRESULT CLineSegmentDetector::GetColorInformationOfImageEdges(
    vt::vector<LineSegmentEx>& edges
	)
{
    VT_HR_BEGIN()

    VT_HR_EXIT(edges.resize(4));
	
	const int iWidth = m_imgScaled.Width();
	const int iHeight = m_imgScaled.Height();

	const int gapInit = 3;
	const int gap = gapInit + 5;
	const int gapStep = 1;

	vt::CVec3f mean0 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	vt::CVec3f mean1 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	int count0 = 0, count1 = 0;

	for(int i = 0; i < iHeight; i += 2)
	{
		const vt::RGBAPix* ptr0 = m_imgScaled.Ptr(i);
		for(int j = gapInit; j < gap; j += gapStep)
		{
			const vt::RGBAPix& pix0 = ptr0[j];
			mean0[0] += pix0.r;
			mean0[1] += pix0.g;
			mean0[2] += pix0.b;
			count0 ++;

			const vt::RGBAPix& pix1 = ptr0[iWidth - j];
			mean1[0] += pix1.r;
			mean1[1] += pix1.g;
			mean1[2] += pix1.b;
			count1 ++;
		}
	}

	float flCount0 = 1.0f / __max(count0, 1.0f);
	float flCount1 = 1.0f / __max(count1, 1.0f);
	for(int i = 0; i < 3; i ++)
	{
		mean0[i] *= flCount0;
		mean1[i] *= flCount1;
	}
	edges[0].inPixMean = mean0;
	edges[1].inPixMean = mean1;


	mean0 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	mean1 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	count0 = 0, count1 = 0;

	for(int j = gapInit; j < gap; j += gapStep)
	{
		const vt::RGBAPix* ptr0 = m_imgScaled.Ptr(j);
		const vt::RGBAPix* ptr1 = m_imgScaled.Ptr(iHeight - j);

		for(int i = 0; i < iWidth; i += 2)
		{
			const vt::RGBAPix& pix0 = ptr0[i];
			mean0[0] += pix0.r;
			mean0[1] += pix0.g;
			mean0[2] += pix0.b;
			count0 ++;

			const vt::RGBAPix& pix1 = ptr1[i];
			mean1[0] += pix1.r;
			mean1[1] += pix1.g;
			mean1[2] += pix1.b;
			count1 ++;
		}
	}

	flCount0 = 1.0f / __max(count0, 1.0f);
	flCount1 = 1.0f / __max(count1, 1.0f);
	for(int i = 0; i < 3; i ++)
	{
		mean0[i] *= flCount0;
		mean1[i] *= flCount1;
	}
	edges[2].inPixMean = mean0;	
	edges[3].inPixMean = mean1;	
	
	VT_HR_END()
}

void CLineSegmentDetector::ColorInformationOfVerticalLines(
	LineSegmentEx& line, 
	int xLims, int yLims, 
	bool side	// side = true: leftline, false: rightline
	) 
{
	vt::CVec2d p0, p1; // p0 : top point, p1 : bottom point
	if(line.start.y < line.end.y)
	{
		p0 = line.start; p1 = line.end;
	}
	else
	{
		p0 = line.end; p1 = line.start;
	}

	const float angle = (float)atan2(p1.y - p0.y, p1.x - p0.x);
	const float len = (float)(p1 - p0).Magnitude();
	
	const int gapInit = 3;
	const int gap = gapInit + 5;
	const int gapStep = 1;

	vt::CVec3f mean0 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	vt::CVec3f mean1 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	int count0 = 0, count1 = 0;

	for(int i = 0; i < len; i += 2)
	{
		vt::CVec2d p2;
		p2.x = p0.x + cos(angle) * (float) i;
		p2.y = p0.y + sin(angle) * (float) i;

		float theta = float(angle - PI_HALF);
		for(int j = gapInit; j < gap; j += gapStep)
		{
			vt::CVec2f q0;
			q0.x = (float)p2.x + cos(theta) * (float) j;
			q0.y = (float)p2.y + sin(theta) * (float) j;
			int x = (int) (q0.x + 0.5);
			int y = (int) (q0.y + 0.5);
			if(x < 0 || x > xLims || y < 0 || y > yLims)
				continue;

			const vt::RGBAPix& pix = m_imgScaled(x, y);
			mean0[0] += pix.r;
			mean0[1] += pix.g;
			mean0[2] += pix.b;
			count0 ++;
		}

		for(int j = gapInit; j < gap; j += gapStep)
		{
			vt::CVec2f q0;
			q0.x = (float)p2.x - cos(theta) * (float) j;
			q0.y = (float)p2.y - sin(theta) * (float) j;
			int x = (int) (q0.x + 0.5);
			int y = (int) (q0.y + 0.5);
			if(x < 0 || x > xLims || y < 0 || y > yLims)
				continue;

			const vt::RGBAPix& pix = m_imgScaled(x, y);
			mean1[0] += pix.r;
			mean1[1] += pix.g;
			mean1[2] += pix.b;
			count1 ++;
		}
	}

	float flCount0 = 1.0f / __max(count0, 1.0f);
	float flCount1 = 1.0f / __max(count1, 1.0f);
	for(int i = 0; i < 3; i ++)
	{
		mean0[i] *= flCount0;
		mean1[i] *= flCount1;
	}

	if(side)  // left line
	{
		line.inPixMean = mean0;		
		line.outPixMean = mean1;
		float diff = (mean0 - mean1).MagnitudeSq();
		line.pixContrast = __min(sqrt(diff / 3.0f) / 128.0f, 1.0f); 
	}
	else      // right line
	{
		line.inPixMean = mean1;		
		line.outPixMean = mean0;
		float diff = (mean0 - mean1).MagnitudeSq();
		line.pixContrast = __min(sqrt(diff / 3.0f) / 128.0f, 1.0f); 
	}
}

void CLineSegmentDetector::ColorInformationOfHorizontalLines(
	LineSegmentEx& line, 
	int xLims, int yLims, 
	bool side	// side = true: topline, false: bottomline
	) 
{
	vt::CVec2d p0, p1; // p0 : left point, p1 : right point
	if(line.start.x < line.end.x)
	{
		p0 = line.start; p1 = line.end;
	}
	else
	{
		p0 = line.end; p1 = line.start;
	}
	
	const float angle = (float)atan2(p1.y - p0.y, p1.x - p0.x);
	const float len = (float)(p1 - p0).Magnitude();
	
	const int gapInit = 3;
	const int gap = gapInit + 5;
	const int gapStep = 1;

	vt::CVec3f mean0 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	vt::CVec3f mean1 = vt::CVec3f(0.0f, 0.0f, 0.0f);
	int count0 = 0, count1 = 0;

	for(int i = 0; i < len; i += 2)
	{
		vt::CVec2f p2;
		p2.x = (float)p0.x + cos(angle) * (float) i;
		p2.y = (float)p0.y + sin(angle) * (float) i;

		float theta = float(angle - PI_HALF);
		for(int j = gapInit; j < gap; j += gapStep)
		{
			vt::CVec2f q0;
			q0.x = p2.x + cos(theta) * (float) j;
			q0.y = p2.y + sin(theta) * (float) j;
			int x = (int) (q0.x + 0.5);
			int y = (int) (q0.y + 0.5);
			if(x < 0 || x > xLims || y < 0 || y > yLims)
				continue;

			const vt::RGBAPix& pix = m_imgScaled(x, y);
			mean0[0] += pix.r;
			mean0[1] += pix.g;
			mean0[2] += pix.b;
			count0 ++;
		}

		for(int j = gapInit; j < gap; j += gapStep)
		{
			vt::CVec2f q0;
			q0.x = p2.x - cos(theta) * (float) j;
			q0.y = p2.y - sin(theta) * (float) j;
			int x = (int) (q0.x + 0.5);
			int y = (int) (q0.y + 0.5);
			if(x < 0 || x > xLims || y < 0 || y > yLims)
				continue;

			const vt::RGBAPix& pix = m_imgScaled(x, y);
			mean1[0] += pix.r;
			mean1[1] += pix.g;
			mean1[2] += pix.b;
			count1 ++;
		}
	}

	float flCount0 = 1.0f / __max(count0, 1.0f);
	float flCount1 = 1.0f / __max(count1, 1.0f);
	for(int i = 0; i < 3; i ++)
	{
		mean0[i] *= flCount0;
		mean1[i] *= flCount1;
	}

	if(side)  // top line
	{
		line.inPixMean = mean1;		
		line.outPixMean = mean0;
		float diff = (mean0 - mean1).MagnitudeSq();
		line.pixContrast = __min(sqrt(diff / 3.0f) / 128.0f, 1.0f); 
	}
	else      // bottom line
	{
		line.inPixMean = mean0;		
		line.outPixMean = mean1;	
		float diff = (mean0 - mean1).MagnitudeSq();
		line.pixContrast = __min(sqrt(diff / 3.0f) / 128.0f, 1.0f); 
	}
}


HRESULT CLineSegmentDetector::GetScaledImageSize(
	int& width, 
	int& height
	)
{
	VT_HR_BEGIN()
	
	width = m_imgScaled.Width();
	height = m_imgScaled.Height();

	VT_HR_END()
}

HRESULT CLineSegmentDetector::FilterLinesByTheta(
    const vt::vector<LineSegmentEx> &lines,
    vt::vector<LineSegmentEx> &leftLines,
    vt::vector<LineSegmentEx> &rightLines,
    vt::vector<LineSegmentEx> &topLines,
    vt::vector<LineSegmentEx> &bottomLines
	)
{
	VT_HR_BEGIN()

	vt::vector<float> edgelThetas;
	if (lines.size() > 0) VT_HR_EXIT(edgelThetas.resize(lines.size()));

	size_t numLeftLines = 0;
	size_t numTopLines = 0;
	size_t numRightLines = 0;
	size_t numBottomLines = 0;

	// estimate center
	vt::CVec2f center(m_imgScaled.Width() / 2.0f, m_imgScaled.Height() / 2.0f);

	for (size_t i = 0; i < lines.size(); ++i)
	{
		const LineSegmentEx &line = lines[i];
		float theta = (float) (VT_180_PI * atan2(line.end.y-line.start.y, line.end.x-line.start.x));
		if (theta < -180.0f)
		{
			theta += 360.0f;
		}
		if (theta >= 180.0f)
		{
			theta -= 360.0f;
		}

		if (theta >= -c_angleBuffer && theta <= c_angleBuffer)
		{
			if(line.center.y < center.y) 
				++numTopLines;
			else
				++numBottomLines;
		}
		if (theta >= (180 - c_angleBuffer) || theta <= (-180 + c_angleBuffer))
		{
			if(line.center.y < center.y) 
				++numTopLines;
			else
				++numBottomLines;
		}
		if (theta >= (90 - c_angleBuffer) && theta <= (90 + c_angleBuffer))
		{
			if(line.center.x < center.x) 
				++numLeftLines;
			else
				++numRightLines;
			
		}
		if (theta >= (-90 - c_angleBuffer) && theta <= (-90 + c_angleBuffer))
		{
			if(line.center.x < center.x) 
				++numLeftLines;
			else
				++numRightLines;		
		}

		edgelThetas[i] = theta;
	}

	if (numLeftLines > 0) VT_HR_EXIT(leftLines.reserve(numLeftLines));
	if (numTopLines > 0) VT_HR_EXIT(topLines.reserve(numTopLines));
	if (numRightLines > 0) VT_HR_EXIT(rightLines.reserve(numRightLines));
	if (numBottomLines > 0) VT_HR_EXIT(bottomLines.reserve(numBottomLines));

	for (size_t i = 0; i < lines.size(); ++i)
	{
		const LineSegmentEx &line = lines[i];

		float theta = edgelThetas[i];

		if (theta >= -c_angleBuffer && theta <= c_angleBuffer)
		{
			if(line.center.y < center.y) 
				VT_HR_EXIT(topLines.push_back(line));
			else
				VT_HR_EXIT(bottomLines.push_back(line));
		}
		if (theta >= (180 - c_angleBuffer) || theta <= (-180 + c_angleBuffer))
		{
			if(line.center.y < center.y) 
				VT_HR_EXIT(topLines.push_back(line));
			else
				VT_HR_EXIT(bottomLines.push_back(line));
		}
		if (theta >= (90 - c_angleBuffer) && theta <= (90 + c_angleBuffer))
		{
			if(line.center.x < center.x) 
				VT_HR_EXIT(leftLines.push_back(line));
			else
				VT_HR_EXIT(rightLines.push_back(line));		
		}
		if (theta >= (-90 - c_angleBuffer) && theta <= (-90 + c_angleBuffer))
		{
			if(line.center.x < center.x) 
				VT_HR_EXIT(leftLines.push_back(line));
			else
				VT_HR_EXIT(rightLines.push_back(line));	
		}
	}

	VT_HR_END()
}


HRESULT CLineSegmentDetector::SobelDetector(
    const vt::CByteImg& imgSrc,
    vt::CIntImg& imgGrad,
	bool bHorizontal
	)
{
	VT_HR_BEGIN()

	const int iWid = imgSrc.Width();
	const int iHei = imgSrc.Height();

	VT_HR_EXIT( imgGrad.Create(iWid, iHei) );
    memset(imgGrad.Ptr(0), 0, imgGrad.Width() * sizeof(vt::CIntImg::PixelType));
    memset(imgGrad.Ptr(iHei - 1), 0, imgGrad.Width() * sizeof(vt::CIntImg::PixelType));

    vt::CIntImg imgTemp;
	VT_HR_EXIT( imgTemp.Create(iWid, iHei) );

	if(bHorizontal)
	{
		for(int y = 0; y < iHei; y ++)
		{
			const BYTE* ptrS = imgSrc.BytePtr(y);
			int* ptrD = imgTemp.Ptr(y);
			for(int x = 1; x < iWid - 1; x ++)
			{
				ptrD[x] = (int)ptrS[x + 1] - ptrS[x - 1];
			}
		}
       
		for(int y = 1; y < iHei - 1; y ++)
		{
			int* ptrSPrev = imgTemp.Ptr(y - 1);
			int* ptrSCur = imgTemp.Ptr(y);
			int* ptrSNext = imgTemp.Ptr(y + 1);
			int* ptrD = imgGrad.Ptr(y);

            ptrD[0] = ptrD[iWid - 1] = 0;

			for(int x = 1; x < iWid - 1; x ++)
			{
				ptrD[x] = ptrSPrev[x] + 2 * ptrSCur[x] + ptrSNext[x];
			}
		}
	}
	else
	{
		for(int y = 1; y < iHei - 1; y ++)
		{
			const BYTE* ptrS1 = imgSrc.BytePtr(y - 1);
			const BYTE* ptrS2 = imgSrc.BytePtr(y + 1);
			int* ptrD = imgTemp.Ptr(y);
			for(int x = 0; x < iWid; x ++)
			{
				ptrD[x] = (int)ptrS2[x] - ptrS1[x];
			}
		}

		for(int y = 1; y < iHei - 1; y ++)
		{
			int* ptrS = imgTemp.Ptr(y);
			int* ptrD = imgGrad.Ptr(y);

            ptrD[0] = ptrD[iWid - 1] = 0;

			for(int x = 1; x < iWid - 1; x ++)
			{
				ptrD[x] = ptrS[x - 1] + 2 * ptrS[x] + ptrS[x + 1];
			}
		}
	}

	VT_HR_END()
}

HRESULT CLineSegmentDetector::ComputeGradient_Direction(
    const vt::CIntImg& imgDx,
    const vt::CIntImg& imgDy,
    vt::CIntImg& imgOrgGrad,
    vt::CIntImg& imgGrad,
    vt::CIntImg& imgDirect
	)
{
	VT_HR_BEGIN()

	const int iWid = imgDx.Width();
	const int iHei = imgDx.Height();

	VT_HR_EXIT( imgOrgGrad.Create(iWid, iHei) );
	VT_HR_EXIT( imgGrad.Create(iWid, iHei) );
	VT_HR_EXIT( imgDirect.Create(iWid, iHei) );

	const int thres = m_params.gradientThreshold + 1;

	for(int y = 0; y < iHei; y ++)
	{
		const int* ptrDx = imgDx.Ptr(y);
		const int* ptrDy = imgDy.Ptr(y);
		int* ptrOG = imgOrgGrad.Ptr(y);
		int* ptrG = imgGrad.Ptr(y);
		int* ptrD = imgDirect.Ptr(y);

		for(int x = 0; x < iWid; x ++)
		{
			int dx = abs(ptrDx[x]);
			int dy = abs(ptrDy[x]);
			int sumXY = dx + dy;

			ptrOG[x] = (sumXY >> 2);
			ptrG[x] = sumXY > thres ? ptrOG[x] : 0;
            //ptrG[x] = (~((sumXY - thres) >> 31) & ptrOG[x]);		

			ptrD[x] = dx < dy ? 255 : 0;
            //ptrD[x] = (((dx - dy) >> 31) & 255);
		}
	}

	VT_HR_END()
}

HRESULT CLineSegmentDetector::DetectAnchors(
    const vt::CIntImg& imgGrad,
    const vt::CIntImg& imgDirect,
    vt::vector<CVec2i>& vecAnchor
	)
{
	//vecAnchor.reserve(m_numEdgePixs * 2);

	VT_HR_BEGIN()

	const int iWid = imgGrad.Width();
	const int iHei = imgGrad.Height();
	const int iStep = 2; // scanIntervals_
	const int thres = m_params.anchorThreshold;

	for(int y = 1; y < iHei - 1; y += iStep)
	{
		const int* ptrDirc = imgDirect.Ptr(y);
		const int* ptrGradPrev = imgGrad.Ptr(y - 1);
		const int* ptrGradCur = imgGrad.Ptr(y);
		const int* ptrGradNext = imgGrad.Ptr(y + 1);

		for(int x = 1; x < iWid - 1; x += iStep)
		{
			if(ptrDirc[x] == 255) // horizontal: |dx|<|dy|;
			{
				if( ptrGradCur[x] >= (ptrGradPrev[x] + thres) && 
					ptrGradCur[x] >= (ptrGradNext[x] + thres) )
				{
					vecAnchor.push_back(CVec2i((UINT16)x, (UINT16)y));
				}
			}
			else // vertical: |dy|<=|dx|;
			{
				if( ptrGradCur[x] >= (ptrGradCur[x - 1] + thres) && 
					ptrGradCur[x] >= (ptrGradCur[x + 1] + thres) )
				{
					vecAnchor.push_back(CVec2i((UINT16)x, (UINT16)y));
				}
			}
		}
	}

	VT_HR_END()
}

void CLineSegmentDetector::TraceEdgePixels(const vt::CIntImg& imgGrad, 
    const vt::CIntImg& imgDirect, 
    vt::CByteImg& imgEdge, 
    int iWid, 
    int iHei,
    int& x, 
    int& y, 
    int& prev_x, 
    int& prev_y, 
    EdgeDirec& prevDirec, 
    EdgeDirec& nextDirec, 
    vt::vector<CVec2i>& vecPassEdgePnt, 
    int& iPtOffset)
{
	int iG1, iG2, iG3; //gValue1, gValue2, gValue3;

	int edgePntSize = (int) vecPassEdgePnt.size();

	while( (imgGrad(x, y) > 0) && (imgEdge(x, y) == 0) )
	{
		if(iPtOffset >= edgePntSize)
			break;

		imgEdge(x, y) = 1;        // Mark this pixel as an edge pixel

		vecPassEdgePnt[iPtOffset].x = (UINT16)x;
        vecPassEdgePnt[iPtOffset].y = (UINT16)y;
		iPtOffset ++;

		nextDirec = ED_UNKNOWN;

		if(imgDirect(x, y) == 255)
		{
			if(prevDirec == ED_UP || prevDirec == ED_DOWN)
			{
				if(x > prev_x) { nextDirec = ED_RIGHT; } else { nextDirec = ED_LEFT; } 	
			}
			prev_x = x;	prev_y = y;

			if(prevDirec == ED_RIGHT || nextDirec == ED_RIGHT)
			{
				if(x == iWid-1 || y == 0 || y == iHei-1)
					break;

				iG1 = imgGrad(x+1, y-1);
				iG2 = imgGrad(x+1, y);
				iG3 = imgGrad(x+1, y+1);

				if(iG1 >= iG2 && iG1 >= iG3)	//up-right
				{
					x = x + 1;  y = y - 1;
				}
				else if(iG3 >= iG2 && iG3 >= iG1)	//down-right
				{
					x = x + 1;	y = y + 1;

				}
				else //straight-right
				{
					x = x + 1;
				}
				prevDirec = ED_RIGHT;
			} 
			else if( (prevDirec == ED_LEFT) || (nextDirec == ED_LEFT) )
			{
				if(x == 0 || y == 0 || y == iHei - 1)
					break;

				iG1 = imgGrad(x-1, y-1);
				iG2 = imgGrad(x-1, y);
				iG3 = imgGrad(x-1, y+1);

				if(iG1 >= iG2 && iG1 >= iG3) //up-left
				{
					x = x - 1;	y = y - 1;

				}
                else if(iG3 >= iG2 && iG3 >= iG1) //down-left
				{
					x = x - 1;  y = y + 1;

				}
				else //straight-left
				{
					x = x - 1;
				}
				prevDirec = ED_LEFT;
			}
		}
		else
		{
			if(prevDirec == ED_RIGHT || prevDirec == ED_LEFT)
			{
				if(y > prev_y) { nextDirec = ED_DOWN; } else { nextDirec = ED_UP; }
			}
			prev_x = x;	 prev_y = y;

			if(prevDirec == ED_DOWN || nextDirec == ED_DOWN)
			{
				if(x == 0 || x == iWid-1 || y == iHei-1 )
					break;

				// Look at 3 neighbors to the down and pick the one with the max. gradient value
				iG1 = imgGrad(x+1, y+1); 
				iG2 = imgGrad(x, y+1);
				iG3 = imgGrad(x-1, y+1); 
				if(iG1 >= iG2 && iG1 >= iG3) //down-right
				{
					x = x+1;	y = y+1;	
				}
				else if(iG3 >= iG2 && iG3 >= iG1) //down-left
				{
					x = x-1;	y = y+1;
				}
				else //straight-down
				{
					y = y+1;
				}

				prevDirec = ED_DOWN;

			}
			else if(prevDirec == ED_UP || nextDirec == ED_UP)
			{
				if(x == 0 || x == iWid-1 || y == 0)
					break;

				// Look at 3 neighbors to the up and pick the one with the max. gradient value
				iG1 = imgGrad(x+1, y-1);
				iG2 = imgGrad(x, y-1);
				iG3 = imgGrad(x-1, y-1);
				if(iG1 >= iG2 && iG1 >= iG3) //up-right
				{	
					x = x+1;	y = y-1;
				}
				else if(iG3 >= iG2 && iG3 >= iG1) //up-left
				{
					x = x-1;	y = y-1;
				}
				else //straight-up
				{
					y = y-1;
				}
				prevDirec = ED_UP;
			}
		}
	}//end while go right
}

HRESULT CLineSegmentDetector::LinkEdgesBetweenAnchor(
    const vt::CIntImg& imgGrad,
    const vt::CIntImg& imgDirect,
    const vt::vector<CVec2i>& vecAnchor,
	EdgeChains& edges
	)
{
	VT_HR_BEGIN()

	const int iWid = imgDirect.Width();
	const int iHei = imgDirect.Height();
	const int numEdgePixs = (iWid * iHei) / 5;
	const int maxEdgeNum = numEdgePixs / 20;


    vt::CByteImg imgEdge;
	VT_HR_EXIT( imgEdge.Create(iWid, iHei) );
	imgEdge.Clear();

	int x, y, prev_x, prev_y; //lastX, lastY
	EdgeDirec prevDirec, nextDirec; //lastDirection, shouldGoDirection

    vt::vector<CVec2i> vecPass1EdgePnt, vecPass2EdgePnt;
	vecPass1EdgePnt.resize(numEdgePixs);	//pFirstPartEdgeX_, pFirstPartEdgeY_
	vecPass2EdgePnt.resize(numEdgePixs); //pSecondPartEdgeX_, pSecondPartEdgeY_
    vt::vector<UINT16> vecPass1EdgeS, vecPass2EdgeS;
	vecPass1EdgeS.resize(maxEdgeNum);	//pFirstPartEdgeS_
	vecPass2EdgeS.resize(maxEdgeNum);	//pSecondPartEdgeS_
 
	int iPtOffset1 = 0, iPtOffset2 = 0;  //offsetPFirst=0, offsetPSecond=0;
	int iPSOffset = 0; //offsetPS=0
	
	const int minLineLen = m_params.minLineLen + 1;

	for(int i = 0; i < (int)vecAnchor.size(); i ++)
	{
		x = vecAnchor[i].x;
		y = vecAnchor[i].y;

		if( imgEdge(x, y) )
			continue;
		
        vecPass1EdgeS[iPSOffset] = (UINT16)iPtOffset1;
		if(imgDirect(x, y) == 255) // the direction of this pixel is horizontal, then go right and left.
		{  
			prevDirec = ED_RIGHT;

			TraceEdgePixels(imgGrad, imgDirect, imgEdge, iWid, iHei, x, y, prev_x, prev_y, prevDirec, 
				nextDirec, vecPass1EdgePnt, iPtOffset1);

			prevDirec = ED_LEFT;

			x = vecAnchor[i].x;
			y = vecAnchor[i].y;
			imgEdge(x, y) = 0;
            vecPass2EdgeS[iPSOffset] = (UINT16)iPtOffset2;

			TraceEdgePixels(imgGrad, imgDirect, imgEdge, iWid, iHei, x, y, prev_x, prev_y, prevDirec, 
				nextDirec, vecPass2EdgePnt, iPtOffset2);
		}
		else // the direction of this pixel is vertical, go down and up
		{
			prevDirec = ED_DOWN;

			TraceEdgePixels(imgGrad, imgDirect, imgEdge, iWid, iHei, x, y, prev_x, prev_y, prevDirec, 
				nextDirec, vecPass1EdgePnt, iPtOffset1);

			prevDirec = ED_UP;

			x = vecAnchor[i].x;
			y = vecAnchor[i].y;
			imgEdge(x, y) = 0;
            vecPass2EdgeS[iPSOffset] = (UINT16)iPtOffset2;

			TraceEdgePixels(imgGrad, imgDirect, imgEdge, iWid, iHei, x, y, prev_x, prev_y, prevDirec, 
				nextDirec, vecPass2EdgePnt, iPtOffset2);
		}

		//only keep the edge chains whose length is larger than the minLineLen_;
		int edgeLen1 = iPtOffset1 - vecPass1EdgeS[iPSOffset];
		int edgeLen2 = iPtOffset2 - vecPass2EdgeS[iPSOffset];
		if(edgeLen1 + edgeLen2 < minLineLen) //short edge, drop it
		{
			iPtOffset1 = vecPass1EdgeS[iPSOffset];
			iPtOffset2 = vecPass2EdgeS[iPSOffset];
		}
		else
		{
			iPSOffset++;
		}

        if (iPSOffset >= maxEdgeNum)
        {
            --iPSOffset;
            break;
        }
	}

	//store the last index
    vecPass1EdgeS[iPSOffset] = (UINT16)iPtOffset1;
    vecPass2EdgeS[iPSOffset] = (UINT16)iPtOffset2;

	_ASSERT(iPSOffset < maxEdgeNum);
	_ASSERT(iPtOffset1 < numEdgePixs && iPtOffset2 < numEdgePixs);

	//we should reorganize all the edge information into edgeChains for easily using
	int tempID;
	edges.pts.resize(iPtOffset1 + iPtOffset2);
	edges.sId.resize(iPSOffset + 1);
	iPtOffset1 = 0;
	iPtOffset2= 0;
	int indexInCors = 0;
	int numOfEdges = 0;
	for(int edgeId = 0; edgeId < iPSOffset; edgeId++)
	{
		edges.sId[numOfEdges++] = indexInCors;
		
		int index = vecPass1EdgeS[edgeId];
		iPtOffset1 = vecPass1EdgeS[edgeId + 1];
		for(tempID = iPtOffset1 - 1; tempID >= index; tempID--)
		{
			edges.pts[indexInCors] = vecPass1EdgePnt[tempID];
			indexInCors ++;
		}

		index = vecPass2EdgeS[edgeId];
		iPtOffset2= vecPass2EdgeS[edgeId+1];
		for(tempID = index + 1; tempID < iPtOffset2; tempID++)
		{
			edges.pts[indexInCors] = vecPass2EdgePnt[tempID];
			indexInCors ++;
		}
	}

	for(int i = indexInCors; i < iPtOffset1 + iPtOffset2; i ++)
	{
		edges.pts[i].x = 0;
		edges.pts[i].y = 0;
	}

	edges.sId[numOfEdges] = indexInCors;//the end index of the last edge
	edges.numOfEdges = numOfEdges;

	VT_HR_END()
}

double CLineSegmentDetector::FitLineByLeastSquare(
    const vt::vector<CVec2i>& pts,
    const vt::CIntImg& imgDirect,
	int iOffsetS,
    vt::CVec2d& lineEqn
	)
{
	_ASSERT(m_params.minLineLen == (int) m_vecFitMatT.size());
	_ASSERT(m_params.minLineLen == (int) m_vecFitVec.size());

	const int iDirect = imgDirect(pts[iOffsetS].x, pts[iOffsetS].y);
	const int minLineLen = m_params.minLineLen;
	double fitError = 0;
	bool bDirect = false;

	// if the first point of line is horizontal:
	if(iDirect == 255)
	{
		for(int i = 0; i < minLineLen; i ++)
		{
			m_vecFitMatT[i].x = pts[iOffsetS + i].x;
			m_vecFitVec[i] = pts[iOffsetS + i].y;
		}
		bDirect = true;
	}

	if(iDirect == 0)
	{
		for(int i = 0; i < minLineLen; i ++)
		{
			m_vecFitMatT[i].x = pts[iOffsetS + i].y;
			m_vecFitVec[i] = pts[iOffsetS + i].x;
		}
		bDirect = true;
	}

	if(bDirect)
	{
		double Sx2 = 0.0, Sx1 = 0.0, Sxy = 0.0, Sy1 = 0.0;
		for(int i = 0; i < minLineLen; i ++)
		{
			double flX = m_vecFitMatT[i].x;
			Sx2 += (flX * flX);
			Sx1 += flX;

			double flY = m_vecFitVec[i];
			Sxy += (flX * flY);
			Sy1 += flY;
		}
		
		// m_mtxATA = fitMatT * {fitMatT}^T;
		m_mtxATA(0, 0) = Sx2;		m_mtxATA(0, 1) = Sx1;
		m_mtxATA(1, 0) = Sx1;		m_mtxATA(1, 1) = (double)minLineLen;

		// m_mtxATV = fitMatT * {fitVec}^T;
		m_mtxATV[0] = Sxy;		
        m_mtxATV[1] = Sy1;

		double coef = 1.0 / (m_mtxATA(0,0)*m_mtxATA(1,1) - m_mtxATA(0,1)*m_mtxATA(1,0));

		// lineEquation = svd.Invert(mtxATA) * matT * vec;
		lineEqn[0] = coef *( m_mtxATA(1,1)*m_mtxATV[0] - m_mtxATA(0,1)*m_mtxATV[1] );
		lineEqn[1] = coef *( m_mtxATA(0,0)*m_mtxATV[1] - m_mtxATA(1,0)*m_mtxATV[0] );
	}

	if(iDirect == 255)
	{
		// compute fitting error of line
		for(int i = 0; i < minLineLen; i++)
		{
			double flX = pts[iOffsetS + i].x;
			double flY = pts[iOffsetS + i].y;

			double coef = flY - flX * lineEqn[0] - lineEqn[1];
			fitError += coef * coef;
		}
		return sqrt(fitError);
	}

	if(iDirect == 0)
	{
		// compute fitting error of line
		for(int i = 0; i < minLineLen; i++)
		{
			double flX = pts[iOffsetS + i].x;
			double flY = pts[iOffsetS + i].y;

			double coef = flX - flY * lineEqn[0] - lineEqn[1];
			fitError += coef * coef;
		}
		return sqrt(fitError);
	}
	
	return 0;
}

double CLineSegmentDetector::FitLineByLeastSquare(
    const vt::vector<CVec2i>& pts,
    const vt::CIntImg& imgDirect,
	int iOffsetS, int iOffsetSNew,
	int iOffsetE,
    vt::CVec2d& lineEqn
	)
{
	_ASSERT(m_params.minLineLen == (int) m_vecFitMatT.size());
	_ASSERT(m_params.minLineLen == (int) m_vecFitVec.size());

	int iNewLen = iOffsetE - iOffsetSNew;

	_ASSERT(iNewLen > 0);

    vt::vector<vt::CVec2d> matT;
	matT.resize(iNewLen);
    vt::vector<double> vec;
	vec.resize(iNewLen);
	
	const int iDirect = imgDirect(pts[iOffsetS].x, pts[iOffsetS].y);
	
	// if the first point of line is horizontal:
	if(iDirect == 255)
	{
		for(int i = 0; i < iNewLen; i ++)
		{
			matT[i].x = pts[iOffsetSNew + i].x;
			matT[i].y = 1; 
			vec[i] = pts[iOffsetSNew + i].y;
		}
	}

	// if the first point of line is vertical:
	if(iDirect == 0)
	{
		for(int i = 0; i < iNewLen; i ++)
		{
			matT[i].x = pts[iOffsetSNew + i].y;
			matT[i].y = 1; 
			vec[i] = pts[iOffsetSNew + i].x;
		}
	}

	if(iDirect == 0 || iDirect == 255)
	{
		double Sx2 = 0.0, Sx1 = 0.0, Sxy = 0.0, Sy1 = 0.0;
		for(int i = 0; i < iNewLen; i ++)
		{
			double flX = matT[i].x;
			Sx2 += (flX * flX);
			Sx1 += flX;

			double flY = vec[i];
			Sxy += (flX * flY);
			Sy1 += flY;
		}

        vt::CMtx2x2d tmpMat;  // tmpMat = fitMatT * {fitMatT}^T;
		tmpMat(0, 0) = Sx2;		tmpMat(0, 1) = Sx1;
		tmpMat(1, 0) = Sx1;		tmpMat(1, 1) = (double)iNewLen;

        vt::CVec2d tmpVec;   // tmpVec = fitMatT * {fitVec}^T;
		tmpVec[0] = Sxy;		tmpVec[1] = Sy1;

		m_mtxATA = m_mtxATA + tmpMat;
		m_mtxATV = m_mtxATV + tmpVec;

		double coef = 1.0 / (m_mtxATA(0,0)*m_mtxATA(1,1) - m_mtxATA(0,1)*m_mtxATA(1,0));

		// lineEquation = svd.Invert(mtxATA) * matT * vec;
		lineEqn[0] = coef *( m_mtxATA(1,1)*m_mtxATV[0] - m_mtxATA(0,1)*m_mtxATV[1] );
		lineEqn[1] = coef *( m_mtxATA(0,0)*m_mtxATV[1] - m_mtxATA(1,0)*m_mtxATV[0] );
	}

	return 0;
}

bool CLineSegmentDetector::LineValidation(
    const vt::vector<CVec2i>& pts,
    const vt::CIntImg& imgDx,
    const vt::CIntImg& imgDy,
    int iOffsetS, int iOffsetE,
    vt::CVec3d& lineEqn,
    float& direct,
    double logNT)
{
    int iSize = iOffsetE - iOffsetS;
    vt::vector<float> ptsDirect;
    ptsDirect.resize(iSize);

    int meanGradX = 0, meanGradY = 0;
    for (int i = 0; i < iSize; i++)
    {
        int x = pts[iOffsetS + i].x;
        int y = pts[iOffsetS + i].y;
        int dx_ = imgDx(x, y);
        int dy_ = imgDy(x, y);
        meanGradX += dx_;
        meanGradY += dy_;

        ptsDirect[i] = atan2((float)(-dx_), (float)dy_);
    }

    float dx = fabs((float)lineEqn[1]);
    float dy = fabs((float)lineEqn[0]);
    if (meanGradX == 0 && meanGradY == 0) //not possible, if happens, it must be a wrong line,
    {
        return false;
    }
    if (meanGradX > 0 && meanGradY >= 0)  //first quadrant, and positive direction of X axis.
    {
        direct = atan2(-dy, dx);	//line direction is in fourth quadrant
    }
    if (meanGradX <= 0 && meanGradY > 0)  //second quadrant, and positive direction of Y axis.
    {
        direct = atan2(dy, dx);	//line direction is in first quadrant
    }
    if (meanGradX < 0 && meanGradY <= 0)	 //third quadrant, and negative direction of X axis.
    {
        direct = atan2(dy, -dx);	//line direction is in second quadrant
    }
    if (meanGradX >= 0 && meanGradY < 0)	//fourth quadrant, and negative direction of Y axis.
    {
        direct = atan2(-dy, -dx);	//line direction is in third quadrant
    }

    // then check whether the line is on the border of the image. We don't keep the border line.
    const int iWid = imgDx.Width();
    const int iHei = imgDx.Height();

    if (fabs(direct) < 0.15 || (PI - fabs(direct)) < 0.15) //Horizontal line
    {
        if (fabs(lineEqn[2]) < 0.1 || fabs(iHei - fabs(lineEqn[2])) < 0.1) //upper border or lower border
            return false;
    }

    if (fabs(fabs(direct) - PI*0.5) < 0.15) //Vertical line
    {
        if (fabs(lineEqn[2]) < 0.1 || fabs(iWid - fabs(lineEqn[2])) < 0.1) //left border or right border
            return false;
    }

    //if(fabs(direct) < 0.15|| (PI - fabs(direct)) < 0.15) //Horizontal line
    //{
    //	if(fabs(lineEqn[2]) < 10 || fabs(iHei - fabs(lineEqn[2])) < 10) //upper border or lower border
    //		return false;
    //}

    //if(fabs(fabs(direct) - PI*0.5) < 0.15) //Vertical line
    //{
    //	if(fabs(lineEqn[2]) < 10 || fabs(iWid - fabs(lineEqn[2])) < 10) //left border or right border
    //		return false;
    //}

    //count the aligned points on the line which have the same direction as the line.
    float disDirect;
    int k = 0;
    const float disThres = (float)PI / 8;
    for (int i = 0; i < iSize; i++)
    {
        disDirect = fabs(direct - ptsDirect[i]);
        if (fabs(2 * PI - disDirect) < disThres || disDirect < disThres) //same direction
            k++;
    }

    //now compute Number of False Alarms
    double ret = number_of_false_alarm(iSize, k, 0.125, logNT);

    return (ret > 0); //0 corresponds to 1 mean false alarm
}

HRESULT CLineSegmentDetector::DectectLinesFromEdges(
    const vt::CIntImg& imgDirect,
    const vt::CIntImg& imgDx,
    const vt::CIntImg& imgDy,
    const vt::CIntImg& imgOrgGrad,
	const EdgeChains& edges,
	LineChains& lines
	)
{
	VT_HR_BEGIN()

	int linePixelID = edges.sId[edges.numOfEdges];
	lines.linePts.resize(linePixelID);
	lines.sId.resize(5 * edges.numOfEdges);

	double logNT = 2.0 * ( log10((double)imgDirect.Width() * imgDirect.Height()) );
	
    int iStartOffset = 0;
    int iEndOffset = 0;
    int iNewOffsetS = 0; 

	double lineFitErr = 0;	//the line fit error;
    vt::CVec2d lineEqn;
    vt::CVec3d lineEqnUpdate;

	m_lineEqns.clear();
	m_lineEndPts.clear();
	m_lineDirect.clear();

	int numOfLines = 0;
	int offsetInLineArray = 0;
	float direction; //line direction
	const int minLineLen = m_params.minLineLen;
	const float lineFitErrThreshold = m_params.lineFitErrThreshold;

	m_vecFitMatT.resize(minLineLen);
	m_vecFitVec.resize(minLineLen);
	for(int i = 0; i < minLineLen; i ++)
		m_vecFitMatT[i].y = 1;

	for(int edgeID = 0; edgeID < edges.numOfEdges; edgeID++)
	{
		iStartOffset = edges.sId[edgeID];
		iEndOffset = edges.sId[edgeID + 1];

		while(iEndOffset > (iStartOffset + minLineLen)) //extract line segments from an edge, may find more than one segments
		{
			//find an initial line segment
			while(iEndOffset > (iStartOffset + minLineLen))
			{
				lineFitErr = FitLineByLeastSquare(edges.pts, imgDirect, iStartOffset, lineEqn);
				if(lineFitErr <= lineFitErrThreshold) //ok, an initial line segment detected
					break;
				iStartOffset += 2; //skip the first two pixel in the chain and try with the remaining pixels
			}

			if(lineFitErr > lineFitErrThreshold) //no line is detected
				break; 

			//An initial line segment is detected. Try to extend this line segment
			lines.sId[numOfLines] = offsetInLineArray;

			double coef1 = 0;	//for a line ax+by+c=0, coef1 = 1/sqrt(a^2+b^2);
			double pointToLineDis;	//for a line ax+by+c=0 and a point(xi, yi), pointToLineDis = coef1*|a*xi+b*yi+c|
			bool bExtended = true;
			bool bFirstTry = true;
			int numOfOutlier;	//to against noise, we accept a few outlier of a line.
			int tryTimes = 0;
			if( imgDirect(edges.pts[iStartOffset].x, edges.pts[iStartOffset].y) == 255 )  //y=ax+b, i.e. ax-y+b=0
			{
				while(bExtended)
				{
					tryTimes++;
					if(bFirstTry)
					{
						bFirstTry = false;
						for(int i = 0; i < minLineLen; i++) //First add the initial line segment to the line array
						{
							lines.linePts[offsetInLineArray++] = edges.pts[iStartOffset++];
						}
					}
					else	//after each try, line is extended, line equation should be re-estimated
					{
						//adjust the line equation
						lineFitErr = FitLineByLeastSquare(lines.linePts, imgDirect, lines.sId[numOfLines], iNewOffsetS, offsetInLineArray, lineEqn);
					}

					coef1 = 1 / sqrt(lineEqn[0] * lineEqn[0] + 1);
					numOfOutlier = 0;
					iNewOffsetS = offsetInLineArray;
					while(iEndOffset > iStartOffset)
					{
						pointToLineDis = fabs(lineEqn[0]*edges.pts[iStartOffset].x - edges.pts[iStartOffset].y + lineEqn[1])*coef1;
						lines.linePts[offsetInLineArray++] = edges.pts[iStartOffset++];

						if(pointToLineDis > lineFitErrThreshold)
						{
							numOfOutlier++;
							if(numOfOutlier > 3) break;
						}
						else  //we count number of connective outliers.
						{
							numOfOutlier = 0;
						}
					}
					//pop back the last few outliers from lines and return them to edge chain
					offsetInLineArray -= numOfOutlier;
					iStartOffset -= numOfOutlier;
					if( (offsetInLineArray - iNewOffsetS) <= 0 || tryTimes >= 6 ) //no new pixels are added.
					{
						bExtended = false;
					}
				}
				
				//the line equation coefficients,for line w1x+w2y+w3 =0, we normalize it to make w1^2+w2^2 = 1.
				lineEqnUpdate[0] = lineEqn[0]*coef1;
				lineEqnUpdate[1] = -coef1;
				lineEqnUpdate[2] = lineEqn[1] * coef1;
			}
			else
			{
				while(bExtended)
				{
					tryTimes++;
					if(bFirstTry)
					{
						bFirstTry = false;
						for(int i = 0; i < minLineLen; i++) //First add the initial line segment to the line array
						{
							lines.linePts[offsetInLineArray++] = edges.pts[iStartOffset++];
						}
					}
					else	//after each try, line is extended, line equation should be re-estimated
					{
						//adjust the line equation
						lineFitErr = FitLineByLeastSquare(lines.linePts, imgDirect, lines.sId[numOfLines], iNewOffsetS, offsetInLineArray, lineEqn);
					}

					coef1 = 1 / sqrt(lineEqn[0] * lineEqn[0] + 1);
					numOfOutlier = 0;
					iNewOffsetS = offsetInLineArray;
					while(iEndOffset > iStartOffset)
					{
						pointToLineDis = fabs(edges.pts[iStartOffset].x - lineEqn[0]*edges.pts[iStartOffset].y - lineEqn[1])*coef1;

						lines.linePts[offsetInLineArray++] = edges.pts[iStartOffset++];

						if(pointToLineDis > lineFitErrThreshold)
						{
							numOfOutlier++;
							if(numOfOutlier > 3) break;
						}
						else  //we count number of connective outliers.
						{
							numOfOutlier = 0;
						}
					}
					//pop back the last few outliers from lines and return them to edge chain
					offsetInLineArray -= numOfOutlier;
					iStartOffset -= numOfOutlier;
					if( (offsetInLineArray - iNewOffsetS) <= 0 || tryTimes >= 6 ) //no new pixels are added.
					{
						bExtended = false;
					}
				}

				//the line equation coefficients,for line w1x+w2y+w3 =0, we normalize it to make w1^2+w2^2 = 1.
				lineEqnUpdate[0] = coef1;
				lineEqnUpdate[1] = -lineEqn[0] * coef1;
				lineEqnUpdate[2] = -lineEqn[1] * coef1;
			}
			

			if(LineValidation(lines.linePts, imgDx, imgDy, lines.sId[numOfLines], offsetInLineArray, lineEqnUpdate, direction, logNT)) //check the line
			{
				//store the line equation coefficients
				m_lineEqns.push_back(lineEqnUpdate);
                vt::CVec4f lineEndPts;//line endpoints
				double a1 = lineEqnUpdate[1] * lineEqnUpdate[1];
				double a2 = lineEqnUpdate[0] * lineEqnUpdate[0];
				double a3 = lineEqnUpdate[0] * lineEqnUpdate[1];
				double a4 = lineEqnUpdate[2] * lineEqnUpdate[0];
				double a5 = lineEqnUpdate[2] * lineEqnUpdate[1];
				int Px = lines.linePts[ lines.sId[numOfLines] ].x;//first pixel
				int Py = lines.linePts[ lines.sId[numOfLines] ].y;
				lineEndPts[0] = float(a1 * Px - a3 * Py - a4);//x
				lineEndPts[1] = float(a2 * Py - a3 * Px - a5);//y
				Px = lines.linePts[offsetInLineArray - 1].x;//last pixel
				Py = lines.linePts[offsetInLineArray - 1].y;
				lineEndPts[2] = float(a1 * Px - a3 * Py - a4);//x
				lineEndPts[3] = float(a2 * Py - a3 * Px - a5);//y
				m_lineEndPts.push_back(lineEndPts);
				m_lineDirect.push_back(direction);
				numOfLines++;
			}
			else
			{
				offsetInLineArray = lines.sId[numOfLines];// line was not accepted, the offset is set back
			}
		}
	}  //end for(unsigned int edgeID=0; edgeID<edges.numOfEdges; edgeID++)
	
	lines.sId[numOfLines] = offsetInLineArray;
	lines.numOfLines = numOfLines;

	VT_HR_END()
}

void CLineSegmentDetector::RecognizeKeyLines(
	const LineChains& lines,
	float lengthThreshold,
    vt::vector<LineSegmentEx>& vecAllLines
	)
{
	const int numOfFinalLine = lines.numOfLines;

	float s1, e1, s2, e2, dx, dy;
	bool shouldChange;

	vecAllLines.clear();

	const float lenThres2 = lengthThreshold * lengthThreshold;

	for(int  i = 0; i < numOfFinalLine; i++)
	{
		float direction = m_lineDirect[i];
		shouldChange = false;
		s1 = m_lineEndPts[i][0]; //sx
		s2 = m_lineEndPts[i][1]; //sy
		e1 = m_lineEndPts[i][2]; //ex
		e2 = m_lineEndPts[i][3]; //ey
		dx = e1 - s1; //ex-sx
		dy = e2 - s2; //ey-sy
		if( (direction >= -0.75 * PI) && (direction < -0.25 * PI) )
		{
			if(dy > 0) { shouldChange = true; }
		}
		if( (direction >= -0.25 * PI) && (direction < 0.25 * PI) )
		{
			if(dx < 0) { shouldChange = true; }
		}
		if( (direction >= 0.25 * PI) && (direction < 0.75 * PI) )
		{
			if(dy < 0) { shouldChange = true; }
		}
		if( ((direction >= 0.75 * PI) && (direction < PI)) || 
			((direction >= -PI) && (direction < -0.75 * PI)) )
		{
			if(dx > 0) { shouldChange = true; }
		}

		float x1, y1, x2, y2;
		if(shouldChange)
		{
			x1 = e1;	y1 = e2;	x2 = s1;	y2 = s2;
		}
		else
		{
			x1 = s1;	y1 = s2;	x2 = e1;	y2 = e2;
		}

		float len = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);

		if(len > lenThres2)
		{
			LineSegmentEx l;
			l.start.x = x1;
			l.start.y = y1;
			l.end.x = x2;
			l.end.y = y2;
			l.length = sqrt(len);
			l.center.x = (x1 + x2) * 0.5;
			l.center.y = (y1 + y2) * 0.5;
			l.para[0] = y1 - y2;
			l.para[1] = -(x1 - x2);
			l.para[2] = x1 * y2 - x2 * y1;

			vecAllLines.push_back(l);
		}
	}
}

//+-----------------------------------------------------------------------
//  Convert image from RGB -> Gray
// 
//HRESULT CLineSegmentDetector::ConvertRGBtoGray(
//	CByteImg& imgDst, 
//	const CRGBAByteImg &imgSrc
//	)
//{
//	VT_HR_BEGIN()
//
//	const int iWid = imgSrc.Width();
//	const int iHei = imgSrc.Height();
//
//	VT_HR_EXIT( imgDst.Create(iWid, iHei) );
//
//	for(int y = 0; y < iHei; y ++)
//	{
//		const BYTE* ptrS = imgSrc.BytePtr(y);
//		BYTE* ptrD = imgDst.BytePtr(y);
//
//		for(int x = 0; x < iWid; x ++, ptrS+=4, ptrD++)
//		{
//			//*ptrD = ptrS[0];
//			//int iB = ptrS[0];
//			//int iG = ptrS[1];
//			//int iR = ptrS[2];
//
//			UINT32 iB = ptrS[0];
//			UINT32 iG = ptrS[1];
//			UINT32 iR = ptrS[2];
//
//			UINT32 iV = ((iR * 1225 + iG * 2404 + iB * 467 + 2048) >> 12);
//			*ptrD = (BYTE) iV;
//
//			//int iGray = iR;//(((117*iR + 601*iG + 306*iB + 512) >> 10));
//			//*ptrD = (BYTE) iGray;
//		}
//	}
//
//	VT_HR_END()
//}

HRESULT CLineSegmentDetector::SeperateRGBtoSingleChannel(
    vt::vector<vt::CByteImg>& imgChannels,
    const vt::CRGBAByteImg &imgSrc
	)
{
	VT_HR_BEGIN()

	const int iWid = imgSrc.Width();
	const int iHei = imgSrc.Height();

	VT_HR_EXIT( imgChannels.resize(3) );
	VT_HR_EXIT( imgChannels[0].Create(iWid, iHei) );
	VT_HR_EXIT( imgChannels[1].Create(iWid, iHei) );
	VT_HR_EXIT( imgChannels[2].Create(iWid, iHei) );

	for(int y = 0; y < iHei; y ++)
	{
		const BYTE* ptrS = imgSrc.BytePtr(y);
		BYTE* ptrR = imgChannels[0].BytePtr(y);
		BYTE* ptrG = imgChannels[1].BytePtr(y);
		BYTE* ptrB = imgChannels[2].BytePtr(y);

		for(int x = 0; x < iWid; x ++, ptrS+=4, ptrR++, ptrG++, ptrB++)
		{
			*ptrB = ptrS[0];
			*ptrG = ptrS[1];
			*ptrR = ptrS[2];
		}
	}

	VT_HR_END()
}

HRESULT CLineSegmentDetector::LineSegmentsAreCollinearAndClose(bool &value, float angleThreshold, float distThreshold, LineSegmentEx &a, LineSegmentEx &b)
{
    VT_HR_BEGIN()

    value = false;

    // estimate the angle between two line segments
    vt::CVec2d va = (a.start - a.end);
    vt::CVec2d vb = (b.start - b.end);

    a.length = (float)(a.start - a.end).Magnitude();
    b.length = (float)(b.start - b.end).Magnitude();

    double dotp = fabs((va*vb) / double(a.length*b.length));

    if (dotp > angleThreshold)
    {
        vt::CVec3d aeq = (a.start.Hom()).Cross(a.end.Hom());
        vt::CVec3d beq = (b.start.Hom()).Cross(b.end.Hom());
        vt::CVec2d asp, aep, bsp, bep;
        VT_HR_EXIT(ProjectOnLine(asp, beq, a.start));
        VT_HR_EXIT(ProjectOnLine(aep, beq, a.end));
        VT_HR_EXIT(ProjectOnLine(bsp, aeq, b.start));
        VT_HR_EXIT(ProjectOnLine(bep, aeq, b.end));

        // next, compute relationship between two line segments
        // type 0: ------				OR				------
        //					---------		 ---------
        // type 1: ------			OR		   ------
        //			   ---------		 ---------
        // type 2:		------		OR		---------
        //			   ---------			  ------

        int aspInb = PointBetweenLineSegments(asp, b.start, b.end);  // 0 - outside; 1 - inside
        int aepInb = PointBetweenLineSegments(aep, b.start, b.end);
        int bspIna = PointBetweenLineSegments(bsp, a.start, a.end);
        int bepIna = PointBetweenLineSegments(bep, a.start, a.end);

        int iRelationType = 2;
        int aSum = aspInb + aepInb;
        int bSum = bspIna + bepIna;

        if (aSum == 0 && bSum == 0)  { iRelationType = 0; }
        if (aSum == 1 && bSum == 1)  { iRelationType = 1; }
        if ((aSum == 2 && bSum <= 1) || (aSum <= 1 && bSum == 2))  { iRelationType = 2; }

        const double len = (double)vt::VtMax(a.length, b.length);
        double vd1 = (a.start - asp).MagnitudeSq();
        double vd2 = (a.end - aep).MagnitudeSq();
        double vd3 = (b.start - bsp).MagnitudeSq();
        double vd4 = (b.end - bep).MagnitudeSq();
        const double minPerpDist = 3;

        if (iRelationType != 0 && vt::VtMax(vd1, vd2) < minPerpDist &&
            vt::VtMax(vd3, vd4) < minPerpDist)
        {
            value = true;
        }
        else
        {
            if (iRelationType == 0 || iRelationType == 1)
            {
                double aMaxPerpendicular = Sqr(a.length*c_maxLineMergingDistancePerpendicular);
                double bMaxPerpendicular = Sqr(b.length*c_maxLineMergingDistancePerpendicular);

                if (vt::VtMax(vd1, vd2) < aMaxPerpendicular && vt::VtMax(vd3, vd4) < bMaxPerpendicular)
                {
                    const double vdmax = (a.length < b.length) ? vt::VtMax(vd1, vd2) : vt::VtMax(vd3, vd4);

                    if (iRelationType == 0)
                    {
                        double d1 = (b.start - asp).MagnitudeSq();
                        double d2 = (b.end - asp).MagnitudeSq();
                        double d3 = (b.start - aep).MagnitudeSq();
                        double d4 = (b.end - aep).MagnitudeSq();

                        double dmin = vt::VtMin(d1, vt::VtMin(d2, vt::VtMin(d3, d4)));

                        if (dmin < Sqr(vt::VtMax(len*c_maxLineMergingDistanceColinear - vdmax, 2.0)) &&
                            vdmax < Sqr(distThreshold))
                            value = true;
                    }
                    else
                    {
                        double dmin = vdmax;
                        if (dmin < Sqr(distThreshold*0.75))
                            value = true;
                    }
                }
            }

            if (iRelationType == 2)
            {
                double dmin = (a.length < b.length) ? vt::VtMax(vd1, vd2) : vt::VtMax(vd3, vd4);
                if (dmin < Sqr(vt::VtMin(len * c_maxLineMergingDistancePerpendicular, (double)distThreshold*0.75)))
                    value = true;

                /*double ds = (a.length < b.length) ? vd1 : vd3;
                double de = (a.length < b.length) ? vd2 : vd4;
                double dmax = vt::VtMax(ds, de);

                if(dmax < Sqr( vt::VtMin(len * c_maxLineMergingDistancePerpendicular, (double)minPerpThres * 0.5) ))
                value = true;*/
            }
        }


    }

    VT_HR_END()
}

HRESULT DisjointSets::Init(int n)
{
    VT_HR_BEGIN()
        VT_HR_EXIT(p.resize(n));
    VT_HR_EXIT(rank.resize(n));
    for (int i = 0; i<n; i++)
    {
        rank[i] = 0;
        p[i] = i;
    }
    VT_HR_END()
};

int DisjointSets::FindSet(int i)
{
    if (i != p[i])
        p[i] = FindSet(p[i]);
    return p[i];
}

void DisjointSets::Union(int i, int j)
{
    Link(FindSet(i), FindSet(j));
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
    for (int i = 0; i<(int)p.size(); i++)
        ids[i] = FindSet(i);
    VT_HR_END()
}

}