#pragma once

#include "vtcore.h"
#include <array>
#include <vector>

namespace WhiteboardCleanup
{

enum EdgeDirec {ED_UNKNOWN, ED_UP, ED_RIGHT, ED_DOWN, ED_LEFT};
typedef vt::CVec2<UINT16> CVec2i;

struct EdgeChains
{
	vt::vector<CVec2i> pts; // all coordinates of edge points
    vt::vector<int> sId;	 //the start index of each edge in the coordinate arrays
	int numOfEdges;			//the number of edges whose length are larger than minLineLen; numOfEdges < sId.size;
};

struct LineChains{
    vt::vector<CVec2i> linePts;	//all coordinates of line points
    vt::vector<int> sId;	//the start index of each line in the coordinate arrays
	int numOfLines;		//the number of lines whose length are larger than minLineLen; numOfLines < sId.size;
};

struct LineParams{
	int gradientThreshold;
	int anchorThreshold;
	int minLineLen;
	float lineFitErrThreshold;
	int maxImageDims;

	/// <summary> Distance Threshold (in pixels) used to decide whether two line segments can be merged.  </summary>
	float	lineSegmentMergingDistThreshold;  
	
    /// <summary> Angle Threshold used to decide if two line segments can be merged </summary>
	float	lineSegmentMergingAngleThreshold; 
};

struct LineSegmentEx{
    vt::CVec2d center;
    vt::CVec2d start;
    vt::CVec2d end;
	float length;
    vt::CVec3d para; // (a, b, c) of ax+by+c=0;

    vt::CVec3f outPixMean;	// color average of pixels outside of line segments
    vt::CVec3f inPixMean;	// color average of pixels inside of line segments
	float pixContrast;  // color contrast of inside and outside pixels
};

class CLineSegmentDetector
{
public:
	CLineSegmentDetector(void);

public: 	
    HRESULT DetectLines(const vt::CImg& imgSrc, 
        bool smoothInput, 
        vt::vector<LineSegmentEx>& vecAllLines);

    HRESULT DetectEdges(const vt::CByteImg& imgSrc,
        bool smoothInput,
        std::vector<POINT>& edgePts);

    HRESULT MergeCoLines(vt::vector<LineSegmentEx>& vecAllLines, vt::vector<LineSegmentEx>& vecMergedLines);

    HRESULT FilterLinesByTheta(const vt::vector<LineSegmentEx> &lines, 
        vt::vector<LineSegmentEx> &leftLines,
        vt::vector<LineSegmentEx> &rightLines, 
        vt::vector<LineSegmentEx> &topLines, 
        vt::vector<LineSegmentEx> &bottomLines);

	HRESULT GetScaledImageSize(int& width, int& height);

    HRESULT GetColorInformationOfAllLines(vt::vector<LineSegmentEx>& leftLines, 
        vt::vector<LineSegmentEx>& rightLines,
        vt::vector<LineSegmentEx>& topLines, 
        vt::vector<LineSegmentEx>& bottomLines);

    HRESULT GetColorInformationOfImageEdges(vt::vector<LineSegmentEx>& edges);

    HRESULT MergeLineSegments(vt::vector<LineSegmentEx> &mergedLines, 
        vt::vector<LineSegmentEx> &lines, 
        int &nmerged, 
        const LineParams &ldParams);

    HRESULT ComputeDSets(vt::vector<LineSegmentEx> &data, 
        vt::vector<vt::vector<int>> &sets, 
        float angleThreshold, 
        float distThreshold);

private:
    HRESULT DetectEdgesImpl(const vt::CByteImg& imgSrc,
        bool smooth,
        std::vector<POINT>& edgePts);

    HRESULT PrepareInternalScaledImage(const vt::CImg& imgScaled);

    HRESULT LineDetection(const vt::CByteImg& imgSrc, 
        bool smooth, 
        vt::vector<LineSegmentEx>& vecAllLines);

    HRESULT SobelDetector(const vt::CByteImg& imgSrc, 
        vt::CIntImg& imgGrad, 
        bool bHorizontal);

    HRESULT ComputeGradient_Direction(const vt::CIntImg& imgDx, 
        const vt::CIntImg& imgDy, 
        vt::CIntImg& imgOrgGrad, 
        vt::CIntImg& imgGrad, 
        vt::CIntImg& imgDirect);

    HRESULT DetectAnchors(const vt::CIntImg& imgGrad, 
        const vt::CIntImg& imgDirect, 
        vt::vector<CVec2i>& vecAnchor);

    HRESULT LinkEdgesBetweenAnchor(const vt::CIntImg& imgGrad, 
        const vt::CIntImg& imgDirect, 
        const vt::vector<CVec2i>& vecAnchor, 
        EdgeChains& edges);

    HRESULT DectectLinesFromEdges(const vt::CIntImg& imgDirect, 
        const vt::CIntImg& imgDx, 
        const vt::CIntImg& imgDy, 
        const vt::CIntImg& imgOrgGrad,
		const EdgeChains& edges, 
        LineChains& lines);

	void RecognizeKeyLines(const LineChains& lines, 
        float lengthThreshold, 
        vt::vector<LineSegmentEx>& vecAllLines);

    HRESULT ConvertRGBtoGray(vt::CByteImg& imgDst, 
        const vt::CRGBAByteImg& imgSrc);

    HRESULT SeperateRGBtoSingleChannel(vt::vector<vt::CByteImg>& imgChannels, 
        const vt::CRGBAByteImg& imgSrc);

    double FitLineByLeastSquare(const vt::vector<CVec2i>& pts, 
        const vt::CIntImg& imgDirect, 
        int iOffsetS, 
        vt::CVec2d& lineEqn);

    double FitLineByLeastSquare(const vt::vector<CVec2i>& pts, 
        const vt::CIntImg& imgDirect, 
        int iOffsetS, 
        int iOffsetSNew, 
        int iOffsetE, 
        vt::CVec2d& lineEqn);

    bool LineValidation(const vt::vector<CVec2i>& pts, 
        const vt::CIntImg& imgDx, 
        const vt::CIntImg& imgDy, 
        int iOffsetS, 
        int iOffsetE, 
        vt::CVec3d& lineEqn,
		float& direct, 
        double logNT);

	void ColorInformationOfVerticalLines(LineSegmentEx& line, 
        int xLims, 
        int yLims, 
        bool side);
	
    void ColorInformationOfHorizontalLines(LineSegmentEx& line, 
        int xLims, 
        int yLims, 
        bool side);

    inline void TraceEdgePixels(const vt::CIntImg& imgGrad, 
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
        int& iPtOffset);

    HRESULT LineSegmentsAreCollinearAndClose(bool &value, float angleThreshold, float distThreshold, LineSegmentEx &a, LineSegmentEx &b);

    HRESULT ProjectOnLine(vt::CVec2d &proj, vt::CVecd &leq, vt::CVec2d &p);

    HRESULT ProjectOnLine(vt::CVec2d &proj, vt::CVec3d &leq, vt::CVec2d &p);

    int PointBetweenLineSegments(vt::CVec2d& a,
        vt::CVec2d& lineStart,
        vt::CVec2d& lineEnd);

    inline double Sqr(double a)
    {
        return a * a;
    }

	inline double log_gamma(double x)
	{
		if(x > 15.0) //log_gamma_windschitl
		{
			double a = log( x*sinh(1/x) + 1/(810.0*pow(x,6.0)) );
			return (0.918938533204673 + (x-0.5)*log(x) - x + 0.5*x*a);
		}
		else //log_gamma_lanczos
		{
			static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477,
				8687.24529705, 1168.92649479, 83.8676043424, 2.50662827511 };
			double a = (x+0.5) * log(x+5.5) - (x+5.5);
			double b = 0.0;
			int n;
			for(n=0;n<7;n++){
				a -= log( x + (double) n );
				b += q[n] * pow( x, (double) n );
			}
			return a + log(b);
		}
	}

	inline bool almost_equal(double a, double b)
	{
		double abs_diff,aa,bb,abs_max;
		if( a == b ) return true;
		abs_diff = fabs(a-b);
		aa = fabs(a);
		bb = fabs(b);
		abs_max = aa > bb ? aa : bb;
		if( abs_max < DBL_MIN ) abs_max = DBL_MIN;
		return (abs_diff / abs_max) <= (100.0 * DBL_EPSILON);
	}

	inline double number_of_false_alarm(int n, int k, double p, double logNT)
	{
		double tolerance = 0.1;       // an error of 10% in the result is accepted
		double log1term,term,bin_term,mult_term,bin_tail,err,p_term;
		int i;

		if( n<0 || k<0 || k>n || p<=0.0 || p>=1.0 )
		{
			printf("nfa: wrong n, k or p values.\n");
			exit(0);
		}

		if( n==0 || k==0 ) return -logNT;
		if( n==k ) return -logNT - (double) n * log10(p);

		p_term = p / (1.0-p);

		log1term = log_gamma( (double)n+1.0 ) - log_gamma( (double)k+1.0 ) - log_gamma( (double)(n-k)+1.0 )
					+ (double)k * log(p) + (double)(n-k) * log(1.0-p);
		term = exp(log1term);

		const double LN10 = 2.30258509299404568402;

		if( almost_equal(term, 0.0) ) // the first term is almost zero
		{  
			if( (double) k > (double) n * p )    
				return -log1term / LN10 - logNT; 
			else
				return -logNT;                     
		}

		bin_tail = term;
		for(i=k+1;i<=n;i++)
		{
			bin_term = (double) (n-i+1) / (double) i;
			mult_term = bin_term * p_term;
			term *= mult_term;
			bin_tail += term;
			if(bin_term<1.0)
			{                           
				err = term * ( ( 1.0 - pow( mult_term, (double) (n-i+1) ) ) / (1.0-mult_term) - 1.0 );
				if( err < tolerance * fabs(-log10(bin_tail)-logNT) * bin_tail ) break;
			}
		}

		return (-log10(bin_tail) - logNT);
	}

private:
    vt::vector<CVec2i>	m_vecFitMatT;
    vt::vector<UINT16>	m_vecFitVec;
    vt::CMtx2x2d	m_mtxATA;
    vt::CVec2d		m_mtxATV;
	
    vt::vector<vt::CVec3d> m_lineEqns;		//store the line Equation coefficients, vec3=[w1,w2,w3] for line w1*x + w2*y + w3=0;
    vt::vector<vt::CVec4f> m_lineEndPts;	//store the line endpoints, [x1,y1,x2,y3]
    vt::vector<float>  m_lineDirect;	//store the line direction

	float	m_flLineScalar;   // detection on scaled images, endpoints are mapped by the scalar
    vt::CRGBAImg m_imgScaled;

public:
	LineParams	m_params;
};

// Helper class 
struct DisjointSets
{
    HRESULT Init(int n);
    int     FindSet(int i);
    void    Union(int i, int j);
    void    Link(int i, int j);
    HRESULT SetIds(vt::vector<int>& ids);
    vt::vector<int> p;
    vt::vector<int> rank;
};
}
