#include "stdafx.h"

#include "vt_morphology.h"
#include <limits>

using namespace vt;


// Functors that are used for computing min/max
// and for initializing with the maxVal for a 
// particular type
template<typename T>
struct FindMin {
	T DoOp(T a, T b) const { return ( a > b) ? b : a; }
    T InitVal() const { return std::numeric_limits<T>::max(); }
};

template<typename T>
struct FindMax {
	T DoOp(T a, T b) const { return ( a > b) ? a : b; }
    T InitVal() const { return std::numeric_limits<T>::min(); }
};


// creates a structuring element
HRESULT
	vt::VtCreateStructuringElement( CMtx<Byte>& se, int rows, int cols, 
	eStructuringElementType eSet )
{
	VT_HR_BEGIN()

	// make them odd
	if(!(rows %2)) rows++;
	if(!(cols %2)) cols++;

	int rowCenter = (rows-1)/2;
	int colCenter = (cols-1)/2;

	VT_HR_EXIT( se.Create(rows, cols) );

	switch(eSet)
	{
	// Rectangular case, fill the matrix
	case eStructuringElementTypeRectangle:
		for(int i = 0; i < rows; ++i)
		{
			for(int j = 0; j < cols; ++j)
			{
				se(i,j) =  (Byte) 1;
			}
		}
		break;
	// create a cross
	case eStructuringElementTypeCross:

		for(int i = 0; i < rows; ++i)
		{
			for(int j = 0; j < cols; ++j)
			{
				if((i == rowCenter) || (j == colCenter))
					se(i,j) = (Byte) 1;
				else
					se(i,j) = (Byte) 0;
			}
		}
		break;
	// fill a circular region
	case eStructuringElementTypeDisk:
		if(rows != cols)
			VT_HR_EXIT(E_INVALIDARG);

		// sq of the radius
		int radSq = (rowCenter)*(rowCenter);

		for(int i = 0; i < rows; ++i)
		{
			for(int j = 0; j < rows; ++j)
			{
				if(((i-rowCenter)*(i-rowCenter) + (j-rowCenter)*(j-rowCenter)) <= radSq)
					se(i,j) = (Byte) 1;
				else
					se(i,j) = (Byte) 0;
			}
		}
        break;
    }

    VT_HR_END()
}


// Templated function that is called by the erode and dilate functions
// T specifies the type of the image, OpType is the type of the 
// class used for the comparison functions
template<typename T, typename OpType>
HRESULT
	MorphologyInternal(CImg& dst, const CRect& rctDst, const CTypedImg<T>& src, CPoint ptSrcOrigin,
	const CMtx<Byte>& stElem, int iters, const OpType& mathOp, ExtendMode ex)
{
	VT_HR_BEGIN();

    // pad size along each dimension
    int stelW = stElem.Cols();
    int stelH = stElem.Rows();
    int padW = (stelW - 1)/2;
    int padH = (stelH - 1)/2;

    int bands = src.Bands();

    // create image extend
	IMAGE_EXTEND imgEx;
    VT_HR_EXIT( imgEx.Initialize(ex, ex) );

    // Compute source required region
    CRect rctSrcBlk = rctDst;
	rctSrcBlk.InflateRect(padW, padH);
    rctSrcBlk -= ptSrcOrigin;
	
    int h = rctDst.Height();
    int w = rctDst.Width();

	// Temporary image used for processing
	CTypedImg<T> tmp;
    CRect rctSrc = src.Rect();
    if(rctSrc.RectInRect(rctSrcBlk))
    {
        // Src image already padded and setup
        VT_HR_EXIT( src.CopyTo(tmp, &rctSrcBlk) );
    }
    else
    {
	    // crop pad the temp image from the src
	    VT_HR_EXIT( VtCropPadImage(tmp, rctSrcBlk, src, imgEx) );
    }

	// create a rolling buffer which will be used for
	// doing the operation in-place
	CTypedImg<T> buffer;
	VT_HR_EXIT( buffer.Create(w, padH+1, src.Bands()) );

	// store
	int ps = tmp.PixSize();
	int sb = tmp.StrideBytes();

	// for number of iterations
	for(int itr = 0; itr < iters; ++itr)
	{
		// initialize
		int bufLoc = 0;
		bool buffFull = false;
		int itmsInBuff = 0;

		for(int y = 0; y < h; ++y)
		{
			for(int b = 0; b < bands; ++b)
			{
				// if the buffer is full roll the counter
				if(bufLoc == padH+1)  bufLoc = 0;
				// get ptrs
				T* pBuf = buffer.Ptr(0, bufLoc, b);
				Byte* pTmp1 = tmp.BytePtr(padW, y+padH, b);

				for(int x = 0; x < w; ++x, pBuf+=bands)
				{
					Byte* pTmp2 = (pTmp1 + static_cast<ptrdiff_t>(x*ps - padH*sb));
					// perform the min or max comparison for
					// on pixels in the stElem
					T val =  mathOp.InitVal();
					for(int yse = -padH; yse <= padH; ++yse, pTmp2+=sb)
					{
						for(int xse = -padW; xse <= padW; ++xse)
						{
							if(stElem(yse+padH, xse+padW))
							{
								T curval = *((T*) (pTmp2 + static_cast<ptrdiff_t>(xse*ps))); 
								val = mathOp.DoOp(val, curval);
							}
						}
						
					}
					// store the value in the buffer
					*pBuf = val;
				}

			}
			// increment number of items in buffer
			itmsInBuff++;
            bufLoc++;

            // check if buffer is full
            if(!buffFull && itmsInBuff == padH+1)
                buffFull = true;

			// once buffer is full, copy a line back to temp
			if(buffFull)
			{
				T* pbfTmp = tmp.Ptr(padW, y, 0);
				int ind = (bufLoc == padH+1) ? 0 : bufLoc;
				T* pbfBuf = buffer.Ptr(0, ind, 0);
                VtMemcpy((void *)pbfTmp, (void*) pbfBuf, w*ps);
				itmsInBuff--;
			}
		}

		// flush any remaining items from the buffer
		if(itmsInBuff)
		{
			int ind = bufLoc+1;
			for(; itmsInBuff > 0; ++ind, --itmsInBuff)
			{
				if(ind >= padH+1) ind = ind-padH-1;
				T* pbfTmp = tmp.Ptr(padW, h+padH-itmsInBuff, 0);
				T* pbfBuf = buffer.Ptr(0, ind, 0);
				VtMemcpy((void *)pbfTmp, (void*) pbfBuf, w*ps);
			}
		}
	}

	// copy the tmp image to dst image without the padding
	CRect tmpRect(padW, padH, w+padW, h+padH);
	CImg rectImg;
	VT_HR_EXIT( tmp.Share(rectImg, &tmpRect) );
	VT_HR_EXIT( VtConvertImage(dst, rectImg) );

	VT_HR_END()
}


template< template<typename> class OpType>
HRESULT OpTypeSwitcher(CImg& dst, const CRect& rctDst, const CImg &src, 
					   CPoint ptSrcOrigin, const CMtx<Byte>& stElem, 
					   int iterations, ExtendMode ex)
{
    
	VT_HR_BEGIN()

	if (!src.IsValid())
		VT_HR_EXIT(E_INVALIDSRC);

	// create the destination image if it needs to be
	VT_HR_EXIT( CreateImageForTransform(dst, rctDst.Width(), rctDst.Height(), 
										src.GetType()) );

	// determine which specialization to call
	switch( EL_FORMAT(src.GetType()) )
	{
	case EL_FORMAT_BYTE:
		// class
		OpType<Byte> fmb;
		// call
		MorphologyInternal<Byte, OpType<Byte> >(dst, rctDst, (CByteImg&) src, 
												ptSrcOrigin, stElem, iterations, 
												fmb,  ex);
		break;
	case EL_FORMAT_SHORT:
		OpType<UInt16> fms;
		MorphologyInternal<UInt16, OpType<UInt16> >(dst, rctDst, (CShortImg&) src, 
													ptSrcOrigin, stElem, iterations, 
													fms,  ex) ;
		break;
    case EL_FORMAT_INT:
		OpType<int> fmi;
		MorphologyInternal<int, OpType<int> >(dst, rctDst, (CIntImg&) src, ptSrcOrigin,
											  stElem, iterations, fmi,  ex) ;
		break;
	case EL_FORMAT_FLOAT:
		OpType<float> fmf;
		MorphologyInternal<float, OpType<float> >(dst, rctDst, (CFloatImg&) src, 
												  ptSrcOrigin, stElem, iterations, 
												  fmf,  ex);
		break;
	case EL_FORMAT_HALF_FLOAT:
		{
			// Convert half float to float, process and convert back
			CFloatImg floatSrcImg;
			CFloatImg floatDstImg;
			VT_HR_EXIT( VtConvertImage(floatSrcImg, src) );
			MorphologyInternal<float, OpType<float> >(
				floatDstImg, rctDst, (CFloatImg&) floatSrcImg, ptSrcOrigin, 
				stElem, iterations, fmf,   ex) ;
			VT_HR_EXIT( VtConvertImage(dst, floatDstImg) );
		}
		break;
	default:
		VT_HR_EXIT(E_INVALIDARG);
	}

	VT_HR_END()
}

// erosion
HRESULT
    vt::VtErode(CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
				const CMtx<Byte>& stElem, ExtendMode ex, int iterations)
{
    return OpTypeSwitcher<FindMin>(dst, rctDst, src, ptSrcOrigin, stElem, iterations, ex);
}

// dilation
HRESULT
	vt::VtDilate(CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
				 const CMtx<Byte>& stElem, ExtendMode ex, int iterations)
{
	 return OpTypeSwitcher<FindMax>(dst, rctDst, src, ptSrcOrigin, stElem, iterations, ex);
}

// open
HRESULT
	vt::VtMorphologyOpen(CImg& dst, const CRect& rctDst, const CImg& src, 
						 CPoint ptSrcOrigin, const CMtx<Byte>& se, int iterations)
{
	VT_HR_BEGIN()

    for(int i = 0; i < iterations; ++i)
    {
        CImg tmp;
	    // erode
	    VT_HR_EXIT( VtErode(tmp, rctDst, src, ptSrcOrigin, se, TypeMax, 1) );
	    // dilate
	    VT_HR_EXIT( VtDilate(dst, rctDst, tmp, ptSrcOrigin, se, TypeMin, 1) );
    }

	VT_HR_END()
}



HRESULT
    vt::VtMorphologyClose(CImg& dst, const CRect& rctDst, const CImg& src, 
						  CPoint ptSrcOrigin, const CMtx<Byte>& se, int iterations)
{
    VT_HR_BEGIN()

    for(int i = 0; i < iterations; ++i)
    {
        CImg tmp;
        // dilate
        VT_HR_EXIT( VtDilate(tmp, rctDst, src, ptSrcOrigin, se, TypeMin, 1) );
        // erode
        VT_HR_EXIT( VtErode(dst, rctDst, tmp, ptSrcOrigin, se, TypeMax, 1) );
    }
    VT_HR_END()
}



#ifndef VT_NO_XFORMS


// Tranform definitions

HRESULT CMorphologyTransform::Initialize(int dstType, int rows, int cols,
										 eStructuringElementType eSet,             
										 eMorphologyOperation eMorOp)
{
    m_dstType = dstType;
    m_eMorOp  = eMorOp;
    return vt::VtCreateStructuringElement(m_strElem, rows, cols, eSet);
}

HRESULT CMorphologyTransform::Initialize(int dstType, 
										 const CMtx<Byte>& structuringElement,
										 eMorphologyOperation eMorOp)
{
    m_eMorOp  = eMorOp;
    m_dstType = dstType;
    m_strElem = structuringElement;
    return m_strElem.GetError();
}

HRESULT CMorphologyTransform::Clone(ITaskState **ppState)
{
	return CloneTaskState<CMorphologyTransform>(ppState, 
		[this](CMorphologyTransform* pN)
    { return pN->Initialize(m_dstType, m_strElem, m_eMorOp); });
}

void CMorphologyTransform::GetDstPixFormat(OUT int& frmtDst,
                                           IN  const int* pfrmtSrcs, 
                                           IN  UInt32  /*uSrcCnt*/)
{   
	frmtDst = UpdateMutableTypeFields(m_dstType, *pfrmtSrcs);
	if( !IsValidConvertPair(*pfrmtSrcs, frmtDst) )
    {
		frmtDst = OBJ_UNDEFINED;
	}
}

CRect CMorphologyTransform::GetRequiredSrcRect(const vt::CRect& rctDst)
{ 
    // inflat based on size of the structuring element
    CRect rctRqdSrc = rctDst;
    int stelW = m_strElem.Cols();
    int stelH = m_strElem.Rows();
    int padW = (stelW - 1)/2;
    int padH = (stelH - 1)/2;
    rctRqdSrc.InflateRect(padW, padH);
    return rctRqdSrc;
}

CRect CMorphologyTransform::GetAffectedDstRect(const vt::CRect& rctSrc)
{
    // inflat based on size of the structuring element
    CRect rctAffDst = rctSrc;
    int stelW = m_strElem.Cols();
    int stelH = m_strElem.Rows();
    int padW = (stelW - 1)/2;
    int padH = (stelH - 1)/2;
    rctAffDst.InflateRect(padW,padH);
    return rctAffDst;
}

CRect CMorphologyTransform::GetResultingDstRect(const vt::CRect& rctSrc)
{
    return GetAffectedDstRect(rctSrc);
}

HRESULT CMorphologyTransform::Transform(CImg* pimgDst, IN  const CRect& rctDst,
                                        const CImg& imgSrc, const CPoint& ptSrc)
{
    VT_HR_BEGIN()
    
    // Call corresponding operations
    switch(m_eMorOp)
    {
    case eMorphologyErode:
        VT_HR_EXIT( VtErode(*pimgDst, rctDst, imgSrc, ptSrc, m_strElem) );
        break;
    case eMorphologyDilate:
        VT_HR_EXIT( VtDilate(*pimgDst, rctDst, imgSrc, ptSrc, m_strElem) );
        break;
    case eMorphologyOpen:
        VT_HR_EXIT( VtMorphologyOpen(*pimgDst, rctDst, imgSrc, ptSrc, m_strElem) );
        break;
    case eMorphologyClose:
        VT_HR_EXIT( VtMorphologyClose(*pimgDst, rctDst, imgSrc, ptSrc, m_strElem) );
        break;
    }

    VT_HR_END()
}

#endif