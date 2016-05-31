//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for performing morphological operations
//
//  History:
//      2011/12/2-kramnath
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"
#include "vt_image.h"
#include "vt_pad.h"
#include "vt_convert.h"

namespace vt {

#ifdef MAKE_DOXYGEN_WORK
void foo();
#endif

/// \ingroup morph
/// <summary> This defines the choices in shape of the structuring element
/// that is used to perform the morphological erode and dilate
/// operations. Some of the common types are listed here.
/// </summary>
enum eStructuringElementType
{
    
    eStructuringElementTypeRectangle=0,  ///< rectangular structuring element
   
	eStructuringElementTypeCross=1,      ///< cross structuring element 
   
	eStructuringElementTypeDisk=2       ///< disk shaped structuring element
};

/// \ingroup morph
/// <summary> This defines the set of morphological operations that
/// are currently supported.
/// </summary>
enum eMorphologyOperation
{
   
    eMorphologyErode=0,      ///< Morphological erosion 
   
    eMorphologyDilate=1,     ///< Morphological dilation 
     
    eMorphologyOpen=2,      ///< Morphological open 
     
    eMorphologyClose=3      ///< Morphological close 
};


/// \ingroup morph
/// <summary> A structuring element defines the set of 'on' pixels that will
/// be considered for performing various morphological operations. Use this function
/// to create a structuring element that is one of the structuring element types. 
/// This will create a matrix of one of the predefined structuring elements for you.
/// </summary>
/// <param name="stElem"> Output structuring element matrix. </param>
/// <param name="rows"> Required number of rows in the structuring element. </param>
/// <param name="cols"> Required number of cols in the structuring element. </param>
/// <param name="eSet"> Required structuring element enum. </param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the structuring element type requested is eStructuringElementTypeDisk 
/// and the rows != cols (for disk structuring element type the region should be a square).
/// - E_OUTOFMEMORY if stElem allocation was needed and failed. 
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - For structuring element eStructuringElementTypeDisk , 'rows' should 
/// be equal to 'cols'.
/// - You do not have to create a structuring element by using this function call. You 
/// can create a custom structuring element by creating a Byte matrix of your choice.
/// This function is provided for convenience to pre-build some of the most common
/// structuring element types for you.
HRESULT
VtCreateStructuringElement( 
	OUT CMtx<Byte>& stElem, int rows, int cols, eStructuringElementType eSet );

/// \ingroup morph
/// <summary> Performs an erosion operation on the source image by using the 
/// specified structuring element. The erosion operation per pixel is defined 
/// as assigning the minimum pixel value of the neighbors as indicated by the 
/// structuring element. This operation can be repeated multiple times
/// and can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <param name="ex"> Required extend mode for padding the image based on ExtendMode.
/// Default value is TypeMax which pads the image with '+inf' for erosion.
/// Zero padding and Extend mode are also supported.</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported. 
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
HRESULT
VtErode( CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
        const CMtx<Byte>& stElem, ExtendMode ex, int iterations=1);


/// \ingroup morph
/// <summary> Performs an erosion operation on the source image by using the 
/// specified structuring element. The erosion operation per pixel is defined 
/// as assigning the minimum pixel value of the neighbors as indicated by the 
/// structuring element. This operation can be repeated multiple times
/// and can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported. 
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This function pads the source image with the default values (+inf) for erosion.
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
inline HRESULT
VtErode( CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
        const CMtx<Byte>& stElem, int iterations=1)
{
    return VtErode(dst, rctDst, src, ptSrcOrigin, stElem, TypeMax, iterations);
}


/// \ingroup morph
/// <summary> Performs a dilation operation on the source image by using the 
/// specified structuring element. The dilation operation per pixel is defined 
/// as assigning the maximum pixel value of the neighbors as indicated by the 
/// structuring element. This operation can be repeated multiple times
/// and can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <param name="ex"> Required extend mode for padding the image based on ExtendMode. 
/// Default value is TypeMin which pads the image with '-inf' for dilation.
/// Zero padding and Extend mode are also supported.</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported.
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
HRESULT
VtDilate( CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
        const CMtx<Byte>& stElem, ExtendMode ex , int iterations=1 );


/// \ingroup morph
/// <summary> Performs a dilation operation on the source image by using the 
/// specified structuring element. The dilation operation per pixel is defined 
/// as assigning the maximum pixel value of the neighbors as indicated by the 
/// structuring element. This operation can be repeated multiple times
/// and can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported.
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This function pads the source image with the default values (-inf) for dilation.
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
inline HRESULT
VtDilate( CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
        const CMtx<Byte>& stElem, int iterations=1)
{
    return VtDilate(dst, rctDst, src, ptSrcOrigin, stElem, TypeMin, iterations);
}

/// \ingroup morph
/// <summary> Performs a morphological open operation which is the dilation of
/// the erosion of the source image. It can be used to fill small holes. This operation
/// can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported.
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This function pads the source image with the default values (-inf) for dilation 
/// and (+inf) for erosion.
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
HRESULT
VtMorphologyOpen(
	CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
	  const CMtx<Byte>& stElem, int iterations=1);

/// \ingroup morph
/// <summary> Performs a morphological close operation which is the erosion of
/// the dilation of the source image. It can be used to remove small elements. This 
/// operation can be done in-place in this function.
/// </summary>
/// <param name="dst"> Output destination image. Will be created with the same type 
/// and size of the source image if un-initialized. dst can be equal to src (in-place).</param>
/// <param name="rctDst"> Indicates the location and dimensions of the output image. </param>
/// <param name="src"> Input source image </param>
/// <param name="ptSrcOrigin"> The origin of the input image </param>
/// <param name="stElem"> Input structuring element to perform the operation. </param>
/// <param name="iterations"> Required number of times the operation is to be applied
/// (default = 1).</param>
/// <returns> 
/// - S_OK on success.
/// - E_INVALIDARG if the source image is not initialized.
/// - E_INVALIDARG if the source image type is not supported.
/// - E_OUTOFMEMORY if dst image allocation was needed and failed.
/// - E_OUTOFMEMORY if temp memory allocation failed.
/// </returns>
/// <DL><DT> Remarks: </DT></DL>
/// - This function pads the source image with the default values (-inf) for dilation 
/// and (+inf) for erosion.
/// - This operation can be done 'in-place' so the src and dst images can be the same.
/// - This operation is supported for binary images and gray-scale images. For color images
/// or multi-channel images each channel is processed independently.
HRESULT
VtMorphologyClose(
CImg& dst, const CRect& rctDst, const CImg &src, CPoint ptSrcOrigin, 
	  const CMtx<Byte>& stElem, int iterations=1);


#ifndef VT_NO_XFORMS


/// \ingroup filtertransforms
/// <summary>Implementation of \link IImageTransform image transform \endlink that performs
/// a morphological operation.
/// </summary> 

class CMorphologyTransform: 
    public CImageTransformUnaryGeo<CMorphologyTransform, false>
{
	// IImageTransform implementation
public:
    void GetDstPixFormat(OUT int& frmtDst,
                         IN  const int* pfrmtSrcs, 
                         IN  UInt32  /*uSrcCnt*/);

	CRect GetRequiredSrcRect(const vt::CRect& rctDst);

	CRect GetAffectedDstRect(const vt::CRect& rctSrc);
    
	CRect GetResultingDstRect(const vt::CRect& rctSrc);

	HRESULT Transform(CImg* pimgDst, IN  const CRect& rctDst,
					  const CImg& imgSrc, const CPoint& ptSrc);

	virtual HRESULT Clone(ITaskState **ppState);

public:
    CMorphologyTransform()
    {}

	/// <summary> Initialize transform with a structuring element type </summary> 
	/// <param name="dstType">The image type that the transform will generate.</param> 
	/// <param name="rows">The number of rows in the structuring element.</param> 
	/// <param name="cols">The number of columns in the structuring element.</param> 
    /// <param name="eSet">The structuring element type.</param>
    /// <param name="eMorOp">The morphological operation. 
    /// that needs to be performed using the transform.</param>
	/// <returns> 
	///		- S_OK on success 
	///		- E_OUTOFMEMORY if memory allocation failed 
	/// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    /// - For eStructuringElementTypeDisk , 'rows' should be equal to 'cols'.
    HRESULT Initialize(int dstType, int rows, int cols, 
					   eStructuringElementType eSet, eMorphologyOperation eMorOp);

    /// <summary> Initialize transform with a structuring element matrix </summary> 
	/// <param name="dstType">The image type that the transform will generate.</param> 
	/// <param name="structuringElement">The structuring element matrix.</param> 
    /// <param name="eMorOp">The morphological operation.
    /// that needs to be performed using the transform.</param>
	/// <returns> 
	///		- S_OK on success 
	///		- E_OUTOFMEMORY if memory allocation failed 
	/// </returns>
    /// <DL><DT> Remarks: </DT></DL>
    HRESULT Initialize(int dstType, const CMtx<Byte>& structuringElement, 
					   eMorphologyOperation eMorOp);

protected:
    // structuring element used to perform morphology
    CMtx<Byte> m_strElem;
    // the type of the morphological operation desired (enum)
    eMorphologyOperation m_eMorOp;
    // the dst image type the transform will generate
    int m_dstType;
};

#endif

}; // end namespace vt
