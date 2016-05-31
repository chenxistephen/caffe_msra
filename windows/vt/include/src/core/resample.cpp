//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Routines for image warping
//
//  History:
//      2004/11/08-mattu
//          Created
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "vt_image.h"
#include "vt_filter.h"
#include "vt_convert.h"
#if !defined(VT_GCC)
#include "vt_layer.h"
#include "vt_readerwriter.h"
#include "vt_taskmanager.h"
#include "vt_transform.h"
#endif
#include "vt_resample.h"

using namespace vt;

//+-----------------------------------------------------------------------
// 
// Functions to point sample an image and return an interpolated value
//
// These functions all take an image with one or more bands and write the 
// sampled value into the array pRtn which has size equal to the number of bands 
// in the source image. The default value is used when the point (fX, fY) lies 
// outside the source image. The default pDefault can be NULL in which case the 
// elements of the return value are set to zero.
//
//------------------------------------------------------------------------
template <class T>
void vt::VtSampleNearest(const CTypedImg<T> &imgSrc, float fX, float fY, 
                         T *pDefault, T *pRtn)
{
    if(!imgSrc.IsValid() || pRtn==NULL)
        return;

    int iW = imgSrc.Width();
    int iH = imgSrc.Height();

    if(fX<0 || fY<0 || fX>iW-1 || fY>iH-1)
    {
        // outside the image
        if(pDefault)
            memcpy(pRtn, pDefault, imgSrc.PixSize());
        else
            memset(pRtn, 0, imgSrc.PixSize());
        return;
    }

    int iX = F2I(fX); // rounds to nearest under normal fpu mode
    int iY = F2I(fY);
    memcpy(pRtn, imgSrc.BytePtr(iX, iY), imgSrc.PixSize());
}

template <class T>
void vt::VtSampleBilinear(const CTypedImg<T> &imgSrc, float fX, float fY, 
                          T *pDefault, T *pRtn)
{
    if(!imgSrc.IsValid() || pRtn==NULL)
        return;

    int iW = imgSrc.Width();
    int iH = imgSrc.Height();

    if(fX<0 || fY<0 || fX>iW-1 || fY>iH-1)
    {
        // outside the image
        if(pDefault)
            memcpy(pRtn, pDefault, imgSrc.PixSize());
        else
            memset(pRtn, 0, imgSrc.PixSize());
        return;
    }

    int iXL = (int)fX; // get upper left coordinate
    int iYT = (int)fY;
    int iXR = iXL + 1;
    int iYB = iYT + 1;
    if(iXR==iW) // allow interpolation along right, bottom edge
        iXR--;
    if(iYB==iH)
        iYB--;
    float fAX = fX - iXL; // bilin factors
    float fAY = fY - iYT;
    float fAXAY = fAX * fAY;

    int iB;
    const T *pPixA = imgSrc.Ptr(iXL, iYT);
    const T *pPixB = imgSrc.Ptr(iXR, iYT);
    const T *pPixC = imgSrc.Ptr(iXL, iYB);
    const T *pPixD = imgSrc.Ptr(iXR, iYB);

    for(iB = 0; iB<imgSrc.Bands(); iB++)
    {
        float fA = (float)pPixA[iB];
        float fB = (float)pPixB[iB];
        float fC = (float)pPixC[iB];
        float fD = (float)pPixD[iB];
        float fQ = fA + fAX * (fB-fA) + fAY * (fC-fA) + fAXAY * (fA-fB-fC+fD);
        VtClip(pRtn + iB, fQ);
    }
}

template <class T>
void vt::VtSampleBilinearNoRangeTests(const CTypedImg<T> &imgSrc, float fX, float fY, T *pRtn)
{
    if(!imgSrc.IsValid() || pRtn==NULL)
        return;

    int iXL = (int)fX; // get upper left coordinate
    int iYT = (int)fY;
    int iXR = iXL + 1;
    int iYB = iYT + 1;
    float fAX = fX - iXL; // bilin factors
    float fAY = fY - iYT;
    float fAXAY = fAX * fAY;

    int iB;
    const T *pPixA = imgSrc.Ptr(iXL, iYT);
    const T *pPixB = imgSrc.Ptr(iXR, iYT);
    const T *pPixC = imgSrc.Ptr(iXL, iYB);
    const T *pPixD = imgSrc.Ptr(iXR, iYB);

    for(iB = 0; iB<imgSrc.Bands(); iB++)
    {
        float fA = (float)pPixA[iB];
        float fB = (float)pPixB[iB];
        float fC = (float)pPixC[iB];
        float fD = (float)pPixD[iB];
        float fQ = fA + fAX * (fB-fA) + fAY * (fC-fA) + fAXAY * (fA-fB-fC+fD);
        VtClip(pRtn + iB, fQ);
    }
}


// ensure that a set of values lie in the range (lower,upper) inclusive
static void VtRangeLimitSet(int &iX0, int &iX1, int &iX2, int iXLower, int iXUpper)
{
    if(iX0<iXLower)
        iX0 = iXLower;
    else if(iX0>iXUpper)
        iX0 = iXUpper;
    if(iX1<iXLower)
        iX1 = iXLower;
    else if(iX1>iXUpper)
        iX1 = iXUpper;
    if(iX2<iXLower)
        iX2 = iXLower;
    else if(iX2>iXUpper)
        iX2 = iXUpper;
}

template <class T>
void vt::VtSampleBicubic(const CTypedImg<T> &imgSrc, float fX, float fY, 
                         T *pDefault, T *pRtn)
{
    if(!imgSrc.IsValid() || pRtn==NULL)
        return;

    int iW = imgSrc.Width();
    int iH = imgSrc.Height();

    if(fX<0 || fY<0 || fX>iW-1 || fY>iH-1)
    {
        // outside the image
        if(pDefault)
            memcpy(pRtn, pDefault, imgSrc.PixSize());
        else
            memset(pRtn, 0, imgSrc.PixSize());
        return;
    }

    int iX1 = (int)fX; // get upper left coordinate of inner square
    int iY1 = (int)fY;
    int iX0 = iX1 - 1;
    int iY0 = iY1 - 1;
    int iX2 = iX1 + 1;
    int iY2 = iY1 + 1;
    int iX3 = iX1 + 2;
    int iY3 = iY1 + 2;

    // allow interpolation up to the image edges
    if(iX1<1 || iX1>iW-3)
        VtRangeLimitSet(iX0, iX2, iX3, 0, iW-1);
    if(iY1<1 || iY1>iH-3)
        VtRangeLimitSet(iY0, iY2, iY3, 0, iH-1);

    // bilin factors
    float fAX = fX - iX1;
    float fAY = fY - iY1;
    float fMAX = (1-fAX);
    float fMAY = (1-fAY);

    // bicubic weights for x
    float fPA = (1/6.0f) * (fMAX * fMAX - 1.0f) * fMAX;
    float fPB = 0.5f * (fAX * fMAX + 2.0f) * fMAX;
    float fPC = 0.5f * (fMAX * fAX + 2.0f) * fAX;
    float fPD = (1/6.0f) * (fAX * fAX - 1.0f) * fAX;

    // bicubic weights for y
    float fQA = (1/6.0f) * (fMAY * fMAY - 1.0f) * fMAY;
    float fQB = 0.5f * (fAY * fMAY + 2.0f) * fMAY;
    float fQC = 0.5f * (fMAY * fAY + 2.0f) * fAY;
    float fQD = (1/6.0f) * (fAY * fAY - 1.0f) * fAY;

    const T *pSrc0 = imgSrc.Ptr(iY0);
    const T *pSrc1 = imgSrc.Ptr(iY1);
    const T *pSrc2 = imgSrc.Ptr(iY2);
    const T *pSrc3 = imgSrc.Ptr(iY3);

    int iBands = imgSrc.Bands();

    iX0 *= iBands;
    iX1 *= iBands;
    iX2 *= iBands;
    iX3 *= iBands;
    int iX1End = iX1 + iBands - 1;

    for (;;)
    {
        float fY0 = fPA * pSrc0[iX0] + fPB * pSrc0[iX1] + fPC * pSrc0[iX2] + 
                    fPD * pSrc0[iX3];
        float fY1 = fPA * pSrc1[iX0] + fPB * pSrc1[iX1] + fPC * pSrc1[iX2] + 
                    fPD * pSrc1[iX3];
        float fY2 = fPA * pSrc2[iX0] + fPB * pSrc2[iX1] + fPC * pSrc2[iX2] + 
                    fPD * pSrc2[iX3];
        float fY3 = fPA * pSrc3[iX0] + fPB * pSrc3[iX1] + fPC * pSrc3[iX2] + 
                    fPD * pSrc3[iX3];

        VtClip(pRtn, fQA * fY0 + fQB * fY1 + fQC * fY2 + fQD * fY3);

        if(iX1==iX1End)
            break;

        iX0++; iX1++; iX2++; iX3++; pRtn++;
    }
}

template <class T>
void vt::VtSampleBicubicNoRangeTests(const CTypedImg<T> &imgSrc, float fX, float fY, T *pRtn)
{
    if(!imgSrc.IsValid() || pRtn==NULL)
        return;

    int iX1 = (int)fX; // get upper left coordinate of inner square
    int iY1 = (int)fY;
    int iX0 = iX1 - 1;
    int iY0 = iY1 - 1;
    int iX2 = iX1 + 1;
    int iY2 = iY1 + 1;
    int iX3 = iX1 + 2;
    int iY3 = iY1 + 2;

    // bilin factors
    float fAX = fX - iX1;
    float fAY = fY - iY1;
    float fMAX = (1-fAX);
    float fMAY = (1-fAY);

    // bicubic weights for x
    float fPA = (1/6.0f) * (fMAX * fMAX - 1.0f) * fMAX;
    float fPB = 0.5f * (fAX * fMAX + 2.0f) * fMAX;
    float fPC = 0.5f * (fMAX * fAX + 2.0f) * fAX;
    float fPD = (1/6.0f) * (fAX * fAX - 1.0f) * fAX;

    // bicubic weights for y
    float fQA = (1/6.0f) * (fMAY * fMAY - 1.0f) * fMAY;
    float fQB = 0.5f * (fAY * fMAY + 2.0f) * fMAY;
    float fQC = 0.5f * (fMAY * fAY + 2.0f) * fAY;
    float fQD = (1/6.0f) * (fAY * fAY - 1.0f) * fAY;

    const T *pSrc0 = imgSrc.Ptr(iY0);
    const T *pSrc1 = imgSrc.Ptr(iY1);
    const T *pSrc2 = imgSrc.Ptr(iY2);
    const T *pSrc3 = imgSrc.Ptr(iY3);

    int iBands = imgSrc.Bands();

    iX0 *= iBands;
    iX1 *= iBands;
    iX2 *= iBands;
    iX3 *= iBands;
    int iX1End = iX1 + iBands - 1;

    for (;;)
    {
        float fY0 = fPA * pSrc0[iX0] + fPB * pSrc0[iX1] + fPC * pSrc0[iX2] + 
            fPD * pSrc0[iX3];
        float fY1 = fPA * pSrc1[iX0] + fPB * pSrc1[iX1] + fPC * pSrc1[iX2] + 
            fPD * pSrc1[iX3];
        float fY2 = fPA * pSrc2[iX0] + fPB * pSrc2[iX1] + fPC * pSrc2[iX2] + 
            fPD * pSrc2[iX3];
        float fY3 = fPA * pSrc3[iX0] + fPB * pSrc3[iX1] + fPC * pSrc3[iX2] + 
            fPD * pSrc3[iX3];

        VtClip(pRtn, fQA * fY0 + fQB * fY1 + fQC * fY2 + fQD * fY3);

        if(iX1==iX1End)
            break;

        iX0++; iX1++; iX2++; iX3++; pRtn++;
    }
}

#define INST_POINTSAMPLE(type) \
template \
void vt::VtSampleNearest(const CTypedImg<type> &imgSrc, float fX, float fY,\
                         type *pDefault, type *pRtn);\
template \
void vt::VtSampleBilinear(const CTypedImg<type> &imgSrc, float fX, float fY,\
                          type *pDefault, type *pRtn);\
template \
void vt::VtSampleBicubic(const CTypedImg<type> &imgSrc, float fX, float fY,\
                         type *pDefault, type *pRtn);\
template \
void vt::VtSampleBilinearNoRangeTests(const CTypedImg<type> &imgSrc, float fX, float fY,\
                          type *pRtn);\
template \
void vt::VtSampleBicubicNoRangeTests(const CTypedImg<type> &imgSrc, float fX, float fY,\
                          type *pRtn);

INST_POINTSAMPLE(float)
INST_POINTSAMPLE(Byte)
INST_POINTSAMPLE(UInt16)

template <class T>
void vt::VtSampleRowBilinear(const CTypedImg<T> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount,
					T *pDefault, T *pRtn)
{
	int w = imgSrc.Width();
	int h = imgSrc.Height();
    int iX1 = (int)floor(fX);
    int iY1 = (int)floor(fY);
	float fXEnd = fX + fDX * (iCount-1);
	float fYEnd = fY + fDY * (iCount-1);
    int iX2 = (int)floor(fXEnd);
    int iY2 = (int)floor(fYEnd);
	if(iX1<1 || iX1>w-3 || iY1<1 || iY1>h-3 || iX2<1 || iX2>w-3 || iY2<1 || iY2>h-3)
	{
		// some outside the image
		int i;
		for(i=0; i<iCount; i++)
		{
			VtSampleBilinear(imgSrc, fX, fY, pDefault, pRtn + i * imgSrc.Bands());
			fX += fDX;
			fY += fDY;
		}
	}
	else
	{
		// avoid calling the function and doing per pixel tests
		int i;
		for(i=0; i<iCount; i++)
		{
			int iXL = (int)fX; // get upper left coordinate of inner square
			int iYT = (int)fY;
			int iXR = iXL + 1;
			int iYB = iYT + 1;

			float fAX = fX - iXL; // bilin factors
			float fAY = fY - iYT;
			float fAXAY = fAX * fAY;

			int iB;
			const T *pPixA = imgSrc.Ptr(iXL, iYT);
			const T *pPixB = imgSrc.Ptr(iXR, iYT);
			const T *pPixC = imgSrc.Ptr(iXL, iYB);
			const T *pPixD = imgSrc.Ptr(iXR, iYB);

			for(iB = 0; iB<imgSrc.Bands(); iB++, pRtn++)
			{
				float fA = (float)pPixA[iB];
				float fB = (float)pPixB[iB];
				float fC = (float)pPixC[iB];
				float fD = (float)pPixD[iB];
				float fQ = fA + fAX * (fB-fA) + fAY * (fC-fA) + fAXAY * (fA-fB-fC+fD);
				VtClip(pRtn, fQ);
			}

			fX += fDX;
			fY += fDY;
		}
	}
}

template <class T>
void vt::VtSampleRowBicubic(const CTypedImg<T> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount,
					T *pDefault, T *pRtn)
{
	int w = imgSrc.Width();
	int h = imgSrc.Height();
    int iXs = (int)floor(fX);
    int iYs = (int)floor(fY);
	float fXEnd = fX + fDX * (iCount-1);
	float fYEnd = fY + fDY * (iCount-1);
    int iXe = (int)floor(fXEnd);
    int iYe = (int)floor(fYEnd);
	if(iXs<2 || iXs>w-4 || iYs<2 || iYs>h-4 || iXe<2 || iXe>w-4 || iYe<2 || iYe>h-4)
	{
		// some potentially outside the image (including margin for floating point rounding errors)
		int i;
		for(i=0; i<iCount; i++)
		{
			VtSampleBicubic(imgSrc, fX, fY, pDefault, pRtn + i * imgSrc.Bands());
			fX += fDX;
			fY += fDY;
		}
	}
	else
	{
		// avoid calling the function and doing per pixel tests
		int i;
		for(i=0; i<iCount; i++)
		{
			int iX1 = (int)fX; // get upper left coordinate of inner square
			int iY1 = (int)fY;
			int iX0 = iX1 - 1;
			int iY0 = iY1 - 1;
			int iX2 = iX1 + 1;
			int iY2 = iY1 + 1;
			int iX3 = iX1 + 2;
			int iY3 = iY1 + 2;

			// bilin factors
			float fAX = fX - iX1;
			float fAY = fY - iY1;
			float fMAX = (1-fAX);
			float fMAY = (1-fAY);

			// bicubic weights for x
			float fPA = (1/6.0f) * (fMAX * fMAX - 1.0f) * fMAX;
			float fPB = 0.5f * (fAX * fMAX + 2.0f) * fMAX;
			float fPC = 0.5f * (fMAX * fAX + 2.0f) * fAX;
			float fPD = (1/6.0f) * (fAX * fAX - 1.0f) * fAX;

			// bicubic weights for y
			float fQA = (1/6.0f) * (fMAY * fMAY - 1.0f) * fMAY;
			float fQB = 0.5f * (fAY * fMAY + 2.0f) * fMAY;
			float fQC = 0.5f * (fMAY * fAY + 2.0f) * fAY;
			float fQD = (1/6.0f) * (fAY * fAY - 1.0f) * fAY;

			const T *pSrc0 = imgSrc.Ptr(iY0);
			const T *pSrc1 = imgSrc.Ptr(iY1);
			const T *pSrc2 = imgSrc.Ptr(iY2);
			const T *pSrc3 = imgSrc.Ptr(iY3);

			int iBands = imgSrc.Bands();

			iX0 *= iBands;
			iX1 *= iBands;
			iX2 *= iBands;
			iX3 *= iBands;
			int iX1End = iX1 + iBands - 1;

			for (;;)
			{
				// TODO: SSE
				float fY0 = fPA * pSrc0[iX0] + fPB * pSrc0[iX1] + fPC * pSrc0[iX2] + 
							fPD * pSrc0[iX3];
				float fY1 = fPA * pSrc1[iX0] + fPB * pSrc1[iX1] + fPC * pSrc1[iX2] + 
							fPD * pSrc1[iX3];
				float fY2 = fPA * pSrc2[iX0] + fPB * pSrc2[iX1] + fPC * pSrc2[iX2] + 
							fPD * pSrc2[iX3];
				float fY3 = fPA * pSrc3[iX0] + fPB * pSrc3[iX1] + fPC * pSrc3[iX2] + 
							fPD * pSrc3[iX3];

				VtClip(pRtn, fQA * fY0 + fQB * fY1 + fQC * fY2 + fQD * fY3);

				pRtn++;

				if(iX1==iX1End)
					break;

				iX0++; iX1++; iX2++; iX3++; 
			}

	     	fX += fDX;
			fY += fDY;
		}
	}
}

#define INST_ROWSAMPLE(type) \
template \
void vt::VtSampleRowBilinear(const CTypedImg<type> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount, \
	type *pDefault, type *pRtn);\
template \
void vt::VtSampleRowBicubic(const CTypedImg<type> &imgSrc, float fX, float fY, float fDX, float fDY, int iCount, \
	type *pDefault, type *pRtn);

INST_ROWSAMPLE(float)
INST_ROWSAMPLE(Byte)


template <class T>
void vt::VtSamplePatchBilinear(const CTypedImg<T> &imgSrc, const CMtx3x3f &mAff, CTypedImg<T> &imgDest, T *pDefault)
{
	int w = imgDest.Width();
	int h = imgDest.Height();
	int y;
	float fX = mAff(0,2);
	float fY = mAff(1,2);
	for(y=0; y<h; y++)
	{
		VtSampleRowBilinear(imgSrc, fX, fY, mAff(0,0), mAff(1,0), w, pDefault, imgDest.Ptr(y));
		fX += mAff(0,1);
		fY += mAff(1,1);
	}
}

template <class T>
void vt::VtSamplePatchBicubic(const CTypedImg<T> &imgSrc, const CMtx3x3f &mAff, CTypedImg<T> &imgDest, T *pDefault)
{
	int w = imgDest.Width();
	int h = imgDest.Height();
	int y;
	float fX = mAff(0,2);
	float fY = mAff(1,2);
	for(y=0; y<h; y++)
	{
		VtSampleRowBicubic(imgSrc, fX, fY, mAff(0,0), mAff(1,0), w, pDefault, imgDest.Ptr(y));
		fX += mAff(0,1);
		fY += mAff(1,1);
	}
}

#define INST_PATCHSAMPLE(type) \
template \
void vt::VtSamplePatchBilinear(const CTypedImg<type> &imgSrc, const CMtx3x3f &mHom, \
                               CTypedImg<type> &imgDest, type *pDefault);\
template \
void vt::VtSamplePatchBicubic(const CTypedImg<type> &imgSrc, const CMtx3x3f &mHom, \
                              CTypedImg<type> &imgDest, type *pDefault);

INST_PATCHSAMPLE(float)
INST_PATCHSAMPLE(Byte)

//+-----------------------------------------------------------------------------
// 
// Function: VtResampleNearest
// 
// Synopsis: resample an image using nearest-neighbor resampling.
// 
//------------------------------------------------------------------------------
void
vt::VtResampleNearest(const CImg& imgSrc, 
					  float sx, float tx, float sy, float ty, CImg& imgDst)
{
    VT_ASSERT( imgSrc.Bands()  == imgDst.Bands() );

	size_t iPixSize = imgSrc.PixSize();

	float fY = ty;
    for( int y = 0; y < imgDst.Height(); y++, fY+=sy )
    {
        int iSrcY = F2I( fY );
        iSrcY = VtMin(iSrcY, imgSrc.Height()-1);
        iSrcY = VtMax(iSrcY, 0);

        Byte* pDst = imgDst.BytePtr(y);
        const Byte* pSrc = imgSrc.BytePtr(iSrcY);

        float fX = tx;
        for( int x = 0; x < imgDst.Width(); x++, fX+=sx, pDst+=iPixSize )
        {
            // TODO: could pre-compute when the boundary conditions occur
            //       and move that into separate loops
			// TODO: multiply in inner loop is not good
            int iSrcX = F2I( fX );
            iSrcX = VtMin(iSrcX, imgSrc.Width()-1);
            iSrcX = VtMax(iSrcX, 0);

			memcpy(pDst, pSrc+iSrcX*iPixSize, iPixSize);
        }
    }
}
