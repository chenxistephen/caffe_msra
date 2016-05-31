#include "Geometry.h"

#include <math.h>
#include <Windows.h>
#include <iostream>

namespace DNNTestLib
{

	float GetPointDistance(PointF pt1, PointF pt2)
	{
		float dx = pt1.X - pt2.X;
		float dy = pt1.Y - pt2.Y;
		return sqrt(dx*dx + dy*dy);
	}


	void Inverse(float mo[9], const float mi[9])
	{
		float f1 = mi[4] * mi[8] - mi[5] * mi[7];
		float f2 = mi[5] * mi[6] - mi[3] * mi[8];
		float f3 = mi[3] * mi[7] - mi[4] * mi[6];
		float del = mi[0] * f1 + mi[1] * f2 + mi[2] * f3;
		mo[0] = f1 / del;
		mo[1] = (mi[7] * mi[2] - mi[1] * mi[8]) / del;
		mo[2] = (mi[1] * mi[5] - mi[4] * mi[2]) / del;
		mo[3] = f2 / del;
		mo[4] = (mi[0] * mi[8] - mi[6] * mi[2]) / del;
		mo[5] = (mi[2] * mi[3] - mi[0] * mi[5]) / del;
		mo[6] = f3 / del;
		mo[7] = (mi[6] * mi[1] - mi[0] * mi[7]) / del;
		mo[8] = (mi[0] * mi[4] - mi[3] * mi[1]) / del;
	}


	void ComputeSimTransCoeff(float coeff[9], const float *pfShapeOri, const float *pfShapeNew, int iShapeSize)
	{
		int iPointNum = iShapeSize / 2;
		float a1 = 0, a3 = 0, a4 = 0;
		float c1 = 0, c2 = 0, c3 = 0, c4 = 0;
		int id = 0;
		for (int i = 0; i < iPointNum; i++)
		{
			float x = pfShapeOri[id];
			float y = pfShapeOri[id + 1];
			float xn = pfShapeNew[id];
			float yn = pfShapeNew[id + 1];
			a1 += x*x + y*y;
			a3 += x;
			a4 += y;
			c1 += x*xn + y*yn;
			c2 += y*xn - x*yn;
			c3 += xn;
			c4 += yn;
			id += 2;
		}
		float fPointNum = (float)iPointNum;
		a1 /= fPointNum;
		a3 /= fPointNum;
		a4 /= fPointNum;
		c1 /= fPointNum;
		c2 /= fPointNum;
		c3 /= fPointNum;
		c4 /= fPointNum;

		float t = a1 - a3*a3 - a4*a4;
		float b1 = (float)1 / t;
		float b2 = a1 / t;
		float b3 = -a3 / t;
		float b4 = -a4 / t;

		float a = b1*c1 + b3*c3 + b4*c4;
		float b = b1*c2 + b4*c3 - b3*c4;
		float dx = b3*c1 + b4*c2 + b2*c3;
		float dy = b4*c1 - b3*c2 + b2*c4;

		coeff[0] = a;
		coeff[1] = b;
		coeff[2] = dx;
		coeff[3] = -b;
		coeff[4] = a;
		coeff[5] = dy;
		coeff[6] = 0;
		coeff[7] = 0;
		coeff[8] = 1;
	}


	void ComputeSimTransCoeff(float coeff[9], const PointF *pfShapeOri, const PointF *pfShapeNew, int iPointNum)
	{
		ComputeSimTransCoeff(coeff, (float*)pfShapeOri, (float*)pfShapeNew, iPointNum * 2);
	}


	void AffineTransform(float *pfShapeNew, int iShapeSize, const float coeff[9], const float *pfShapeOri)
	{
		int nPoints = iShapeSize / 2;
		int id = 0;
		for (int i = 0; i < nPoints; i++)
		{
			pfShapeNew[id] = coeff[0] * pfShapeOri[id] + coeff[1] * pfShapeOri[id + 1] + coeff[2];
			pfShapeNew[id + 1] = coeff[3] * pfShapeOri[id] + coeff[4] * pfShapeOri[id + 1] + coeff[5];
			id += 2;
		}
	}


	void AffineTransform(PointF *pfShapeNew, int iPointNum, const float coeff[9], const PointF *pfShapeOri)
	{
		AffineTransform((float*)pfShapeNew, iPointNum * 2, coeff, (float*)pfShapeOri);
	}


	unsigned char Bilinear(const unsigned char *pInput, int iWidth, int iHeight, int iStride, float fX, float fY)
	{
		unsigned char val = 0;
		Bilinear(pInput, iWidth, iHeight, 1, iStride, fX, fY, &val);
		return val;
	}


	void Bilinear(const unsigned char *pInput, int iWidth, int iHeight, int iChannels, int iStride, 
		float fX, float fY, unsigned char *pOutput)
	{
		// out of range
		if (fX < 0.f || fX >(float)(iWidth - 1) || fY < 0 || fY >(float)(iHeight - 1))
		{
			memset(pOutput, 0, sizeof(unsigned char)*iChannels);
			return;
		}

		// int pixels
		if (fX == (float)(int)(fX) && fY == (float)(int)(fY))
		{
			memcpy(pOutput, pInput + iStride*(int)fY + iChannels*(int)fX, sizeof(unsigned char)*iChannels);
			return;
		}

		// otherwise
		int iLeft = (int)floor(fX);
		int iRight = iLeft + 1;
		int iTop = (int)floor(fY);
		int iBottom = iTop + 1;

		float dX0 = fX - iLeft;
		float dX1 = iRight - fX;
		float dY0 = fY - iTop;
		float dY1 = iBottom - fY;

		iRight = iRight >= iWidth ? iWidth - 1 : iRight;
		iBottom = iBottom >= iHeight ? iHeight - 1 : iBottom;

		const unsigned char *v00 = pInput + iTop*iStride + iLeft*iChannels;
		const unsigned char *v01 = pInput + iTop*iStride + iRight*iChannels;
		const unsigned char *v10 = pInput + iBottom*iStride + iLeft*iChannels;
		const unsigned char *v11 = pInput + iBottom*iStride + iRight*iChannels;
		for (int i = 0; i < iChannels; i++)
			pOutput[i] = (unsigned char)(dX1*dY1*v00[i] + dX1*dY0*v10[i] + dX0*dY0*v11[i] + dX0*dY1*v01[i] + 0.5f);
	}


	void ImageTransform(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcStride, const float coeffSrcToDst[9],
		unsigned char *pDstImage, int iDstWidth, int iDstHeight, int iDstStride)
	{
		float coeffDstToSrc[9];
		Inverse(coeffDstToSrc, coeffSrcToDst);
		for (int dy = 0; dy < iDstHeight; dy++)
		{
			for (int dx = 0; dx < iDstWidth; dx++)
			{
				float sx = coeffDstToSrc[0] * dx + coeffDstToSrc[1] * dy + coeffDstToSrc[2];
				float sy = coeffDstToSrc[3] * dx + coeffDstToSrc[4] * dy + coeffDstToSrc[5];
				pDstImage[dy*iDstStride + dx] = Bilinear(pSrcImage, iSrcWidth, iSrcHeight, iSrcStride, sx, sy);
			}
		}
	}

	void ImageTransform(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcChannels, int iSrcStride,
		const float coeffSrcToDst[9], unsigned char *pDstImage, int iDstWidth, int iDstHeight, int iDstStride)
	{
		float coeffDstToSrc[9];
		Inverse(coeffDstToSrc, coeffSrcToDst);
		for (int dy = 0; dy < iDstHeight; dy++)
		{
			unsigned char *pDst = pDstImage + dy * iDstStride;
			for (int dx = 0; dx < iDstWidth; dx++)
			{
				float sx = coeffDstToSrc[0] * dx + coeffDstToSrc[1] * dy + coeffDstToSrc[2];
				float sy = coeffDstToSrc[3] * dx + coeffDstToSrc[4] * dy + coeffDstToSrc[5];
				Bilinear(pSrcImage, iSrcWidth, iSrcHeight, iSrcChannels, iSrcStride, sx, sy, pDst);
				pDst += iSrcChannels;
			}
		}
	}


	bool ImageResize(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcStride, float fScale, 
		unsigned char **pDstImage, int *pDstWidth, int *pDstHeight)
	{
		const int FLOAT2INT_SCALE = 10;
		if (pSrcImage == NULL || pDstImage == NULL || pDstWidth == NULL || pDstHeight == NULL)
			return false;
		int iDstWidth = (int)((float)iSrcWidth*fScale+0.5f);
		int iDstHeight = (int)((float)iSrcHeight*fScale+0.5f);
		*pDstWidth = iDstWidth;
		*pDstHeight = iDstHeight;
		*pDstImage = NULL;
		unsigned char *pDst = new (std::nothrow) unsigned char [iDstWidth*iDstHeight];
		if (pDst == NULL)
			return false;
		int *pTempInt = new (std::nothrow) int [2*(iDstWidth+iDstHeight)];
		if (pTempInt == NULL)
		{
			delete[] pDst;
			pDst = NULL;
			return false;
		}
		*pDstImage = pDst;

		float   fScaleX = 1.f/fScale;
		float   fScaleY = 1.f/fScale;
		int     *pTableX = pTempInt;                    // size: 2*nDesWidth
		int     *pTableY = pTempInt + 2*iDstWidth;      // size: 2*nDesHeight

		// construct LUT
		for (int i = 0; i < iDstWidth; i++)
		{
			float x = (i + 0.5f) * fScaleX - 0.5f;
			int ix = (int)floor(x);     // integer part of x
			x -= ix;                        // fractional part of x
			if (ix >= iSrcWidth - 1)
			{
				x = 0;
				ix = iSrcWidth - 1;
			}
			else if (ix < 0)
			{
				x = 0;
				ix = 0;
			}
			pTableX[(i<<1)] = ix;                               // integer part of source x-coordinate
			pTableX[(i<<1)+1] = (int)(x * (1 << FLOAT2INT_SCALE)); // fractional part of source x-coordinate (in form of int)
		}
		for (int i = 0; i < iDstHeight; i++)
		{
			float   y = (i + 0.5f) * fScaleY - 0.5f;
			int     iy = (int)floor(y);

			y -= iy;
			if (iy >= iSrcHeight - 1)
			{
				y = 0;
				iy = iSrcHeight - 1;
			}
			else if (iy < 0)
			{
				y = 0;
				iy = 0;
			}
			pTableY[(i<<1)] = iy;                               // integer part of source y-coordinate
			pTableY[(i<<1)+1] = (int)(y * (1 << FLOAT2INT_SCALE)); // fractional part of source y-coordinate (in form of int)
		}

		// Resize
		for (int i = 0; i < iDstHeight; i++)
		{
			int     iy = pTableY[(i << 1)];       // integer part
			int     fy = pTableY[(i << 1) + 1];     // fractional part
			int     nDesOffset = i * iDstWidth;
			int     sy0 = iy;
			int     sy1 = sy0 + (fy > 0 && sy0 < iSrcHeight - 1);
			int     nSrcOffset0 = sy0  * iSrcStride;
			int     nSrcOffset1 = sy1 * iSrcStride;

			for (int j = 0; j < iDstWidth; j++)
			{
				int     ix = pTableX[(j << 1)];
				int     fx = pTableX[(j << 1) + 1];
				int     axy = fx * fy;
				int     w00 = (1 << (FLOAT2INT_SCALE << 1)) - (fx << FLOAT2INT_SCALE) - (fy << FLOAT2INT_SCALE) + axy;
				int     w01 = (fx << FLOAT2INT_SCALE) - axy;
				int     w10 = (fy << FLOAT2INT_SCALE) - axy;
				int     w11 = axy;
				int     sx0 = ix;
				int     sx1 = ix + (fx > 0 && ix < iSrcWidth -1);

				int nVal = (
					(w00 * pSrcImage[nSrcOffset0 + sx0] +
					w01 * pSrcImage[nSrcOffset0 + sx1] +
					w10 * pSrcImage[nSrcOffset1 + sx0] +
					w11 * pSrcImage[nSrcOffset1 + sx1])
					>> (FLOAT2INT_SCALE << 1));

				nVal = (nVal < 0? 0: nVal);
				nVal = (nVal > 255? 255: nVal);
				pDst[nDesOffset+j] = (BYTE)nVal;
			}
		}
		if (pTempInt != NULL)
		{
			delete[] pTempInt;
			pTempInt = NULL;
		}
		return true;
	}


	void ComputeScaleOffsetTransCoeff(float coeff[9], const float *pfShapeOri, const float *pfShapeNew, int iShapeSize)
	{
		int n = iShapeSize/2;
		float sx = 0, sy = 0, sxn = 0, syn = 0, sxxn = 0, syyn = 0, sxx = 0, syy = 0;
		for (int i = 0; i < n; i++)
		{
			float x = pfShapeOri[2*i];
			float y = pfShapeOri[2*i+1];
			float xn = pfShapeNew[2*i];
			float yn = pfShapeNew[2*i+1];
			sx += x;
			sy += y;
			sxn += xn;
			syn += yn;
			sxx += x*x;
			syy += y*y;
			sxxn += x*xn;
			syyn += y*yn;
		}
		float scale = (n*(sxxn+syyn)-sx*sxn-sy*syn)/(n*(sxx+syy)-sx*sx-sy*sy);
		float dx = (sxn - scale*sx)/n;
		float dy = (syn - scale*sy)/n;
		coeff[0] = scale;
		coeff[1] = 0;
		coeff[2] = dx;
		coeff[3] = 0;
		coeff[4] = scale;
		coeff[5] = dy;
		coeff[6] = 0;
		coeff[7] = 0;
		coeff[8] = 1.f;
	}


	void GetRotationCoeff(float coeff[9], float fTheta, float fCenterX, float fCenterY)
	{
		float fc = cos(fTheta);
		float fs = sin(fTheta);
		coeff[0] = fc;
		coeff[1] = -fs;
		coeff[2] = fCenterX-fCenterX*fc+fCenterY*fs;
		coeff[3] = fs;
		coeff[4] = fc;
		coeff[5] = fCenterY-fCenterX*fs-fCenterY*fc;
		coeff[6] = 0;
		coeff[7] = 0;
		coeff[8] = 1;
	}


	void Multiply(float m[9], const float m1[9], const float m2[9])
	{
		m[0] = m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6];
		m[1] = m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7];
		m[2] = m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8];

		m[3] = m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6];
		m[4] = m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7];
		m[5] = m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8];

		m[6] = m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6];
		m[7] = m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7];
		m[8] = m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8];

	}


	void GetRectFromPoints(RECT &rect, const float *pPoints, int iLen)
	{
		if (iLen < 10)
			return;

		const float pStdPoints[] = {
			0.290484f, 0.263211f,
			0.709854f, 0.262557f,
			0.500565f, 0.525255f,
			0.317176f, 0.751148f,
			0.684654f, 0.750561f};
		float pfCoeff[9];
		ComputeSimTransCoeff(pfCoeff, pStdPoints, pPoints, 10);
		float fCenterX = pfCoeff[0] * 0.5f + pfCoeff[1] * 0.5f + pfCoeff[2];
		float fCenterY = pfCoeff[3] * 0.5f + pfCoeff[4] * 0.5f + pfCoeff[5];
		float fTopLeftX = pfCoeff[2];
		float fTopLeftY = pfCoeff[5];
		float fTopRightX = pfCoeff[0] + pfCoeff[2];
		float fTopRightY = pfCoeff[3] + pfCoeff[5];
		float fWidth = sqrt((fTopLeftX - fTopRightX)*(fTopLeftX - fTopRightX) + (fTopLeftY - fTopRightY)*(fTopLeftY - fTopRightY));
		rect.left = int(fCenterX - fWidth / 2.f + 0.5f);
		rect.right = int(fCenterX + fWidth / 2.f + 0.5f);
		rect.top = int(fCenterY - fWidth / 2.f + 0.5f);
		rect.bottom = int(fCenterY + fWidth / 2.f + 0.5f);
	}
};
