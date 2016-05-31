#pragma once

#include <math.h>
#include <Windows.h>

namespace DNNTestLib
{
	struct PointF
	{
		float X;
		float Y;
	};

	float GetPointDistance(PointF pt1, PointF pt2);
	void Inverse(float mo[9], const float mi[9]);
	void Multiply(float m[9], const float m1[9], const float m2[9]);
	void ComputeSimTransCoeff(float coeff[9], const float *pfShapeOri, const float *pfShapeNew, int iShapeSize);
	void ComputeSimTransCoeff(float coeff[9], const PointF *pfShapeOri, const PointF *pfShapeNew, int iPointNum);
	void ComputeScaleOffsetTransCoeff(float coeff[9], const float *pfShapeOri, const float *pfShapeNew, int iShapeSize);
	void AffineTransform(float *pfShapeNew, int iShapeSize, const float coeff[9], const float *pfShapeOri);
	void AffineTransform(PointF *pfShapeNew, int iPointNum, const float coeff[9], const PointF *pfShapeOri);
	unsigned char Bilinear(const unsigned char *pInput, int iWidth, int iHeight, int iStride, float fX, float fY);
	void Bilinear(const unsigned char *pInput, int iWidth, int iHeight, int iChannels, int iStride,
		float fX, float fY, unsigned char *pOutput);
	void ImageTransform(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcStride, const float coeffSrcToDst[9],
		unsigned char *pDstImage, int iDstWidth, int iDstHeight, int iDstStride);
	void ImageTransform(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcChannels, int iSrcStride,
		const float coeffSrcToDst[9], unsigned char *pDstImage, int iDstWidth, int iDstHeight, int iDstStride);
	bool ImageResize(const unsigned char *pSrcImage, int iSrcWidth, int iSrcHeight, int iSrcStride, float fScale, 
		unsigned char **pDstImage, int *pDstWidth, int *pDstHeight);
	void GetRotationCoeff(float coeff[9], float fTheta, float fCenterX, float fCenterY);
	void GetRectFromPoints(RECT &rect, const float *pPoints, int iLen);
};