#pragma once

#include <Windows.h>
#include <vector>
#include <string>

#include "vtcore.h"
#include "../DNNTestLib/DNNTester.h"
#include "../DNNTestLib/NetUtils.h"
#include "../DNNTestLib/CPUData.h"
#include "../detection/spm_pool/spm_pool.h"
#include "ImageFeatureMap.h"
#include "ImageClassifier.h"
#include "ImageObjectDetector.h"
#include "ImageSegmentation.h"
#include "ImageRecognitionApi.h"

// for extract feature map, the default value will be load in CImageFeatureMap and better not modified if not necessary
struct ImageSizeOpts
{
	bool bConstrainLongEdge;				// resize short/long edge to iImgSizeConstraint
	int iImgSizeConstraint;					// 
};

// wrapper class
class CImageRecognition
{
public:
	CImageRecognition() : m_iDetectedObjectCount(-1), m_initedTasks(0)
	{}
	~CImageRecognition()
	{
		Release();
	}

	void Release();

	bool IsTaskEnable(const int task) const;

	HRESULT LoadModel(const std::wstring &modelConfigFile, const int tasks);
	HRESULT LoadModel(const std::string &modelConfigFile, const int tasks);

	HRESULT DoTasks(
		const BYTE* pImage,
		const int width,
		const int height,
		const int stride,
		const int channel,
		const int tasks,
		const IUOption opts);

	int ClassificationGetCategoryNumber() const;

	HRESULT ClassificationGetTopNResults(__out_ecount_full(count) IUClassificationInfo *ptrResults,
		const int count);

	HRESULT ClassificationGetLayerName(const int layerIdx, const char *&layerName) const;

	HRESULT ClassificationGetLayerDim(const char *layerName, int *channels, int *width, int *height);

	HRESULT ClassificationGetLayerResponse(const char *layerName, float *buffer, int buffer_size);

	int ClassificationGetLayerNum() const;

    HRESULT DetectionGetDetectedObjectCount(int& count);

	HRESULT DetectionFillDetectedObjectInfos(__out_ecount_full(count) IUObjectInfo* pObjInfos,
        const int count);

	HRESULT DetectionSetPerClassThresholds(const float *pThresholds,
		const int count);

	HRESULT DetectionSetMaxNumPreFilter(const int num);
	HRESULT DetectionSetMaxNumPostFilter(const int num);
	HRESULT DetectionSetNMSThreshold(const float thres);

	int DetectionGetCategoryNumber() const;

	HRESULT SegmentationGetResponseMapSize(int &channels, int &width, int &height);

	HRESULT SegmentationGetResponseMap(__out_ecount_full(channels * width * height) float *pResponseMap,
		const int channels, const int width, const int height);

	int SegmentationGetCategoryNumber() const;

	HRESULT SetNumThreads(const int numThreads);

private:
	HRESULT FillVtRgbImage(const BYTE* pImage,
		const int width,
		const int height,
		const int stride,
		const int channel,
		vt::CRGBImg& rgbImage);

	static bool wstring2string(__in const std::wstring &wstr, __out std::string &str);

private:
	ImageSizeOpts m_imageSizeOps;

	CImageFeatureMap m_featureMap;
	CImageClassifier m_classifier;
	CImageObjectDetector m_detector;
	CImageSegmentation m_segmentation;

	int m_initedTasks;

	// intermediate data
	DNNTestLib::CPUData m_convData;   // conv feature map

	std::vector<IUClassificationInfo> m_classificationInfos;

	int m_iDetectedObjectCount;
	std::vector<IUObjectInfo> m_detectedObjects;

	DNNTestLib::CPUData m_segmentationResponseMap;  // response map for segmentation

};