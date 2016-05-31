#pragma once

#ifdef IMAGERECOGNITIONAPI_EXPORTS
#define IMAGERECOGNITIONAPI extern "C" __declspec(dllexport)
#else
#define IMAGERECOGNITIONAPI extern "C" __declspec(dllimport)
#endif

#include <Windows.h>

// support tasks (multiple choice)
#define IUTask_Classification 1
#define IUTask_Detection 2
#define IUTask_Segmentation 4

// classification view option
#define IUClassification_FullView  1							// faster but with slightly lower precision
#define IUClassification_NineView  2							// higher precision but slower

struct IUClassificationOption
{
	int classificationView = IUClassification_FullView;
	float cropRatio = 0.875;									// crop rectangle area ( cropRatio * short edge of image ) for NineViews
};

struct IUDetectionOption
{
	int maxNumPreFilter = 2000;									// max number of candidates for box proposal evaluation, faster with smaller value but lower performance
	int maxNumPostFilter = 500;									// max number of candidates for detection evaluation, faster with smaller value but lower performance
	float nmsThreshold = 0.3f;									// threshold for Non-Maximum Suppression, smaller value brings less output
};

struct IUOption
{
	IUClassificationOption classificationOption;
	IUDetectionOption detectionOption;
};

// for classification
struct IUClassificationInfo
{
	int categoryId;
	float confidence;
};

// for detection
struct IURECTReal
{
	float    left;
	float    top;
	float    right;
	float    bottom;
};

struct IUObjectInfo
{
	IURECTReal position;
    int categoryId;
    float confidence;
};

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUCreateHandle(
    __out intptr_t& handle);

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUReleaseHandle(
    __in intptr_t handle);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IULoadModel(
	__in intptr_t handle,
	__in int tasks,
	__in const wchar_t* modelFolder);

/// free memory up, can load model again afterwards
IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUUnloadModel(
    __in intptr_t handle);

/// no image resizing is required before calling the API,
/// as the image will be resized internally 
IMAGERECOGNITIONAPI
__checkReturn HRESULT IUDoTasks(
	__in intptr_t handle,
	__in_bcount(height * stride) const BYTE* pImage,
	__in const int width,
	__in const int height,
	__in const int stride,
	__in const int channel,
	__in int tasks, 
	__in IUOption opts);


IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationCategoryNumber(
__in intptr_t handle,
__out_ecount_full(1) int* pCategoryNumber);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationTopNResults(
__in intptr_t handle,
__in int count,
__out_ecount_full(count) IUClassificationInfo* pCategories);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerNum(
__in intptr_t handle,
__out_ecount_full(1) int *num);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerName(
__in intptr_t handle,
__in const int layerIdx,
__out const char **layerName);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerDim(
__in intptr_t handle,
__in const char *layerName,
__out_ecount_full(1) int *channels,
__out_ecount_full(1) int *width,
__out_ecount_full(1) int *height);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerResponse(
__in intptr_t handle,
__in const char *layerName,
__in const int bufferSize,
__out_ecount_full(bufferSize) float *buffer);

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUGetDetectionCategoryNumber(
__in intptr_t handle, 
__out_ecount_full(1) int* pCategoryNumber);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetDetectedObjectCount(
__in intptr_t handle,
__out_ecount_full(1) int* pCount);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetDetectedObjects(
__in intptr_t handle,
__in int count,
__out_ecount_full(count) IUObjectInfo* pObjects);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUSetDetectionPerClassThresholds(
__in intptr_t handle,
__in int count,
__in_ecount(count) float* pThresholds);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationCategoryNumber(
__in intptr_t handle,
__out_ecount_full(1) int* pCategoryNumber);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationResponseMapSize(
__in intptr_t handle,
__out_ecount_full(1) int *channels,
__out_ecount_full(1) int *width,
__out_ecount_full(1) int *height);

// channel is fastest, then width, then height
// idx = c + w * channels + h * width * channels
IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationResponseMap(
__in intptr_t handle,
__in int channels,
__in int width,
__in int height,
__out_ecount_full(channels*width*height) float* pResponseMap);

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUSetNumThreads(
__in intptr_t handle,
__in int numThreads);