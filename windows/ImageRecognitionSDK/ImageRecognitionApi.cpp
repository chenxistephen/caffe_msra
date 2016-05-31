#include <assert.h>
#include <Windows.h>
#include <string>

BOOL APIENTRY DllMain( HMODULE /*hModule*/,
                       DWORD  ul_reason_for_call,
                       LPVOID /*lpReserved*/
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

#define IMAGERECOGNITIONAPI_EXPORTS
#include "ImageRecognitionApi.h"

#include "ImageRecognition.h"

#define CHECK_HANDLE_PARAMETER(ptr) { assert(ptr != NULL); if(ptr == NULL) return E_INVALIDARG; }

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUCreateHandle(
    __out intptr_t& handle)
{
		handle = (intptr_t)(new CImageRecognition());

    return S_OK;
}

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUReleaseHandle(
    __in intptr_t handle)
{
    CHECK_HANDLE_PARAMETER(handle);

	delete ((CImageRecognition*)handle);

    return S_OK;
}

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IULoadModel(
	__in intptr_t handle,
	__in int tasks,
	__in const wchar_t* modelConfigFile)
{
    CHECK_HANDLE_PARAMETER(handle);

	return ((CImageRecognition*)handle)->LoadModel(std::wstring(modelConfigFile), tasks);
}

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUUnloadModel(
    __in intptr_t handle)
{
    CHECK_HANDLE_PARAMETER(handle);

	((CImageRecognition*)handle)->Release();

    return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationCategoryNumber(
__in intptr_t handle,
__out_ecount_full(1) int* pCategoryNumber)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pCategoryNumber);

	*pCategoryNumber = ((CImageRecognition*)handle)->ClassificationGetCategoryNumber();

	return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationTopNResults(
__in intptr_t handle,
__in int count,
__out_ecount_full(count) IUClassificationInfo* pCategories)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pCategories);

	return ((CImageRecognition*)handle)->ClassificationGetTopNResults(pCategories, count);
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerNum(
__in intptr_t handle,
__out_ecount_full(1) int *num)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(num);

	*num = ((CImageRecognition*)handle)->ClassificationGetLayerNum();

	return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerName(
__in intptr_t handle,
__in const int layerIdx,
__inout const char **layerName)
{
	CHECK_HANDLE_PARAMETER(handle);

	return ((CImageRecognition*)handle)->ClassificationGetLayerName(layerIdx, *layerName);
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerDim(
__in intptr_t handle,
__in const char *layerName,
__out_ecount_full(1) int *channels,
__out_ecount_full(1) int *width,
__out_ecount_full(1) int *height)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(channels);
	CHECK_HANDLE_PARAMETER(width);
	CHECK_HANDLE_PARAMETER(height);

	return ((CImageRecognition*)handle)->ClassificationGetLayerDim(layerName, channels, width, height);
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetClassificationLayerResponse(
__in intptr_t handle,
__in const char *layerName,
__in const int bufferSize,
__out_ecount_full(bufferSize) float *buffer)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(buffer);

	return ((CImageRecognition*)handle)->ClassificationGetLayerResponse(layerName, buffer, bufferSize);
}

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUGetDetectionCategoryNumber(
__in intptr_t handle, 
__out_ecount_full(1) int* pCategoryNumber)
{
    CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pCategoryNumber);

	*pCategoryNumber = ((CImageRecognition*)handle)->DetectionGetCategoryNumber();

    return S_OK;
}

IMAGERECOGNITIONAPI 
__checkReturn HRESULT IUDoTasks(
__in intptr_t handle,
__in_bcount(height * stride) const BYTE* pImage,
__in const int width,
__in const int height,
__in const int stride,
__in const int channel,
__in int tasks,
__in IUOption opts)
{
    CHECK_HANDLE_PARAMETER(handle);

	return ((CImageRecognition*)handle)->DoTasks(
		pImage,
        width, 
        height, 
        stride, 
        channel,
		tasks, 
		opts);
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetDetectedObjectCount(
__in intptr_t handle,
__out_ecount_full(1) int* pCount)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pCount);

	return ((CImageRecognition*)handle)->DetectionGetDetectedObjectCount(*pCount);

}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetDetectedObjects(
__in intptr_t handle,
__in int count,
__out_ecount_full(count) IUObjectInfo* pObjects)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pObjects);

	return ((CImageRecognition*)handle)->DetectionFillDetectedObjectInfos(pObjects, count);
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUSetDetectionPerClassThresholds(
__in intptr_t handle,
__in int count,
__in_ecount(count) float* pThresholds)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pThresholds);

	((CImageRecognition*)handle)->DetectionSetPerClassThresholds(pThresholds, count);

	return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationCategoryNumber(
__in intptr_t handle,
__out_ecount_full(1) int* pCategoryNumber)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pCategoryNumber);

	*pCategoryNumber = ((CImageRecognition*)handle)->SegmentationGetCategoryNumber();

	return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationResponseMapSize(
__in intptr_t handle,
__out_ecount_full(1) int *channels,
__out_ecount_full(1) int *width,
__out_ecount_full(1) int *height)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(channels);
	CHECK_HANDLE_PARAMETER(width);
	CHECK_HANDLE_PARAMETER(height);

	((CImageRecognition*)handle)->SegmentationGetResponseMapSize(*channels, *width, *height);

	return S_OK;
}

// channels is farest, then width, then height
IMAGERECOGNITIONAPI
__checkReturn HRESULT IUGetSegmentationResponseMap(
__in intptr_t handle,
__in int channels,
__in int width,
__in int height,
__out_ecount_full(channels*width*height) float* pResponseMap)
{
	CHECK_HANDLE_PARAMETER(handle);
	CHECK_HANDLE_PARAMETER(pResponseMap);

	((CImageRecognition*)handle)->SegmentationGetResponseMap(pResponseMap, channels, width, height);

	return S_OK;
}

IMAGERECOGNITIONAPI
__checkReturn HRESULT IUSetNumThreads(
__in intptr_t handle,
__in int numThreads)
{
	CHECK_HANDLE_PARAMETER(handle);

	((CImageRecognition*)handle)->SetNumThreads(numThreads);

	return S_OK;
}