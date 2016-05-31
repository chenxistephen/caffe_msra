#include "ImageRecognition.h"
#include "ReadConfigure.h"
#include "Path.h"
#include "Directory.h"
#include "../DNNTestLib/Blas.h"
#include <algorithm>
#include <iostream>

#define ShowTime

#ifdef ShowTime
#include "HighResolutionTimer.h"
#define TimerStart(x) HighResolutionTimer timer_##x; timer_##x.Reset()
#define TimerShow(x) std::cout << #x##" time = " << timer_##x.Duration() << "ms" << endl
#else
#define TimerStart(x) 
#define TimerShow(x)
#endif

#define IUCheckTaskValidation(t) if (!IsTaskEnable(t)) {std::cout << #t << " must be initialized firstly!" << std::endl; return ERROR_INVALID_FUNCTION;}

void CImageRecognition::Release()
{
	m_featureMap.Release();
	m_classifier.Release();
	m_detector.Release();
	m_segmentation.Release();

	m_initedTasks = 0;

	m_convData.Clear();
	m_classificationInfos.clear();
	m_iDetectedObjectCount = 0;
	m_detectedObjects.clear();
	m_segmentationResponseMap.Clear();
}

HRESULT CImageRecognition::LoadModel(const std::wstring &modelConfigFile, const int tasks)
{
	std::string strModelConfigFile;
	wstring2string(modelConfigFile, strModelConfigFile);
	return LoadModel(strModelConfigFile, tasks);
}

HRESULT CImageRecognition::LoadModel(const std::string &modelConfigFile, const int tasks)
{
	if (tasks == 0)
		return E_INVALIDARG;

	Release();

	std::string strConfigModel = modelConfigFile;
	std::string strSharedConvolutionModel, strClassificationModel, strDetectionModel, strSegmentationModel;

	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "Missing modelConfigFile!", CDirectory::Exist(strConfigModel.c_str()))) return E_INVALIDARG;

	std::string strModelRefPath = CPath::GetDirectoryName(modelConfigFile.c_str());
	char strCurrentDirectoryBack[MAX_PATH];
	GetCurrentDirectoryA(MAX_PATH, strCurrentDirectoryBack);
	SetCurrentDirectoryA(strModelRefPath.c_str());

	// load config file
	BeginReadConfigure(CPath::GetFileName(strConfigModel));
	ReadVariable(m_imageSizeOps.bConstrainLongEdge);
	ReadVariable(m_imageSizeOps.iImgSizeConstraint);
	if (tasks)
	{
		ReadVariable(strSharedConvolutionModel); 
	}
	if ((tasks & IUTask_Classification) == IUTask_Classification)
	{
		ReadVariable(strClassificationModel);
	}
	if ((tasks & IUTask_Detection) == IUTask_Detection)
	{
		ReadVariable(strDetectionModel);
	}
	if ((tasks & IUTask_Segmentation) == IUTask_Segmentation)
	{
		ReadVariable(strSegmentationModel);
	}
	EndReadConfigure();

	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "m_imageSizeOps.iImgSizeConstraint must > 0!", (m_imageSizeOps.iImgSizeConstraint > 0))) return E_INVALIDARG;
	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "Missing SharedConvolution model config!", tasks && !strSharedConvolutionModel.empty())) return E_INVALIDARG;
	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "Missing Classification model config!", ((tasks & IUTask_Classification) != IUTask_Classification) || !strClassificationModel.empty())) return E_INVALIDARG;
	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "Missing Detection model config!", ((tasks & IUTask_Detection) != IUTask_Detection) || !strDetectionModel.empty())) return E_INVALIDARG;
	if (DNNTestLib::ThrowException(DNNTestLib::ConfigFileError, "Missing Segementation model config!", ((tasks & IUTask_Segmentation) != IUTask_Segmentation) || !strSegmentationModel.empty())) return E_INVALIDARG;

	// load models
	if (tasks)
	{
		m_featureMap.LoadModel(strSharedConvolutionModel);
	}
	if ((tasks & IUTask_Classification) == IUTask_Classification)
	{
		m_classifier.LoadModel(strClassificationModel);
	}
	if ((tasks & IUTask_Detection) == IUTask_Detection)
	{
		m_detector.LoadModel(strDetectionModel);
	}
	if ((tasks & IUTask_Segmentation) == IUTask_Segmentation)
	{
		m_segmentation.LoadModel(strSegmentationModel);
	}

	SetCurrentDirectoryA(strCurrentDirectoryBack);

	m_initedTasks = tasks;

    return S_OK;
}

bool CImageRecognition::IsTaskEnable(const int task) const
{
	return (m_initedTasks & task) == task;
}

HRESULT CImageRecognition::DoTasks(
	const BYTE* pImage,
	const int width,
	const int height,
	const int stride,
	const int channel, 
	const int tasks, 
	const IUOption opts)
{
	HRESULT hr = S_OK;
	if (tasks == 0)
		return hr;

	IUCheckTaskValidation(tasks);

	if (pImage == nullptr || channel < 3)
	{
		return E_INVALIDARG;
	}

	if (width == 0 || height == 0)
	{
		return E_INVALIDARG;
	}

	// get image
	TimerStart(FeedImg_Resize);
	vt::CRGBByteImg rgbImg;
	hr = FillVtRgbImage(pImage, width, height, stride, channel, rgbImg);
	if (FAILED(hr)) return hr;

	float scale = 1.f, coordinateScale = 1.f, xShift = 0.f, yShift = 0.f;

	//// if the ratio of image > 2 or < 0.5, crop the center part
	vt::CRGBByteImg croppedImg;
	vt::CRect cropRect = { 0, 0, rgbImg.Width(), rgbImg.Height() };
	const float fMaxRatio = 2.f;
	if (rgbImg.Height() / (float)rgbImg.Width() > fMaxRatio)
	{
		cropRect.top = static_cast<LONG>((rgbImg.Height() - rgbImg.Width() * fMaxRatio) / 2);
		cropRect.bottom = cropRect.top + static_cast<LONG>(rgbImg.Width() * fMaxRatio);
		rgbImg.CopyTo(croppedImg, &cropRect);
		yShift = static_cast<float>(cropRect.top);
	}
	else if (rgbImg.Width() / (float)rgbImg.Height() > fMaxRatio)
	{
		cropRect.left = static_cast<LONG>((rgbImg.Width() - rgbImg.Height() * fMaxRatio) / 2);
		cropRect.right = cropRect.left + static_cast<LONG>(rgbImg.Height() * fMaxRatio);
		rgbImg.CopyTo(croppedImg, &cropRect);
		xShift = static_cast<float>(cropRect.left);
	}
	else
	{
		rgbImg.Share(croppedImg);
	}

	// scale image
	vt::CRGBByteImg scaledImg;
	int imgSize = m_imageSizeOps.bConstrainLongEdge ? __max(croppedImg.Width(), croppedImg.Height()) : __min(croppedImg.Width(), croppedImg.Height());
	int targetImgSize = m_imageSizeOps.iImgSizeConstraint;

	if (imgSize == targetImgSize)
	{
		croppedImg.Share(scaledImg);
	}
	else
	{
		scale = float(imgSize) / targetImgSize;
		coordinateScale = float(imgSize - 1) / (targetImgSize - 1);

		const int sWidth = int(croppedImg.Width() / scale + 0.5);
		const int sHeight = int(croppedImg.Height() / scale + 0.5);

		scaledImg.Create(sWidth, sHeight);

		// use bilinear resize as training
		vt::VtResizeImage(scaledImg, scaledImg.Rect(), croppedImg, vt::eSamplerKernelBilinear);
	}
	TimerShow(FeedImg_Resize);

	// extract convolution feature map
	m_featureMap.GetFeatureMap(scaledImg, m_convData);

	if ((tasks & IUTask_Classification) == IUTask_Classification)
	{
		int categoryNum = static_cast<int>(m_classifier.GetCategoryNumber());
		m_classificationInfos.resize(categoryNum);
		m_classifier.GetTopNResults(scaledImg.Width(), scaledImg.Height(), m_convData, opts.classificationOption, categoryNum, &m_classificationInfos[0]);
	}
	if ((tasks & IUTask_Detection) == IUTask_Detection)
	{
		// set 
		DetectionSetMaxNumPreFilter(opts.detectionOption.maxNumPreFilter);
		DetectionSetMaxNumPostFilter(opts.detectionOption.maxNumPostFilter);
		DetectionSetNMSThreshold(opts.detectionOption.nmsThreshold);

		// run detection and get detection results with convolution feature map on resized image
		m_detector.GetObjectInfoCount(scaledImg.Width(), scaledImg.Height(), m_convData, m_iDetectedObjectCount, coordinateScale);

		if (m_iDetectedObjectCount > 0)
		{
			m_detectedObjects.resize(m_iDetectedObjectCount);
			m_detector.FillObjectInfos(&m_detectedObjects[0], m_iDetectedObjectCount);

			// scaled detected boxes from resized image to raw image
			for (int i = 0; i < m_iDetectedObjectCount; ++i)
			{
				//m_detectedObjects[i].position.left *= coordinateScale;
				m_detectedObjects[i].position.left += xShift;

				//m_detectedObjects[i].position.top *= coordinateScale;
				m_detectedObjects[i].position.top += yShift;

				//m_detectedObjects[i].position.right *= coordinateScale;
				m_detectedObjects[i].position.right += xShift;

				//m_detectedObjects[i].position.bottom *= coordinateScale;
				m_detectedObjects[i].position.bottom += yShift;
			}
		}
		else
			m_detectedObjects.clear();
	}

	if ((tasks & IUTask_Segmentation) == IUTask_Segmentation)
	{
		m_segmentation.GetResponseMap(scaledImg.Width(), scaledImg.Height(), m_convData, m_segmentationResponseMap);
	}

	return hr;
}

int CImageRecognition::ClassificationGetLayerNum() const
{
	IUCheckTaskValidation(IUTask_Classification);

	return m_classifier.GetLayerNum();
}

HRESULT CImageRecognition::ClassificationGetTopNResults(__out_ecount_full(count) IUClassificationInfo *ptrResults,
	const int count)
{
	if (ptrResults == nullptr)
		return E_INVALIDARG;

	IUCheckTaskValidation(IUTask_Classification);

	memcpy(ptrResults, &m_classificationInfos[0], min(count, (int)m_classificationInfos.size()) * sizeof(IUClassificationInfo));

	return S_OK;
}

HRESULT CImageRecognition::ClassificationGetLayerName(const int layerIdx, const char *&layerName) const
{
	IUCheckTaskValidation(IUTask_Classification);

	return m_classifier.GetLayerName(layerIdx, layerName);
}

HRESULT CImageRecognition::ClassificationGetLayerDim(const char *layerName, int *channels, int *width, int *height)
{
	IUCheckTaskValidation(IUTask_Classification);

	return m_classifier.GetLayerDim(layerName, channels, width, height);
}

HRESULT CImageRecognition::ClassificationGetLayerResponse(const char *layerName, float *buffer, int buffer_size)
{
	IUCheckTaskValidation(IUTask_Classification);

	return m_classifier.GetLayerResponse(layerName, buffer, buffer_size);
}

int CImageRecognition::ClassificationGetCategoryNumber() const
{
	IUCheckTaskValidation(IUTask_Classification);

	return m_classifier.GetCategoryNumber();
}

HRESULT CImageRecognition::DetectionGetDetectedObjectCount(int& count)
{
	IUCheckTaskValidation(IUTask_Detection);

	HRESULT hr = S_OK;

	count = m_iDetectedObjectCount;

	return hr;
}

HRESULT CImageRecognition::DetectionFillDetectedObjectInfos(__out_ecount_full(count) IUObjectInfo* pObjInfos,
	const int count)
{
	if (pObjInfos == nullptr)
		return E_INVALIDARG;

	IUCheckTaskValidation(IUTask_Detection);

	memcpy(pObjInfos, &m_detectedObjects[0], min(count, (int)m_detectedObjects.size()) * sizeof(IUObjectInfo));

	return S_OK;
}

HRESULT CImageRecognition::DetectionSetPerClassThresholds(const float *pThresholds,
	const int count)
{
	IUCheckTaskValidation(IUTask_Detection);

	return m_detector.SetThresholds(pThresholds, count);
}

int CImageRecognition::DetectionGetCategoryNumber() const
{
	IUCheckTaskValidation(IUTask_Detection);

	return m_detector.GetCategoryNumber();
}

HRESULT CImageRecognition::DetectionSetMaxNumPreFilter(const int num)
{
	IUCheckTaskValidation(IUTask_Detection);

	return m_detector.SetMaxNumPreFilter(num);
}

HRESULT CImageRecognition::DetectionSetMaxNumPostFilter(const int num)
{
	IUCheckTaskValidation(IUTask_Detection);

	return m_detector.SetMaxNumPostFilter(num);
}

HRESULT CImageRecognition::DetectionSetNMSThreshold(const float thres)
{
	IUCheckTaskValidation(IUTask_Detection);

	return m_detector.SetNMSThreshold(thres);
}

HRESULT CImageRecognition::SegmentationGetResponseMapSize(int &channels, int &width, int &height)
{
	IUCheckTaskValidation(IUTask_Segmentation);

	channels = m_segmentationResponseMap.GetChannels();
	width = m_segmentationResponseMap.GetWidth();
	height = m_segmentationResponseMap.GetHeight();

	return S_OK;
}

// channels is farest, then width, then height
HRESULT CImageRecognition::SegmentationGetResponseMap(__out_ecount_full(channels * width * height) float *pResponseMap,
	const int channels, const int width, const int height)
{
	IUCheckTaskValidation(IUTask_Segmentation);

	assert(m_segmentationResponseMap.GetNum() == 1);
	if (channels != m_segmentationResponseMap.GetChannels())
		return E_INVALIDARG;
	if (width != m_segmentationResponseMap.GetWidth())
		return E_INVALIDARG;
	if (height != m_segmentationResponseMap.GetHeight())
		return E_INVALIDARG;

	HRESULT hr = S_OK;

	DNNTestLib::CPUData denseResponseMap, *pCpuDataResponseMap = &m_segmentationResponseMap;

	if (m_segmentationResponseMap.GetPaddingX() != 0 || m_segmentationResponseMap.GetPaddingY() != 0)
	{
		denseResponseMap.Resize(m_segmentationResponseMap.GetNum(), m_segmentationResponseMap.GetChannels(), m_segmentationResponseMap.GetWidth(), m_segmentationResponseMap.GetHeight(),
			0, 0, DNNTestLib::CPUData::CWHN, false);
		denseResponseMap.CopyDataFrom(&m_segmentationResponseMap);
		pCpuDataResponseMap = &denseResponseMap;
	}

	memcpy(pResponseMap, pCpuDataResponseMap->GetDataPointer(), sizeof(float)*channels*width*height);

	return hr;
}

int CImageRecognition::SegmentationGetCategoryNumber() const
{
	IUCheckTaskValidation(IUTask_Segmentation);

	return m_segmentation.GetCategoryNumber();
}

HRESULT CImageRecognition::SetNumThreads(const int numThreads)
{
	if (numThreads < 1)
	{
		return E_INVALIDARG;
	}
#if USE_MKL
	mkl_set_num_threads(numThreads);
#else
	openblas_set_num_threads(numThreads);
#endif

	return S_OK;
}

HRESULT CImageRecognition::FillVtRgbImage(const BYTE* pImage,
	const int width,
	const int height,
	const int stride,
	const int channel,
	vt::CRGBImg& rgbImage)
{
	if (channel != 3 && channel != 4)
		return E_INVALIDARG;

	assert(channel == 3 || channel == 4);
	if (channel == 3)
	{
		rgbImage.Create((vt::Byte*)pImage, width, height, stride);
	}
	else
	{
		rgbImage.Create(width, height);

		for (int y = 0; y < height; ++y)
		{
			BYTE* dst = rgbImage.BytePtr(y);
			const BYTE* src = pImage + y * stride;

			for (int x = 0; x < width; ++x)
			{
				memcpy(dst, src, 3);    // b, g, r

				dst += 3;
				src += 4;
			}   // x
		}   // y
	}

	return S_OK;
}


bool CImageRecognition::wstring2string(__in const wstring &wstr, __out string &str)
{
	str.clear();
	DWORD dwMinSize = ::WideCharToMultiByte(CP_OEMCP, NULL, wstr.c_str(), -1, NULL, 0, NULL, FALSE);
	char *psz = new char[dwMinSize];
	::WideCharToMultiByte(CP_OEMCP, NULL, wstr.c_str(), -1, psz, dwMinSize, NULL, FALSE);
	str = psz;
	delete[] psz;
	return true;
}