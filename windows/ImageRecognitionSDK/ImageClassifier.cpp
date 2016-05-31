#include "ImageClassifier.h"
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

CImageClassifier::CImageClassifier()
{
}

void CImageClassifier::Release()
{    
}

CImageClassifier::~CImageClassifier()
{
    Release();
}

HRESULT CImageClassifier::LoadModel(const char *buffer, int buffer_size, int &bytes_read)
{
	class mem_stream_buf : public std::streambuf
	{
	public:
		mem_stream_buf(const char *buffer, int buffer_size)
		{
			this->setg(const_cast<char*>(buffer), const_cast<char*>(buffer), const_cast<char*>(buffer)+buffer_size);
		}

		int bytes_read()
		{
			return static_cast<int>(gptr() - eback());
		}
	};

	mem_stream_buf buf(buffer, buffer_size);
	std::istream stream(&buf);
	LoadModel(stream);
	bytes_read = buf.bytes_read();

	return S_OK;
}

HRESULT CImageClassifier::LoadModel(std::istream &stream)
{
	Release();

	m_sppPooler.Load(stream);
	m_netFcs.Load(stream);

	return S_OK;
}

HRESULT CImageClassifier::LoadModel(const std::wstring &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageClassifier::LoadModel(const std::string &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageClassifier::GetLayerDim(const char *layerName, int *channels, int *width, int *height)
{
	std::string strLayerName = layerName;
	m_netFcs.GetLayerDim(strLayerName, channels, width, height, int());

	return S_OK;
}

HRESULT CImageClassifier::GetLayerResponse(const char *layerName, float *buffer, int buffer_size)
{
	std::string strLayerName = layerName;
	const DNNTestLib::CPUData *acts = m_netFcs.GetLayerResponse(strLayerName);
	DNNTestLib::CPUData averageResponses;
	AverageResponses(*acts, averageResponses);
	averageResponses.GetData(buffer, buffer_size);

	return S_OK;
}

HRESULT CImageClassifier::GetLayerName(const int layerIdx, const char *&layerName) const
{
	layerName = m_netFcs.GetLayerName(layerIdx);

	return S_OK;
}

HRESULT CImageClassifier::GetTopNResults(
	const int imgWidth,
	const int imgHeight,
	const DNNTestLib::CPUData &convFeatureMap,
	const IUClassificationOption classificationOption,
	const int topN,
	IUClassificationInfo *ptrResults)
{   
	HRESULT hr = S_OK;

	std::vector<double> rects;
	PrepareClassificationRegions(imgWidth, imgHeight, classificationOption, rects);
	int numRegions = static_cast<int>(rects.size()) / 4;

	int numClasses = 0;
	m_netFcs.GetOutputDim(&numClasses, nullptr, nullptr);

	int fc_input_channels = 0;
	int fc_input_width = 0;
	int fc_input_height = 0;
	m_netFcs.GetInputDim(&fc_input_channels, &fc_input_width, &fc_input_height);

	TimerStart(classification_spooling);
	// spp
	m_sppPooler.SpatialPooling(&m_poolData, &convFeatureMap, numRegions, &rects[0]);
	m_poolData.Resize(numRegions, fc_input_channels, fc_input_width, fc_input_height, 0, 0, DNNTestLib::CPUData::CWHN, false);
	TimerShow(classification_spooling);

	// fc
	TimerStart(classification_fc);
	m_netFcs.TestCPUData(&m_poolData, &m_outputScores);
	assert(m_outputScores.GetPaddingX() == 0 && m_outputScores.GetPaddingY() == 0);
	TimerShow(classification_fc);

	// get average of regions
	DNNTestLib::CPUData averageResponses;
	AverageResponses(m_outputScores, averageResponses);
	float *ptrAverRsp = averageResponses.GetDataPointer();

	m_classificationResults.resize(numClasses);
	for (int i = 0; i < numClasses; ++i)
		m_classificationResults[i] = { i, ptrAverRsp[i] };

	std::sort(m_classificationResults.begin(), m_classificationResults.end(), [](IUClassificationInfo &i, IUClassificationInfo &j){return i.confidence > j.confidence; });
    
	// prepare the output
    for (int i = 0; i < min(topN, m_classificationResults.size()); ++i)
    {
		ptrResults[i] = m_classificationResults[i];
    }   /// i

    return hr;
}

HRESULT CImageClassifier::PrepareClassificationRegions(
	const int imgWidth,
	const int imgHeight,
	const IUClassificationOption classificationOption,
	std::vector<double> &rects)
{
	// prepare regions for classification
	std::vector<IURECTReal> boxes;
	switch (classificationOption.classificationView)
	{
	case IUClassification_FullView:
		boxes.resize(1);
		boxes[0] = { 0, 0, (float)imgWidth - 1, (float)imgHeight - 1 };
		break;
	case IUClassification_NineView:
		{
		boxes.resize(9);
		float cropSize = min(imgWidth, imgHeight) * classificationOption.cropRatio;
		float halfW = (imgWidth - cropSize) / 2;
		float halfH = (imgHeight - cropSize) / 2;
		for (int y = 0; y < 3; y++)
		for (int x = 0; x < 3; x++)
			boxes[y * 3 + x] = { halfW * x, halfH * y, halfW * x + cropSize, halfH * y + cropSize };
		}
		break;
	default:
		break;
	}

	// transform boxes to array
	const int numRects = (int)boxes.size();
	rects.resize(4 * numRects);
	for (int i = 0; i < numRects; ++i)
	{
		rects[4 * i + 0] = (double)boxes[i].left;
		rects[4 * i + 1] = (double)boxes[i].top;
		rects[4 * i + 2] = (double)boxes[i].right;
		rects[4 * i + 3] = (double)boxes[i].bottom;
	}

	return S_OK;
}

HRESULT CImageClassifier::AverageResponses(
	const DNNTestLib::CPUData &responses,
	DNNTestLib::CPUData &averageResponses)
{
	HRESULT hr = S_OK;

	DNNTestLib::ThrowException(DNNTestLib::InitError, "CImageClassifier::AverageResponses:No response is pre-computed", responses.GetNum() >= 1);
	DNNTestLib::ThrowException(DNNTestLib::InitError, "CImageClassifier::AverageResponses:Only support width == 1, height == 1", responses.GetWidth() == 1 && responses.GetHeight() == 1);

	averageResponses.Resize(1, responses.GetChannels(), responses.GetWidth(), responses.GetHeight(), 0, 0, DNNTestLib::CPUData::CWHN, true);

	if (responses.GetNum() == 1)
	{
		averageResponses.CopyDataFrom(&responses);
		return hr;
	}

	assert(responses.GetPaddingX() == 0 && responses.GetPaddingY() == 0);
	assert(responses.GetDataStorageType() == DNNTestLib::CPUData::CWHN);
	int n = responses.GetNum();
	int rstSize = responses.GetChannels() * responses.GetWidth() * responses.GetHeight();
	const float *ptrRsp = responses.GetDataPointer();
	float *ptrAverPsp = averageResponses.GetDataPointer();
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < rstSize; ++j)
		{
			ptrAverPsp[j] += ptrRsp[j];
		}
		ptrRsp += responses.GetChannels();
	}
	for (int j = 0; j < rstSize; ++j)
		ptrAverPsp[j] /= n;

	return hr;
}

#define BOOL2bool(value) ((value) ? true : false)


