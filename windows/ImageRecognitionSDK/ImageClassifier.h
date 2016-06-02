#pragma once
#define NOMINMAX
#include <Windows.h>
#include <vector>

#include "vtcore.h"
#include "BoxProposal.h"
#include "..\\DNNTestLib\\DNNTester.h"
#include "..\\DNNTestLib\\NetUtils.h"
#include "..\\detection\\spm_pool\\spm_pool.h"

#include "ImageRecognitionApi.h"

/// wrapper class around the tester
class CImageClassifier
{
public:
    CImageClassifier();

    ~CImageClassifier();

    void Release();

	HRESULT LoadModel(const std::wstring &file);
	HRESULT LoadModel(const std::string &file);
	HRESULT LoadModel(std::istream &stream);
	HRESULT LoadModel(const char *buffer, int bufferSize, int &bytesRead);

	HRESULT GetTopNResults(
		const int imgWidth,
		const int imgHeight,
		const DNNTestLib::CPUData &convFeatureMap,
		const IUClassificationOption classificationOption,
        const int topN,
		IUClassificationInfo *ptrResults);

	HRESULT GetLayerName(const int layerIdx, const char *&layerName) const;

	HRESULT GetLayerDim(const char *layerName, int *channels, int *width, int *height);

	HRESULT GetLayerResponse(const char *layerName, float *buffer, int buffer_size);

	int GetLayerNum() const
	{
		int num;
		m_netFcs.GetLayerNum(&num);
		return num;
	}

    int GetCategoryNumber() const
    {
		int channels, width, height;
        m_netFcs.GetOutputDim(&channels, &width, &height);
		assert(width == 1);
		assert(height == 1);

		return channels;
    }

private:
	HRESULT PrepareClassificationRegions(
		const int imgWidth,
		const int imgHeight,
		const IUClassificationOption classificationOption,
		std::vector<double> &rects);

	HRESULT AverageResponses(
		const DNNTestLib::CPUData &responses,
		DNNTestLib::CPUData &averageResponses);


	DNNTestLib::CSpmPooler m_sppPooler;
	DNNTestLib::DNNTester m_netFcs;

	// intermediate data
	DNNTestLib::CPUData m_poolData, m_outputScores;
	std::vector<IUClassificationInfo> m_classificationResults;

	std::vector<std::string> m_fcLayerName;
};