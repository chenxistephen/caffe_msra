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

using DNNTestLib::RECTReal;

// wrapper class
class CImageObjectDetector
{
public:
	CImageObjectDetector();
	~CImageObjectDetector();

	void Release();

	HRESULT LoadModel(const std::wstring &file);
	HRESULT LoadModel(const std::string &file);
	HRESULT LoadModel(std::istream &stream);
	HRESULT LoadModel(const char *buffer, int bufferSize, int &bytesRead);


	HRESULT GetObjectInfoCount(
		const int imgWidth,
		const int imgHeight,
		const DNNTestLib::CPUData &convFeatureMap,
		int& count,
		const float coordinateScale);

	HRESULT FillObjectInfos(IUObjectInfo* pObjInfos,
        const int count);

	HRESULT SetThresholds(const float *pThresholds,
		const int count);

	HRESULT SetNumThreads(const int numThreads);

	HRESULT SetMaxNumPreFilter(const int num);
	HRESULT SetMaxNumPostFilter(const int num);
	HRESULT SetNMSThreshold(const float thres);

	int GetCategoryNumber() const
	{
		int channels, width, height;
		m_netFcs.GetOutputDim(&channels, &width, &height);
		assert(width == 1);
		assert(height == 1);
		assert(channels == static_cast<int>(m_classThresholds.size()));

		return channels;
	}

private:
	typedef HRESULT (CImageObjectDetector::*GetObjectInfoCountFun)(
		const int,
		const int,
		const DNNTestLib::CPUData &,
		int&,
		const float);
	HRESULT GetObjectInfoCount_rcnn(
		const int imgWidth,
		const int imgHeight,
		const DNNTestLib::CPUData &convFeatureMap,
		int& count,
		const float coordinateScale);
	HRESULT GetObjectInfoCount_frcn(
		const int imgWidth,
		const int imgHeight,
		const DNNTestLib::CPUData &convFeatureMap,
		int& count,
		const float coordinateScale);

	GetObjectInfoCountFun m_getObjectInfoCountFun;

public://private:
	DNNTestLib::DNNTester m_netFcs;
	DNNTestLib::CSpmPooler m_sppPooler;
	DNNTestLib::CBoxProposal m_boxProposal;

    std::vector<float> m_classThresholds;

	std::vector<IUObjectInfo> m_objInfos;

	float m_fNMSThreshold;

	// intermediate data
	DNNTestLib::CPUData m_poolData, m_outputData, m_outputScores, m_outputBoxes;
};