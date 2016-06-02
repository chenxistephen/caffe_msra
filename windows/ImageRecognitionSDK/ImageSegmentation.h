#pragma once
#define NOMINMAX
#include <Windows.h>
#include <vector>

#include "vtcore.h"
#include "..\\DNNTestLib\\DNNTester.h"
#include "..\\DNNTestLib\\NetUtils.h"

#include "ImageRecognitionApi.h"

// wrapper class
class CImageSegmentation
{
public:
	CImageSegmentation();
	~CImageSegmentation();

	void Release();

	HRESULT LoadModel(const std::wstring &file);
	HRESULT LoadModel(const std::string &file);
	HRESULT LoadModel(std::istream &stream);
	HRESULT LoadModel(const char *buffer, int bufferSize, int &bytesRead);

	HRESULT GetResponseMap(
		const int imgWidth,
		const int imgHeight,
		const DNNTestLib::CPUData &convFeatureMap,
		DNNTestLib::CPUData &responseMap);

	int GetCategoryNumber() const
	{
		int channels, width, height;
		m_net.GetOutputDim(&channels, &width, &height);

		return channels;
	}

private:

private:
	DNNTestLib::DNNTester m_net;

	// intermediate data
	DNNTestLib::CPUData m_poolData, m_outputData, m_outputScores;

};