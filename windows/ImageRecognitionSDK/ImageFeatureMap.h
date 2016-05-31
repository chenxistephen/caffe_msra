#pragma once

#include <Windows.h>
#include <vector>
#include <iostream>

#include "vtcore.h"
#include "..\\DNNTestLib\\DNNTester.h"
#include "..\\DNNTestLib\\NetUtils.h"

#include "ImageRecognitionApi.h"

// wrapper class
class CImageFeatureMap
{
public:
	CImageFeatureMap(){}
	~CImageFeatureMap(){}

	void Release();

	HRESULT LoadModel(const std::wstring &file);
	HRESULT LoadModel(const std::string &file);
	HRESULT LoadModel(std::istream &stream);
	HRESULT LoadModel(const char *buffer, int bufferSize, int &bytesRead);

	HRESULT SetNumThreads(const int numThreads);

	HRESULT GetFeatureMap(const vt::CRGBByteImg &image,
		__out DNNTestLib::CPUData &convFeatureMap);

private:

private:
	DNNTestLib::DNNTester m_netConvs;
};