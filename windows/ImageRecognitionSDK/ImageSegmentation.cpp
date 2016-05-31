#include "ImageSegmentation.h"
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


CImageSegmentation::CImageSegmentation()
{

}

CImageSegmentation::~CImageSegmentation()
{
	Release();
}

void CImageSegmentation::Release()
{
	m_net.Clear();

	m_poolData.Clear();
	m_outputData.Clear();
	m_outputScores.Clear();
}


HRESULT CImageSegmentation::LoadModel(const char *buffer, int buffer_size, int &bytes_read)
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

HRESULT CImageSegmentation::LoadModel(std::istream &stream)
{
	Release();

	m_net.Load(stream);

	// TODO: read in thresholds for each class
	int num_classes = 0;
	m_net.GetOutputDim(&num_classes, nullptr, nullptr);

	return S_OK;
}

HRESULT CImageSegmentation::LoadModel(const std::wstring &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageSegmentation::LoadModel(const std::string &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}


HRESULT CImageSegmentation::GetResponseMap(
	const int imgWidth,
	const int imgHeight,
	const DNNTestLib::CPUData &convFeatureMap,
	DNNTestLib::CPUData &responseMap)
{ 
	HRESULT hr = S_OK;

	m_net.TestCPUData(&convFeatureMap, &responseMap);

    return hr;
}
