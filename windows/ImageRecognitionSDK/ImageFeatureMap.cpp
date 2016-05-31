#include "ImageFeatureMap.h"
#include "../DNNTestLib/DNNTester.h"
#include "../DNNTestLib/CPUData.h"
#include "../DNNTestLib/NetUtils.h"
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

void CImageFeatureMap::Release()
{
	m_netConvs.Clear();
}


HRESULT CImageFeatureMap::SetNumThreads(const int numThreads)
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



HRESULT CImageFeatureMap::LoadModel(const char *buffer, int buffer_size, int &bytes_read)
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

HRESULT CImageFeatureMap::LoadModel(std::istream &stream)
{
    m_netConvs.Load(stream);

    return S_OK;
}

HRESULT CImageFeatureMap::LoadModel(const std::wstring &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageFeatureMap::LoadModel(const std::string &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageFeatureMap::GetFeatureMap(const vt::CRGBByteImg &image,
	DNNTestLib::CPUData &convFeatureMap)
{ 
	TimerStart(Prepare_feature_map);
	m_netConvs.TestImage(image.BytePtr(), image.Width(), image.Height(), image.StrideBytes(), 3, NULL, &convFeatureMap);
	TimerShow(Prepare_feature_map);

    return S_OK;
}
