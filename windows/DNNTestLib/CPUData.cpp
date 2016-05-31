
#include "NetUtils.h"
#include "CPUData.h"
#include "Blas.h"

#include <assert.h>

using namespace DNNTestLib;

#define MEM_ALIGN 32

static void Translate(float *p, int dim_fastest, int dim_second)
{
	ThrowException(ParameterError, "It will effect the program performance");
}

CPUData::CPUData() : m_data(NULL), m_num(0), m_channels(0),
m_width(0), m_height(0), m_padding_x(0), m_padding_y(0), m_capacity(0), m_dst(CPUData::CWHN)
{
}

CPUData::CPUData(const CPUData &)
{
}

const CPUData & CPUData::operator = (const CPUData &)
{
	return *this;
}

CPUData::~CPUData()
{
	Clear();
}

void CPUData::Clear()
{
	SaveAlignedRelease(m_data);
	m_num = 0;
	m_channels = 0;
	m_width = 0;
	m_height = 0;
	m_padding_x = 0;
	m_padding_y = 0;
	m_capacity = 0;
	m_dst = CPUData::CWHN;
}

void CPUData::PermuteData(DataStorageType dst)
{
	if (m_dst == dst)
		return;

	if (m_padding_x != 0 || m_padding_y != 0)
		ThrowException(ParameterError, "Undefined PermuteData with padding!=0");

	switch (m_dst)
	{
	case DNNTestLib::CPUData::CWHN:
		switch (dst)
		{
		case DNNTestLib::CPUData::NCWH:
			Translate(m_data, m_channels*m_width*m_height, m_num);
			break;
		default:
			ThrowException(ParameterError, "Undefined DataStorageType");
			break;
		}	
		break;
	case DNNTestLib::CPUData::NCWH:
		switch (dst)
		{
		case DNNTestLib::CPUData::CWHN:
			Translate(m_data, m_num, m_channels*m_width*m_height);
			break;
		default:
			ThrowException(ParameterError, "Undefined DataStorageType");
			break;
		}
		break;
		break;
	default:
		ThrowException(ParameterError, "Undefined DataStorageType");
		break;
	}

	m_dst = dst;
}

void CPUData::Resize(int num, int channels, int width, int height, int padding_x, int padding_y, DataStorageType dst /*= DataStorageType::CWHN*/, bool clear_mem /*= true*/)
{
	int size = num * channels * (width + 2 * padding_x) * (height + 2 * padding_y);
	if (size > m_capacity)
	{
		if (m_data == NULL)
		{
			m_data = reinterpret_cast<float*>(_aligned_malloc(size*sizeof(float), MEM_ALIGN));
		}
		else
		{
			m_data = reinterpret_cast<float*>(_aligned_realloc(m_data, size*sizeof(float), MEM_ALIGN));
		}
		ThrowException(MemoryAllocationError, "Unable to allocate Layer response", m_data != NULL);
		m_capacity = size;
	}

	if (clear_mem)
		memset(m_data, 0, size*sizeof(float));
	else
		PermuteData(dst);
	
	m_num = num;
	m_channels = channels;
	m_width = width;
	m_height = height;
	m_padding_x = padding_x;
	m_padding_y = padding_y;
}

int CPUData::GetOffset(int num, int channels, int width, int height) const
{
	ThrowException(ParameterError, "Out of range",
		num < m_num && num >= 0 && channels < m_channels && channels >= 0 &&
		width >= -m_padding_x && width < m_width + m_padding_x &&
		height >= -m_padding_y && height < m_height + m_padding_y);

	return Offset(num, channels, width, height);
}

int CPUData::Offset(int num, int channels, int width, int height) const
{
	switch (m_dst)
	{
	case DNNTestLib::CPUData::CWHN:
		return channels + ((width + m_padding_x) + ((height + m_padding_y) + num * (m_height + 2 * m_padding_y)) * (m_width + 2 * m_padding_x)) * m_channels;
	case DNNTestLib::CPUData::NCWH:
		return num + (channels + ((width + m_padding_x) + (height + m_padding_y) * (m_width + 2 * m_padding_x)) * m_channels) * m_num;
	default:
		ThrowException(ParameterError, "Undefined DataStorageType");
		break;
	}

	ThrowException(FATAL, "This sentence should never reached.");
	return -1;
}

CPUData::DataStorageType CPUData::GetDataStorageType() const
{
	return m_dst;
}

float *CPUData::GetDataPointer()
{
	return m_data + GetOffset(0, 0, 0, 0);
}

const float *CPUData::GetDataPointer() const
{
	return m_data + GetOffset(0, 0, 0, 0);
}

float *CPUData::GetDataPointer(int num, int channels, int width, int height)
{
	return m_data + GetOffset(num, channels, width, height);
}

const float *CPUData::GetDataPointer(int num, int channels, int width, int height) const
{
	return m_data + GetOffset(num, channels, width, height);
}

int CPUData::GetDataSize() const
{
	return m_num * m_channels * m_width * m_height;
}

int CPUData::GetDataSizeIncludePadding() const
{
	return m_num * m_channels * (m_width + 2 * m_padding_x) * (m_height + 2 * m_padding_y);
}

int CPUData::GetStrideH() const
{
	return Offset(0, 0, 0, 1) - Offset(0, 0, 0, 0);
}

void CPUData::SetData(const float *buffer, int buffer_size)
{
	ThrowException(ParameterError, "Buffer size is too small", buffer_size >= GetDataSize());
	if (m_padding_x == 0 && m_padding_y == 0)
	{
		memcpy(m_data, buffer, GetDataSize()*sizeof(float));
	}
	else
	{
		float *dst = GetDataPointer();
		int stride = GetStrideH();
		int size = m_width * m_channels * m_num;
		for (int y = 0; y < m_height; y++)
		{
			memcpy(dst + y*stride, buffer + y*size, size*sizeof(float));
		}
	}
}

void CPUData::GetData(float *buffer, int buffer_size) const
{
	ThrowException(ParameterError, "Buffer size is too small", buffer_size >= GetDataSize());
	if (m_padding_x == 0 && m_padding_y == 0)
	{
		memcpy(buffer, m_data, GetDataSize()*sizeof(float));
	}
	else
	{
		const float *src = GetDataPointer();
		int stride = GetStrideH();
		int size = m_width * m_channels * m_num;
		for (int y = 0; y < m_height; y++)
		{
			memcpy(buffer + y*size, src + y*stride, size*sizeof(float));
		}
	}
}

void CPUData::CopyDataFrom(const CPUData *p)
{
	assert(p->GetNum() == GetNum() && p->GetChannels() == GetChannels()
		&& p->GetWidth() == GetWidth() && p->GetHeight() == GetHeight());

	if (p->m_padding_x == 0 && p->m_padding_y == 0 && m_padding_x == 0 && m_padding_y == 0)
	{
		CPUData::DataStorageType dst_back = m_dst;
		m_dst = p->m_dst;
		memcpy(m_data, p->m_data, p->GetDataSize()*sizeof(float));
		PermuteData(dst_back);
	}
	else if (p->m_dst == CPUData::NCWH && m_dst == CPUData::NCWH)
	{
		const float *src = p->GetDataPointer();
		float *dst = GetDataPointer();
		const int src_stride = p->GetStrideH();
		const int dst_stride = GetStrideH();
		int size = m_width * m_channels * m_num;
		for (int y = 0; y < m_height; y++)
		{
			memcpy(dst + y*dst_stride, src + y*src_stride, size*sizeof(float));
		}
	}
	else if (p->m_dst == CPUData::CWHN && m_dst == CPUData::CWHN)
	{
		const int size = m_channels * m_width;
		for (int n = 0; n < m_num; n++)
		{
			const float *src = p->GetDataPointer(n, 0, 0, 0);
			float *dst = GetDataPointer(n, 0, 0, 0);
			const int src_stride = p->GetStrideH();
			const int dst_stride = GetStrideH();
			for (int y = 0; y < m_height; y++)
			{
				memcpy(dst + y * dst_stride, src + y * src_stride, sizeof(float)*size);
			}
		}
	}

	else
	{
		for (int n = 0; n < m_num; n++)
			for (int h = 0; h < m_height; h++)
				for (int w = 0; w < m_width; w++)
					for (int c = 0; c < m_channels; c++)
						*(m_data + Offset(n, c, w, h)) = *(p->m_data + p->Offset(n, c, w, h));
	}
}


int CPUData::GetNum() const
{
	return m_num;
}


int CPUData::GetChannels() const
{
	return m_channels;
}


int CPUData::GetWidth() const
{
	return m_width;
}


int CPUData::GetHeight() const
{
	return m_height;
}

int CPUData::GetPaddingX() const
{
	return m_padding_x;
}

int CPUData::GetPaddingY() const
{
	return m_padding_y;
}

void CPUData::Dump(std::ostream &stream) const
{
	WriteData(stream, m_num);
	WriteData(stream, m_channels);
	WriteData(stream, m_width);
	WriteData(stream, m_height);
	WriteData(stream, m_padding_x);
	WriteData(stream, m_padding_y);
	const float *src = GetDataPointer();
	int stride = GetStrideH();
	int size = m_width * m_channels * m_num;
	for (int y = 0; y < m_height; y++)
	{
		stream.write((const char *)src, size*sizeof(float));
		src += stride;
	}
}


CPUData* CPUData::CreateFromDump(std::istream &stream)
{
	int num = ReadData<int>(stream);
	int channels = ReadData<int>(stream);
	int width = ReadData<int>(stream);
	int height = ReadData<int>(stream);
	int padding_x = ReadData<int>(stream);
	int padding_y = ReadData<int>(stream);
	CPUData *p = CPUData::Create(num, channels, width, height, padding_x, padding_y, CPUData::NCWH);
	float *dst = p->GetDataPointer();
	int stride = p->GetStrideH();
	int size = p->m_width * p->m_channels * p->m_num;
	
	for (int y = 0; y < p->m_height; y++)
	{
		stream.read((char*)dst, size*sizeof(float));
		dst += stride;
	}
	return p;
}


CPUData* CPUData::Create(int num, int channels, int width, int height, int padding_x, int padding_y, DataStorageType dst /*= DataStorageType::CWHN*/)
{
	CPUData *p = new CPUData();
	p->Resize(num, channels, width, height, padding_x, padding_y);
	p->m_dst = dst;
	return p;
}
