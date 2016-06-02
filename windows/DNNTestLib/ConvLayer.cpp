#include "ConvLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"
#include "Blas.h"

using namespace DNNTestLib;

#if USE_BLAS
CPUData ConvLayer::m_im2col;
CPUData ConvLayer::m_output_buffer;
#endif

ConvLayer::ConvLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_padding_x = ReadData<int>(stream);
	m_padding_y = ReadData<int>(stream);
	m_stride_x = ReadData<int>(stream);
	m_stride_y = ReadData<int>(stream);
	m_weight = CPUData::CreateFromDump(stream);
	m_biases = CPUData::CreateFromDump(stream);

	// check parameters
#if USE_SSE
	assert(m_weight->GetNum() % 4 == 0);
#endif
	assert(m_biases->GetNum() == m_weight->GetNum());
	assert(m_biases->GetChannels() == 1);
	assert(m_biases->GetWidth() == 1);
	assert(m_biases->GetHeight() == 1);
}

void ConvLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);
	WriteData(stream, m_padding_x);
	WriteData(stream, m_padding_y);
	WriteData(stream, m_stride_x);
	WriteData(stream, m_stride_y);
	m_weight->PermuteData(CPUData::NCWH);
	m_weight->Dump(stream);
	m_biases->PermuteData(CPUData::NCWH);
	m_biases->Dump(stream);
}

ConvLayer::~ConvLayer()
{
	SaveRelease(m_weight);
	SaveRelease(m_biases);
}

int ConvLayer::GetInputPadX()
{
	return m_padding_x;
}

int ConvLayer::GetInputPadY()
{
	return m_padding_y;
}

void ConvLayer::InitActs()
{
	const CPUData *input = m_inputs[0]->GetActs();
	assert(input->GetDataStorageType() == CPUData::CWHN);
	assert(input->GetChannels() == m_weight->GetChannels());
	int output_width = (input->GetWidth() + 2 * m_padding_x - m_weight->GetWidth()) / m_stride_x + 1;
	int output_height = (input->GetHeight() + 2 * m_padding_y - m_weight->GetHeight()) / m_stride_y + 1;
	int num = input->GetNum();
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(num, m_weight->GetNum(), output_width, output_height, m_min_pad_x, m_min_pad_y, CPUData::CWHN);
	else
		m_layerActs->Resize(num, m_weight->GetNum(), output_width, output_height, m_min_pad_x, m_min_pad_y, CPUData::CWHN);
}

/////////////////////////////////////////
CPUData* ConvLayer::GetWeight() const
{
	return m_weight;
}

CPUData* ConvLayer::GetBias() const
{
	return m_biases;
}

void ConvLayer::CopyWeight(float *buffer, int buffer_size) const
{
	m_weight->GetData(buffer, buffer_size);
}

void ConvLayer::CopyBias(float *buffer, int buffer_size) const
{
	m_biases->GetData(buffer, buffer_size);
}
////////////////////////////////////////

void ConvLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	int output_width = output->GetWidth();
	int output_height = output->GetHeight();
	const float *weight = m_weight->GetDataPointer();
	const float *biases = m_biases->GetDataPointer();
	const int input_channels = input->GetChannels();
	const int input_stride = input->GetStrideH();
	const int filter_num = m_weight->GetNum();
	const int filter_width = m_weight->GetWidth();
	const int filter_height = m_weight->GetHeight();
	const int output_stride = output->GetStrideH();
#if USE_BLAS
	m_im2col.Resize(1, filter_width*filter_height*input_channels, output_width, output_height, 0, 0, CPUData::CWHN, false);
	float *im2col = m_im2col.GetDataPointer();
	const int im2col_stride = filter_width*filter_height*input_channels*output_width;
	m_output_buffer.Resize(1, filter_num, output_width, output_height, 0, 0, CPUData::CWHN, false);
	float *output_buffer = m_output_buffer.GetDataPointer();
	const int output_buffer_stride = m_output_buffer.GetStrideH();
	for (int n = 0; n < input->GetNum(); ++n)
	{
		const float *data = input->GetDataPointer(n, 0, -m_padding_x, -m_padding_y);
		float *target = output->GetDataPointer(n, 0, 0, 0);

		// im2col
		TransIm2Col(data, input_channels, input_stride, filter_width, filter_height, m_stride_x, m_stride_y,
			im2col, output_width, output_height, im2col_stride);

		// multipy
		for (int i = 0; i < output_width * output_height; i++)
			memcpy(output_buffer + i * filter_num, biases, filter_num*sizeof(float));
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			output_width*output_height, filter_num, filter_width*filter_height*input_channels, 1.0f,
			im2col, filter_width*filter_height*input_channels,
			weight, filter_num,
			1.0f, output_buffer, filter_num);

		// copy data
		for (int i = 0; i < output_height; i++)
			memcpy(target + i * output_stride, output_buffer + i * output_buffer_stride, sizeof(float)*filter_num*output_width);
	}


#else
	for (int n = 0; n < input->GetNum(); n++)
	{
		const float *data = input->GetDataPointer(n, 0, -m_padding_x, -m_padding_y);
		float *target = output->GetDataPointer(n, 0, 0, 0);
		for (int y = 0; y < output_height; y++)
		{
			for (int x = 0; x < output_width; x++)
			{
				ConvDotProduct(
					data + x * m_stride_x * input_channels,
					weight,
					biases,
					target + x * filter_num,
					input_channels,
					input_stride,
					filter_num,
					filter_width,
					filter_height);
			}
			data += m_stride_y * input_stride;
			target += output_stride;
		}
	}
	
#endif
}