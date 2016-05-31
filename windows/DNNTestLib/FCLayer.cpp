#include "FCLayer.h"
#include "NetUtils.h"
#include "LayerActs.h"
#include "Blas.h"

using namespace DNNTestLib;


FCLayer::FCLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_weight = CPUData::CreateFromDump(stream);
	m_biases = CPUData::CreateFromDump(stream);

	assert(m_weight->GetNum() == m_biases->GetNum());
	assert(m_weight->GetWidth() == 1);
	assert(m_weight->GetHeight() == 1);
	assert(m_biases->GetChannels() == 1);
	assert(m_biases->GetWidth() == 1);
	assert(m_biases->GetHeight() == 1);
#if USE_SSE
#if USE_BLAS
#else
	assert(m_weight->GetChannels() % 4 == 0);
#endif
#endif
}

void FCLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);
	m_weight->Dump(stream);
	m_biases->Dump(stream);
}

FCLayer::~FCLayer()
{
	SaveRelease(m_weight);
	SaveRelease(m_biases);
}

void FCLayer::InitActs()
{
	assert(GetInputNumber() == 1);
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);
	assert(m_min_pad_x == 0);
	assert(m_min_pad_y == 0);

	CPUData *input = m_inputs[0]->GetActs();
#if USE_BLAS
#else
	assert(input->GetNum() == 1);
#endif
	assert(input->GetChannels() == m_weight->GetChannels());
	assert(input->GetWidth() == 1);
	assert(input->GetHeight() == 1);
	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(input->GetNum(), m_weight->GetNum(), 1, 1, 0, 0);
	else
		m_layerActs->Resize(input->GetNum(), m_weight->GetNum(), 1, 1, 0, 0);
}

void FCLayer::Forward()
{
	CPUData *input = m_inputs[0]->GetActs();
	CPUData *output = m_layerActs;

	const int input_channels = input->GetChannels();
	const int output_channels = m_weight->GetNum();
	const float *data = input->GetDataPointer();
	float *target = output->GetDataPointer();
	float *weight = m_weight->GetDataPointer();
	float *biases = m_biases->GetDataPointer();
	const int input_num = input->GetNum();
	if (input_num < 10 && (output_channels % 4 == 0))  // blas is slower than sse, when input_num is small
	{
		for (int n = 0; n < input_num; ++n)
		{
			weight = m_weight->GetDataPointer();
			memcpy(target, biases, sizeof(float)*output_channels);
#if USE_SSE
			const __m128 *weight_m128 = (const __m128 *)weight;
			__m128 *target_m128 = (__m128 *)target;
			const int expend_length = 8;
			const int expend_loop = (output_channels >> 2) / expend_length;
			const int expend_remain = (output_channels >> 2) % expend_length;
#endif
			for (int i = 0; i < input_channels; i++)
			{
				const float v = data[i];
				if (v != 0)
				{
#if USE_SSE
					__m128 pixel_value = _mm_set_ps1(v);
					__m128 *r_base = target_m128;
					const __m128 *r_weight = weight_m128;
					for (int j = 0; j < expend_loop; j++)
					{
						r_base[0] = _mm_add_ps(r_base[0], _mm_mul_ps(pixel_value, r_weight[0]));
						r_base[1] = _mm_add_ps(r_base[1], _mm_mul_ps(pixel_value, r_weight[1]));
						r_base[2] = _mm_add_ps(r_base[2], _mm_mul_ps(pixel_value, r_weight[2]));
						r_base[3] = _mm_add_ps(r_base[3], _mm_mul_ps(pixel_value, r_weight[3]));
						r_base[4] = _mm_add_ps(r_base[4], _mm_mul_ps(pixel_value, r_weight[4]));
						r_base[5] = _mm_add_ps(r_base[5], _mm_mul_ps(pixel_value, r_weight[5]));
						r_base[6] = _mm_add_ps(r_base[6], _mm_mul_ps(pixel_value, r_weight[6]));
						r_base[7] = _mm_add_ps(r_base[7], _mm_mul_ps(pixel_value, r_weight[7]));
						r_weight += expend_length;
						r_base += expend_length;
					}
					for (int j = 0; j < expend_remain; j++)
						r_base[j] = _mm_add_ps(r_base[j], _mm_mul_ps(pixel_value, r_weight[j]));
#else
					for (int j = 0; j < output_channels; j++)
						target[j] += v * weight[j];
#endif
				}
#if USE_SSE
				weight_m128 += output_channels >> 2;
#else
				weight += output_channels;
#endif
			}
			target += output_channels;
			data += input_channels;
		}
	}
	else
	{
#if USE_BLAS
		for (int n = 0; n < input_num; n++)
			memcpy(target + n * output_channels, biases, sizeof(float)*output_channels);
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, input_num, output_channels, input_channels, 1.0f,
			data, input_channels, weight, output_channels, 1.0f, target, output_channels);
#else
		ThrowException(ConfigFileError, "Please set USE_BLAS to 1!", false);
#endif
	}
}
