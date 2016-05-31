
#pragma once

#include <memory.h>
#include <xmmintrin.h>
#include <math.h>
#include <assert.h>
#include <string>
#include "Blas.h"

//#define USE_SSE				1

class MaxPoolActs
{
public:

	static std::string GetPoolType()
	{
		return "max";
	}

#if USE_SSE
	static void ActsSSE(const __m128 *data, __m128 *target, const int input_channels, const int input_stride,
		const int pool_size_x, const int pool_size_y)
	{
		memcpy(target, data, input_channels*sizeof(float));
		int channel_packs = input_channels >> 2;
		int stride_packs = input_stride >> 2;
		const __m128 *data_line = data;
		for (int y = 0; y < pool_size_y; y++)
		{
			for (int x = 0; x < pool_size_x; x++)
			{
				for (int i = 0; i < channel_packs; i++)
				{
					target[i] = _mm_max_ps(target[i], data_line[i]);
				}
				data_line += channel_packs;
			}
			data += stride_packs;
			data_line = data;
		}
	}

#else
	static void Acts(const float *data, float *target, const int input_channels, const int input_stride,
		const int pool_size_x, const int pool_size_y)
	{
		memcpy(target, data, input_channels*sizeof(float));
		const float *data_line = data;
		for (int y = 0; y < pool_size_y; y++)
		{
			for (int x = 0; x < pool_size_x; x++)
			{
				for (int i = 0; i < input_channels; i++)
				{
					if (target[i] < data_line[i])
					    target[i] = data_line[i];
						
				}
				data_line += input_channels;
			}
			data += input_stride;
			data_line = data;
		}
		//debug
		//if (1)
		//{
		//	for (int i = 0; i < input_channels; ++i)
		//	{
		//		std::cout << target[i] << " ";
		//	}
		//	std::cout << std::endl;
		//}
	}
#endif
};

class AvgPoolActs
{
public:

	static std::string GetPoolType()
	{
		return "avg";
	}

#if USE_SSE
	static void ActsSSE(const __m128 *data, __m128 *target, const int input_channels, const int input_stride,
		const int pool_size_x, const int pool_size_y)
	{
		memset(target, 0, input_channels*sizeof(float));
		int channel_packs = input_channels >> 2;
		int stride_packs = input_stride >> 2;
		const __m128 *data_line = data;
		for (int y = 0; y < pool_size_y; y++)
		{
			for (int x = 0; x < pool_size_x; x++)
			{
				for (int i = 0; i < channel_packs; i++)
				{
					target[i] = _mm_add_ps(target[i], data_line[i]);
				}
				data_line += channel_packs;
			}
			data += stride_packs;
			data_line = data;
		}

		// div by window_pixels
		__m128 div = _mm_set_ps1((float)(pool_size_x*pool_size_y));
		for (int fl_pack = 0; fl_pack<channel_packs; fl_pack++)
		{
			target[fl_pack] = _mm_div_ps(target[fl_pack], div);
		}
	}
#else
	static void Acts(const float *data, float *target, const int input_channels, const int input_stride,
		const int pool_size_x, const int pool_size_y)
	{
		memset(target, 0, input_channels*sizeof(float));
		const float *data_line = data;
		for (int y = 0; y < pool_size_y; y++)
		{
			for (int x = 0; x < pool_size_x; x++)
			{
				for (int i = 0; i < input_channels; i++)
				{
					target[i] += data_line[i];
				}
				data_line += input_channels;
			}
			data += input_stride;
			data_line = data;
		}
		float div = (float)(pool_size_x*pool_size_y);
		for (int i = 0; i < input_channels; i++)
			target[i] /= div;
	}
#endif
};

class ReluNeuronActs
{
public:

	static std::string GetNeuronType()
	{
		return "relu";
	}

#if USE_SSE
	static __m128 ActsSSE(const __m128 &data)
	{
		__m128 zeros = { 0, 0, 0, 0 };
		return _mm_max_ps(data, zeros);
	}
#else
	static float Acts(const float &data)
	{
		return data > 0 ? data : 0;
	}
#endif
};


class SigmoidNeuronActs
{
public:

	static std::string GetNeuronType()
	{
		return "sigmoid";
	}

#if USE_SSE
	static __m128 ActsSSE(const __m128 &data)
	{
		__m128 result;
		result.m128_f32[0] = 1.f / (1.f + exp(-data.m128_f32[0]));
		result.m128_f32[1] = 1.f / (1.f + exp(-data.m128_f32[1]));
		result.m128_f32[2] = 1.f / (1.f + exp(-data.m128_f32[2]));
		result.m128_f32[3] = 1.f / (1.f + exp(-data.m128_f32[3]));
		return result;
	}
#else
	static float Acts(const float &data)
	{
		return 1.f / (1.f + exp(-data));
	}
#endif
};


inline void ConvDotProduct(const float *data,
	const float *weight, const float *biases, float *target,
	const int input_channels, const int input_stride,
	const int filter_num, const int filter_width, const int filter_height)
{
	memcpy(target, biases, filter_num*sizeof(float));
	const int filter_width_channels = filter_width * input_channels;
#if USE_SSE
	const __m128 *weight_m128 = (const __m128 *)weight;
	__m128 *target_m128 = (__m128 *)target;
	const int expend_length = 8;
	const int expend_loop = (filter_num >> 2) / expend_length;
	const int expend_remain = (filter_num >> 2) % expend_length;
#endif
	for (int y = 0; y < filter_height; y++)
	{
		for (int x = 0; x < filter_width_channels; x++)
		{
			const float v = data[x];
			if (v != 0)
			{
#if USE_SSE
				__m128 pixel_value = _mm_set_ps1(v);
				__m128 *r_base = target_m128;
				const __m128 *r_weight = weight_m128;
				for (int i = 0; i < expend_loop; i++)
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
				for (int i = 0; i < expend_remain; i++)
					r_base[i] = _mm_add_ps(r_base[i], _mm_mul_ps(pixel_value, r_weight[i]));
#else
				for (int i = 0; i < filter_num; i++)
				{
					target[i] += v * weight[i];
				}
#endif
			}
#if USE_SSE
			weight_m128 += filter_num >> 2;
#else
			weight += filter_num;
#endif
		}
		data += input_stride;
	}
}


#if USE_SSE

inline void WithinChannelNormActs(const __m128 *data, __m128 *target, const __m128 *center, const int local_x, const int local_y,
	const int channels, const int input_stride, const float alpha, const float beta)
{
	int channel_packs = channels >> 2;
	memset(target, 0, channel_packs*sizeof(float));

	int stride_packs = input_stride >> 2;
	const __m128 *data_line = data;
	for (int y = 0; y < local_y; y++)
	{
		for (int x = 0; x < local_x; x++)
		{
			for (int i = 0; i < channel_packs; i++)
			{
				__m128 t = data_line[i];
				target[i] = _mm_add_ps(target[i], _mm_mul_ps(t, t));
			}
			data_line += channel_packs;
		}
		data += stride_packs;
		data_line = data;
	}

	const __m128 op1 = { 1.0f, 1.0f, 1.0f, 1.0f };
	__m128 mul = _mm_set_ps1(alpha / (float)(local_x*local_y));
	for (int i = 0; i < channel_packs; i++)
	{
		__m128 t = _mm_add_ps(op1, _mm_mul_ps(target[i], mul));
		t.m128_f32[0] = ::pow(t.m128_f32[0], -beta);
		t.m128_f32[1] = ::pow(t.m128_f32[1], -beta);
		t.m128_f32[2] = ::pow(t.m128_f32[2], -beta);
		t.m128_f32[3] = ::pow(t.m128_f32[3], -beta);
		target[i] = _mm_mul_ps(center[i], t);
	}
}

#else

inline void WithinChannelNormActs(const float *data, float *target, const float *center, const int local_x, const int local_y,
	const int channels, const int input_stride, const float alpha, const float beta)
{
	memset(target, 0, channels*sizeof(float));
	const float *data_line = data;
	for (int y = 0; y < local_y; y++)
	{
		for (int x = 0; x < local_x; x++)
		{
			for (int i = 0; i < channels; i++)
			{
				const float t = data_line[i];
				target[i] += t * t;
			}
			data_line += channels;
		}
		data += input_stride;
		data_line = data;
	}

	const float scale = alpha / (float)(local_x*local_y);
	for (int i = 0; i < channels; i++)
	{
		target[i] = ::pow(1.f + target[i] * scale, -beta) * center[i];
	}
}

#endif

inline void AcrossChannelNormActs(const float *data, float *target, float *scale_vec, const int channels,
	const int window_size, const float alpha, const float beta)
{
	const int offset = (window_size - 1) / 2;
	for (int i = 0; i < channels; i++)
		scale_vec[i] = data[i] * data[i];

	float scale_result = 0;
	for (int i = 0; i <= offset; i++)
		scale_result += scale_vec[i];
	target[0] = scale_result * alpha + 1.f;
	for (int i = 1; i < channels; i++)
	{
		if (i + offset < channels)
			scale_result += scale_vec[i + offset];
		if (i - offset > 0)
			scale_result -= scale_vec[i - offset - 1];
		target[i] = scale_result * alpha + 1.f;
	}

	//#if USE_BLAS
	//	vsPowx(channels, target, -beta, target);
	//	vsMul(channels, target, data, target);
	//#else
	for (int i = 0; i < channels; i++)
		target[i] = data[i] * ::pow(target[i], -beta);
	//#endif
}


inline void TransIm2Col(const float *input, const int input_channels, const int input_stride,
	const int filter_width, const int filter_height, const int filter_stride_x, const int filter_stride_y,
	float *output, const int output_width, const int output_height, const int output_stride)
{
	if (filter_width == 1 && filter_height == 1 && filter_stride_x == 1 && filter_stride_y == 1)
	{
		for (int h = 0; h < output_height; h++)
		{
			memcpy(output, input, output_width*input_channels*sizeof(float));
			output += output_stride;
			input += input_stride;
		}
	}
	else
	{
		for (int h = 0; h < output_height; h++)
		{
			float *dst = output;
			for (int w = 0; w < output_width; w++)
			{
				const float *src = input + w * filter_stride_x * input_channels;
				for (int i = 0; i < filter_height; i++)
				{
					memcpy(dst, src, filter_width*input_channels*sizeof(float));
					src += input_stride;
					dst += filter_width*input_channels;
				}
			}
			output += output_stride;
			input += filter_stride_y * input_stride;
		}
	}
}

inline void TransIm2ColSkip(const float *input, const int input_channels, const int input_stride,
	const int filter_width, const int filter_height, const int filter_stride_x, const int filter_stride_y, const int filter_skip_x, const int filter_skip_y,
	float *output, const int output_width, const int output_height, const int output_stride)
{
	if (filter_width == 1 && filter_height == 1 && filter_stride_x == 1 && filter_stride_y == 1 && filter_skip_x == 1 && filter_skip_y == 1)
	{
		for (int h = 0; h < output_height; h++)
		{
			memcpy(output, input, output_width*input_channels*sizeof(float));
			output += output_stride;
			input += input_stride;
		}
	}
	else if (filter_skip_x == 1 && filter_skip_y == 1)
	{
		for (int h = 0; h < output_height; h++)
		{
			float *dst = output;
			for (int w = 0; w < output_width; w++)
			{
				const float *src = input + w * filter_stride_x * input_channels;
				for (int i = 0; i < filter_height; i++)
				{
					memcpy(dst, src, filter_width*input_channels*sizeof(float));
					src += input_stride;
					dst += filter_width*input_channels;
				}
			}
			output += output_stride;
			input += filter_stride_y * input_stride;
		}
	}
	else
	{
		for (int h = 0; h < output_height; h++)
		{
			float *dst = output;
			for (int w = 0; w < output_width; w++)
			{
				for (int i = 0; i < filter_height; i++)
				{
					const float *src = input + w * filter_stride_x * input_channels + i * filter_skip_y * input_stride;

					for (int j = 0; j < filter_width; ++j)
					{
						memcpy(dst, src, input_channels*sizeof(float));
						src += filter_skip_x * input_channels;
						dst += input_channels;
					}
				}
			}
			output += output_stride;
			input += filter_stride_y * input_stride;
		}
	}

}