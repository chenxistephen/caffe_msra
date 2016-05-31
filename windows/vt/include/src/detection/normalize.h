#pragma once

#include "vtcore.h"
#include "dhelpers.h"

namespace vt_internal
{
	template <typename OutputType>
	class RegularizedNormFactorFunc
	{	
	public:
		RegularizedNormFactorFunc(unsigned int numPreShiftBits,
			unsigned int regularization, float scaling, 
			unsigned int numFixShiftBits)
			: preShift(numPreShiftBits), 
			reg(regularization), 
			scale(scaling), 
			shift(numFixShiftBits)
		{
		}

	public:
		OutputType operator()(unsigned int x)
		{
			double f = scale / (x * (1 << preShift) + reg);

			const unsigned int intf = 
				F2I(static_cast<float>(f * (1 << shift)));

			//const OutputType res = static_cast<OutputType>(intf);
			//VT_ASSERT(res == intf);			
			return static_cast<OutputType>(vt::VtClamp<unsigned int>(intf, 
				0, vt::ElTraits<OutputType>::MaxVal()));
		}

	private:
		unsigned int preShift;
		unsigned int reg;
		float scale;
		unsigned int shift;
	};

	template <typename StaticParams>
	class NormalizationEngine
	{
	public:
		NormalizationEngine()
			: isInitialized(false)
		{
		}

	public:
		HRESULT Initialize(unsigned int regularization, float scaling)
		{
			isInitialized = false;

			factorLUT.Initialize(
				RegularizedNormFactorFunc<
				StaticParams::NormFactorLUTBinType>(
				StaticParams::normFactorLUTPreShift,
				regularization, scaling,			
				StaticParams::normFactorFixShift));

			isInitialized = true;

			return S_OK;
		}

		unsigned int Normalize(unsigned int input, unsigned int blurredInput) 
			const
		{
			VT_ASSERT(isInitialized);

			VT_ASSERT((UsesAtMostNumBits<false, 8>(blurredInput)));

			typedef Int32Clamper<false, 
				StaticParams::numNormFactorLUTInputBits + 
				StaticParams::normFactorLUTPreShift, 8> LUTInputClamper;

			const unsigned int index =
				LUTInputClamper::RightShiftAndClamp<
				StaticParams::normFactorLUTPreShift>(blurredInput);

			const unsigned int factor = factorLUT.LookUp(index);

			return (input * factor) >>
				(StaticParams::normFactorFixShift - 8);
		}

#ifdef _M_ARM
		typedef NEON1dUInt8LUT64 NormalizeBatchTmpType;

		HRESULT NormalizeBatchInit(NormalizeBatchTmpType& tmp) const
		{
			if (!isInitialized)
				return E_NOINIT;

			tmp.Initialize(factorLUT.GetTablePtr());

			return S_OK;
		}

		static uint8x8_t NormalizeBatch8(uint8x8_t input, 
			uint8x8_t blurredInput, const NormalizeBatchTmpType& tmp)
		{
			const auto bl = vshr_n_u8(blurredInput, 
				StaticParams::normFactorLUTPreShift);
			const auto fac = tmp.LookUp(bl);
			const auto prod = vmull_u8(input, fac);
			return vqshrn_n_u16(prod, StaticParams::normFactorFixShift - 8);
		}
		
		static uint8x16_t NormalizeBatch16(uint8x16_t input, 
			uint8x16_t blurredInput, const NormalizeBatchTmpType& tmp)
		{
			const auto bl = vshrq_n_u8(blurredInput, 
				StaticParams::normFactorLUTPreShift);
			const auto fac = tmp.LookUp(bl);
			const auto prod_l = 
				vmull_u8(vget_low_u8(input), vget_low_u8(fac));
			const auto prod_h = 
				vmull_u8(vget_high_u8(input), vget_high_u8(fac));
			const auto res_l = vqshrn_n_u16(prod_l, 
				StaticParams::normFactorFixShift - 8);
			const auto res_h = vqshrn_n_u16(prod_h, 
				StaticParams::normFactorFixShift - 8);

			return vcombine_u8(res_l, res_h);
		}
#endif

		HRESULT NormalizeImage(vt::CLumaByteImg& blurredInputAndOutput, 
			const vt::CLumaByteImg& input) const
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(isInitialized ? S_OK : E_NOINIT);

			if (!input.IsValid())
				return E_INVALIDSRC;

			if (!AreSameSize(input, blurredInputAndOutput))
				return E_INVALIDARG;

			for (int y = 0; y < input.Height(); ++y)
			{
				NormalizeSpan<StaticParams::log2NumNormFactorLUTBins>
					(blurredInputAndOutput.Ptr(y), 
					blurredInputAndOutput.StrideBytes(), input.Ptr(y), 
					input.StrideBytes(), input.Width());
			}

			VT_HR_END();
		}

	private:
		template <unsigned int log2LUTSize>
		void NormalizeSpan(unsigned char* blurredInputAndOutput, 
			unsigned int, const unsigned char* input, unsigned int, 
			unsigned int length) const
		{
			VT_ASSERT(isInitialized);

			typedef Int32Clamper<false, 8> clamper;

			for (unsigned int x = 0; x < length; ++x)
			{
				blurredInputAndOutput[x] = static_cast<unsigned char>(
					clamper::Clamp(Normalize(input[x], 
					blurredInputAndOutput[x])));
			}
		}

#if _M_ARM
		template <>
		void NormalizeSpan<6>(unsigned char* blurredInputAndOutput, 
			unsigned int blurredInputStrideBytes, const unsigned char* input, 
			unsigned int inputStrideBytes, unsigned int length) const
		{
			VT_ASSERT(isInitialized);			
			
			// we don't handle boundaries, so check for 16byte stride
			VT_ASSERT(length <= blurredInputStrideBytes);
			VT_ASSERT(length <= inputStrideBytes);
			VT_ASSERT(blurredInputStrideBytes % 16 == 0);
			VT_ASSERT(inputStrideBytes % 16 == 0);

			// only used for asserts
			UNREFERENCED_PARAMETER(blurredInputStrideBytes);
			UNREFERENCED_PARAMETER(inputStrideBytes);
			
			NEON1dUInt8LUT64 lut(factorLUT.GetTablePtr());

			for (unsigned int x = 0; x < length; x += 8)
			{
				const auto in16 = vld1_u8(&input[x]);
				const auto bl16 = vld1_u8(&blurredInputAndOutput[x]);

				vst1_u8(&blurredInputAndOutput[x], 
					NormalizeBatch8(in16, bl16, lut));
			}
		}
#endif

	private:
		bool isInitialized;

		C1DLookUp<false,
			StaticParams::numNormFactorLUTInputBits,
			StaticParams::log2NumNormFactorLUTBins, 
			typename StaticParams::NormFactorLUTBinType> factorLUT;
	};
}