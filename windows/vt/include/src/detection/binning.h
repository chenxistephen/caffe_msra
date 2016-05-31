#pragma once

#include "vtcore.h"
#include "dhelpers.h"

namespace vt_internal
{
	///////////////////////////////////////////////////////////////////////////
	// angle binning

	union InterpolatedAngleBin
	{	
		InterpolatedAngleBin()
			: i(0)
		{
		}

		InterpolatedAngleBin(int init)
			: i(init)
		{
		}

		unsigned int i;
		struct
		{
			unsigned char bin0;
			unsigned char w0;
			// short lookups are slower than int -> pad to 32bit
			unsigned short unused; 
		} v;
	};

	template <unsigned int NumWeightFixShiftBits>
	class AngleBinningFunc
	{
	public:
		AngleBinningFunc(unsigned int numAngBits, unsigned int numBins)
			: numAngleBits(numAngBits), nBins(numBins)
		{
		}

		int operator()(unsigned int a)
		{	
			InterpolatedAngleBin tmp;

			const double fbin = nBins * a / float(1 << numAngleBits) * .999f;
			unsigned int bin0 = static_cast<int>(floor(fbin));
			const double alpha = 1 - (fbin - bin0);

			typedef Int32Clamper<false, NumWeightFixShiftBits> clamper;

			const unsigned int maxVal = 
				MaxInt32Value<false, NumWeightFixShiftBits>::value;

			unsigned int w0 = clamper::Clamp(
				F2I(static_cast<float>(alpha * maxVal)));

			InterpolatedAngleBin res;

			VT_ASSERT(bin0 >= 0 && bin0 < nBins);
			res.v.bin0 = static_cast<unsigned char>(
				VtClamp<unsigned int>(bin0, 0, nBins - 1));

			res.v.w0 = static_cast<unsigned char>(w0);
			VT_ASSERT(res.v.w0 == w0);

			return res.i;
		}

	private:
		unsigned int numAngleBits;
		unsigned int nBins;
	};

	template <typename StaticParams>
	class BinningEngine
	{
	public:
		BinningEngine()
			: isInitialized(false),
			numBins(0)
		{
		}

	public:
		HRESULT Initialize(unsigned int numAngBits, 
			unsigned int numOrientBins)
		{
			VT_HR_BEGIN();

			isInitialized = false;

			VT_HR_EXIT(numOrientBins >= minNumBins &&
				numOrientBins <= maxNumBins ? S_OK : E_NOTIMPL);

			numBins = numOrientBins;

			orientLUT.Initialize(AngleBinningFunc<
				StaticParams::orientBinWeightFixShift>(numAngBits, numBins));

			isInitialized = true;

			VT_HR_END();
		}

		template <typename GradMagNormalizerType>
		HRESULT BinMagAngToChannels(vt::CLumaByteImg& magChannel, 
			vt::CPlanarByteImg& orientChannels, const vt::CLumaByteImg& mag,
			const vt::CLumaByteImg& ang, const vt::CLumaByteImg& blurredMag,
			const GradMagNormalizerType& nrm) const
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(isInitialized ? S_OK : E_NOINIT);	

			if (!magChannel.IsValid() || !orientChannels.IsValid())
				return E_INVALIDDST;
			if (!mag.IsValid() || !ang.IsValid() || !blurredMag.IsValid())
				return E_INVALIDSRC;

			if (!AreSameSize(magChannel, orientChannels) || 
				!AreSameSize(mag, ang) || !AreSameSize(mag, blurredMag))
				return E_INVALIDARG;

			const int binSize = 1 << StaticParams::log2OrientBinSize;
			if (mag.Width() != magChannel.Width() * binSize ||
				mag.Height() != magChannel.Height() * binSize)
				return E_INVALIDARG;

#define BIN_FOR_NUM_BINS(nb) \
	static_assert(nb >= minNumBins && nb <= maxNumBins, \
	"invalid number of bins"); \
	const BinSpanFunc<nb, StaticParams::log2OrientBinSize> \
	binSpanFunc_##nb; \
	for (int y = 0; y < mag.Height(); y += binSize) \
			{ \
			binSpanFunc_##nb(magChannel, orientChannels, mag, ang, \
			blurredMag, nrm, orientLUT, y); \
			} \

			switch (numBins)
			{
			case 6: BIN_FOR_NUM_BINS(6); break;
			default:
				VT_ASSERT(false);
				VT_HR_EXIT(E_NOTIMPL);
			}

#undef BIN_FOR_NUM_BINS

			VT_HR_END();
		}

	private:
		template <unsigned int NumOrientBins, 
			unsigned int Log2BinSize>
		struct BinSpanFunc
		{
			template <typename GradMagNormalizerType, typename OrientLUTType>
			void operator()(vt::CLumaByteImg& magChannel, 
				vt::CPlanarByteImg& orientChannels, 
				const vt::CLumaByteImg& mag, const vt::CLumaByteImg& ang, 
				const vt::CLumaByteImg& blurredMag,
				const GradMagNormalizerType& nrm, 
				const OrientLUTType& orientLUT, unsigned int y) const
			{
				const unsigned int binSize = 1 << Log2BinSize;
				const unsigned int dstY = y >> Log2BinSize;

				unsigned int dstX = 0;			
				for (int x = 0; x < mag.Width(); x += binSize)
				{
					int magAcc = 0;
					int angAcc[NumOrientBins];
					for (unsigned int b = 0; b < NumOrientBins; ++b)
						angAcc[b] = 0;

					for (unsigned int v = 0; v < binSize; ++v)
					{
						const unsigned char* magSrcY = mag.Ptr(y + v);
						const unsigned char* blurredMagSrcY = 
							blurredMag.Ptr(y + v);
						const unsigned char* angSrcY = ang.Ptr(y + v);

						for (int u = 0; u < binSize; ++u)
						{
							auto m = nrm.Normalize(magSrcY[x+u], 
								blurredMagSrcY[x+u]);

							magAcc += m;

							const InterpolatedAngleBin 
								ibin(orientLUT.LookUp(angSrcY[x+u]));

							unsigned int bin1 = ibin.v.bin0 + 1;
							if (bin1 == NumOrientBins)
								bin1 = 0;

							const auto w1 = MaxInt32Value<false, 
								StaticParams::orientBinWeightFixShift>::value 
								- ibin.v.w0;

							angAcc[ibin.v.bin0] += m * ibin.v.w0;
							angAcc[bin1] += m * w1;
						}
					}

					typedef Int32Clamper<false, 8> clamper;

					magChannel(dstX, dstY) = static_cast<unsigned char>(
						clamper::RightShiftAndClamp<2 * Log2BinSize>(magAcc));

					for (unsigned int b = 0; b < NumOrientBins; ++b)
					{
						orientChannels.Band(b)(dstX, dstY) = 
							static_cast<unsigned char>(
							clamper::RightShiftAndClamp<(2 * Log2BinSize
							+ StaticParams::orientBinWeightFixShift)>(
							angAcc[b]));
					}

					++dstX;
				}
			}
		};

#ifdef _M_ARM
		template <unsigned int NumOrientBins, 
			unsigned int OrientBinWeightFixShift, 
			unsigned int NumLookUpBins>
		struct NEONOrientLUT;

		template <>
		struct NEONOrientLUT<6, 4, 32>
		{			
			static const unsigned int numLookUpBins = 32;

			void Initialize(const int* lutData)
			{
				unsigned char buf[numLookUpBins];

				for (unsigned int i = 0; i < numLookUpBins; ++i)
				{
					InterpolatedAngleBin iab(lutData[i]);
					const auto b = iab.v.bin0;
					const auto w = iab.v.w0;
					VT_ASSERT(b < 16);
					VT_ASSERT(w < 16);
					buf[i] = (w << 4) | b;
				}

				lut.Initialize(buf);
			}

			void LookUp(uint8x8_t& bin0, uint8x8_t& w0, uint8x8_t ang)
			{
				const auto v = lut.LookUp(ang);
				const auto binmask = vdup_n_u8(0xf0);
				bin0 = vbic_u8(v, binmask);
				w0 = vshr_n_u8(v, 4);
			}

			NEON1dUInt8LUT32 lut;
		};

		template <>
		struct BinSpanFunc<6, 2>
		{
			static const unsigned int NumOrientBins = 6;
			static const unsigned int Log2BinSize = 2;

			template <typename GradMagNormalizerType, typename OrientLUTType>
			void operator()(vt::CLumaByteImg& magChannel, 
				vt::CPlanarByteImg& orientChannels, 
				const vt::CLumaByteImg& mag, const vt::CLumaByteImg& ang, 
				const vt::CLumaByteImg& blurredMag,
				const GradMagNormalizerType& nrm, 
				const OrientLUTType& orientLUT, unsigned int y) const
			{
				const unsigned int binSize = 1 << Log2BinSize;
				const unsigned int dstY = y >> Log2BinSize;

				// check 16-stride since we don't handle boundaries
				VT_ASSERT(mag.StrideBytes() % 16 == 0);
				VT_ASSERT(ang.StrideBytes() % 16 == 0);
				VT_ASSERT(magChannel.StrideBytes() % 16 == 0);
				VT_ASSERT(orientChannels.StrideBytes() % 16 == 0);

				GradMagNormalizerType::NormalizeBatchTmpType nrmTmp;
				nrm.NormalizeBatchInit(nrmTmp);

				NEONOrientLUT<NumOrientBins, 
					StaticParams::orientBinWeightFixShift, 
					OrientLUTType::numBins> binLut;				
				binLut.Initialize(orientLUT.GetTablePtr());

				unsigned int dstX = 0;			
				for (int x = 0; x < mag.Width(); x += 8, dstX += 2)
				{
#define DECLARE_BIN(binnr)\
	auto bacc##binnr = vdupq_n_u32(0);

#define ACCUMULATE_BIN(binnr) \
	bacc##binnr = vmlal_u8(bacc##binnr, m, vbit_u8(vdup_n_u8(0), \
	w0, vceq_u8(vdup_n_u8(binnr), bin0))); \
	bacc##binnr = vmlal_u8(bacc##binnr, m, vbit_u8(vdup_n_u8(0), \
	w1, vceq_u8(vdup_n_u8(binnr), bin1)));

#define PADD_STORE_BIN(binnr) \
	const auto bres##binnr = vsriq_n_u64(vdupq_n_u32(0), \
	vpaddlq_u32(vpaddlq_u16(bacc##binnr)), \
	2 * Log2BinSize + StaticParams::orientBinWeightFixShift); \
	orientChannels.Band(binnr)(dstX, dstY) = bres##binnr.n128_u8[0]; \
	orientChannels.Band(binnr)(dstX + 1, dstY) = bres##binnr.n128_u8[8];

					auto macc = vdupq_n_u32(0); 
					DECLARE_BIN(0);
					DECLARE_BIN(1);
					DECLARE_BIN(2);
					DECLARE_BIN(3);
					DECLARE_BIN(4);
					DECLARE_BIN(5);

					for (int dy = 0; dy < binSize; ++dy)
					{
						const auto rm = vld1_u8(mag.Ptr(x, y + dy));
						const auto bm = vld1_u8(blurredMag.Ptr(x, y + dy));				
						const auto m = nrm.NormalizeBatch8(rm, bm, nrmTmp);

						const auto a = vld1_u8(ang.Ptr(x, y + dy));
						uint8x8_t bin0, w0;
						binLut.LookUp(bin0, w0, a);
						const auto bin1tmp = vadd_u8(bin0, vdup_n_u8(1));
						const auto bin1 = vbsl_u8(vclt_u8(bin1tmp, 
							vdup_n_u8(NumOrientBins)), bin1tmp, vdup_n_u8(0)); 
						const auto w1 = vsub_u8(vdup_n_u8(16), w0);

						macc = vaddw_u8(macc, m);
						ACCUMULATE_BIN(0);
						ACCUMULATE_BIN(1);
						ACCUMULATE_BIN(2);
						ACCUMULATE_BIN(3);
						ACCUMULATE_BIN(4);
						ACCUMULATE_BIN(5);						
					}

					const auto mres = vsriq_n_u64(vdupq_n_u32(0), 
						vpaddlq_u32(vpaddlq_u16(macc)), 2 * Log2BinSize);

					magChannel(dstX, dstY) = mres.n128_u8[0];
					magChannel(dstX + 1, dstY) = mres.n128_u8[8];

					PADD_STORE_BIN(0);					
					PADD_STORE_BIN(1);					
					PADD_STORE_BIN(2);					
					PADD_STORE_BIN(3);					
					PADD_STORE_BIN(4);					
					PADD_STORE_BIN(5);					

#undef PADD_STORE_BIN
#undef ACCUMULATE_BIN
#undef DECLARE_BIN	

				}
			}
		};
#endif

	private:
		bool isInitialized;

		unsigned int numBins;

		static const unsigned int minNumBins = 4;
		static const unsigned int maxNumBins = 8;
		static const unsigned int maxLog2BinSize = 4;

		static_assert(StaticParams::log2OrientBinSize <= maxLog2BinSize, 
			"bin size too large");

		C1DLookUp<false,
			StaticParams::numOrientBinLUTInputBits,
			StaticParams::log2NumOrientBinLUTBins, 
			typename StaticParams::OrientLUTBinType> orientLUT;		
	};
}