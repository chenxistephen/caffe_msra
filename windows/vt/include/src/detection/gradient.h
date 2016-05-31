#pragma once

#include "vtcore.h"
#include "dhelpers.h"
#include <algorithm>

namespace vt_internal
{
	///////////////////////////////////////////////////////////////////////////
	// 2d gradient computation

	enum PixelPosition
	{
		First,
		Inner,
		Last
	};

	template <PixelPosition PP>
	inline const unsigned char* MakeX0Ptr(const unsigned char* pos, 
		int stepSize);
	template <PixelPosition PP>
	inline const unsigned char* MakeX1Ptr(const unsigned char* pos, 
		int stepSize);

#define DEFINE_MAKE_PTR(ptype, x0delta, x1delta)\
	template <>\
	inline const unsigned char* MakeX0Ptr<ptype>(\
	const unsigned char* pos, int stepSize)\
	{ return pos + x0delta * stepSize; }\
	template <>\
	inline const unsigned char* MakeX1Ptr<ptype>(\
	const unsigned char* pos, int stepSize)\
	{ return pos + x1delta * stepSize; }

	DEFINE_MAKE_PTR(First, 0, 1);
	DEFINE_MAKE_PTR(Inner, -1, 1);
	DEFINE_MAKE_PTR(Last, -1, 0);

#undef DEFINE_MAKE_PTR

	template <PixelPosition PP>
	class DiffFunc
	{
	public:
		static int GetDiffInt8(unsigned char x0, unsigned char x1)
		{
			const int diff = x1 - x0;

			// difference generates one more bit than x0, x1
			typedef Int32Clamper<true, 8, 9> Clamper;

			return Clamper::Clamp(diff);
		}
	};

	template <>
	class DiffFunc<Inner>
	{
	public:
		static int GetDiffInt8(unsigned char x0, unsigned char x1)
		{
			return (x1 - x0) >> 1;
		}
	};

	template <PixelPosition XPosition, 
		PixelPosition YPosition>
	class DiffHelper
	{
	public:
		DiffHelper(const unsigned char *input, int numStrideBytes)
			: pos(0),
			left(MakeX0Ptr<XPosition>(input, 1)),
			right(MakeX1Ptr<XPosition>(input, 1)),
			top(MakeX0Ptr<YPosition>(input, numStrideBytes)),
			bottom(MakeX1Ptr<YPosition>(input, numStrideBytes))
		{
		}

		int GetXDiffInt8() const
		{
			typedef DiffFunc<XPosition> Differ;

			const int res = Differ::GetDiffInt8(left[pos], right[pos]);
			VT_ASSERT((UsesAtMostNumBits<true, 8>(res)));

			return res;
		}

		int GetYDiffInt8() const
		{
			typedef DiffFunc<YPosition> Differ;

			const int res = Differ::GetDiffInt8(top[pos], bottom[pos]);
			VT_ASSERT((UsesAtMostNumBits<true, 8>(res)));

			return res;
		}

		void Advance()
		{
			++pos;
		}

	private:
		size_t pos;

		const unsigned char* left;
		const unsigned char* top;
		const unsigned char* right;
		const unsigned char* bottom;
	};

	///////////////////////////////////////////////////////////////////////////
	// magnitude func

	template <unsigned int NumOutputBits>
	class ClampedMagnitudeFunc
	{
	public:
		ClampedMagnitudeFunc(float scale)
			: factor(scale)
		{
		}

		unsigned int operator()(int x, int y) const
		{
			return Int32Clamper<false, NumOutputBits>::Clamp(
				F2I(factor * sqrtf(static_cast<float>(x*x + y*y))));
		}

	private:		
		float factor;
	};

	///////////////////////////////////////////////////////////////////////////
	// angle func

	template <unsigned int NumAngleBits>
	class HalfAngleFunc
	{
	public:
		HalfAngleFunc(bool flipXY)
			: flip(flipXY)
		{
		}

		unsigned int operator()(int x, int y)
		{
			if (flip)
				std::swap(x, y);

			double a = atan2(static_cast<double>(y), 
				static_cast<double>(x));

			if (a < 0)
				a += VT_PI;	

			a *=  (1 << NumAngleBits) / (VT_PI + 1e-4);

			return static_cast<unsigned int>(floor(a));
		}

	private:
		bool flip;
	};

	///////////////////////////////////////////////////////////////////////////
	// gradient engine
		
	template <typename StaticParams>
	class GradientEngine
	{
	public:
		GradientEngine()
			: isInitialized(false)
		{
		}

	public:
		HRESULT Initialize(bool flipXY)
		{
			isInitialized = false;

			magLUT.Initialize(ClampedMagnitudeFunc<
				StaticParams::numMagBits>(1.f));
						
			angLUT.Initialize(HalfAngleFunc<StaticParams::numAngBits>(
				flipXY));

			isInitialized = true;

			return S_OK;
		}

		HRESULT ComputeMagAngImage(vt::CLumaByteImg& mag, 
			vt::CLumaByteImg& ang, const vt::CLumaByteImg& input) const
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(isInitialized ? S_OK : E_NOINIT);
			
			if (!mag.IsValid() || !AreSameSize(mag, input))
				return E_INVALIDDST;
			if (!ang.IsValid() || !AreSameSize(ang, input))
				return E_INVALIDDST;
			if (!input.IsValid())
				return E_INVALIDSRC;
			
			ComputeMagAngRow<First>(mag.Ptr(0), ang.Ptr(0), input.Ptr(0), 
				input.Width(), input.StrideBytes());

			for (int y = 1; y < input.Height() - 1; ++y)
			{
				ComputeMagAngRow<Inner>(mag.Ptr(y), ang.Ptr(y), input.Ptr(y), 
					input.Width(), input.StrideBytes());
			}

			ComputeMagAngRow<Last>(mag.Ptr(input.Height() - 1), 
				ang.Ptr(input.Height() - 1), input.Ptr(input.Height() - 1), 
				input.Width(), input.StrideBytes());

			VT_HR_END();
		}
		
		HRESULT CombineMagAngImages(vt::CLumaByteImg& mag, 
			vt::CLumaByteImg& ang, const vt::CPlanarByteImg& allMag,
			const vt::CPlanarByteImg& allAng)
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(isInitialized ? S_OK : E_NOINIT);
			
			if (!mag.IsValid() || !ang.IsValid())
				return E_INVALIDDST;

			if (!AreSameSize(mag, ang) ||
				!AreSameSize(mag, allMag) ||				
				!AreSameSize(ang, allAng))
				return E_INVALIDARG;

			if (!allMag.IsValid() || !allAng.IsValid())
				return E_INVALIDSRC;

			if (allMag.Bands() != allAng.Bands())
				return E_INVALIDARG;

#define NUM_BANDS_CASE(nb_) \
	case nb_: \
	VT_HR_EXIT(CombineMagAngImagesInternal<nb_>(mag, ang, allMag, allAng)); \
	break;

			switch (allMag.Bands())
			{
				NUM_BANDS_CASE(1);
				NUM_BANDS_CASE(2);
				NUM_BANDS_CASE(3);
				NUM_BANDS_CASE(4);
				NUM_BANDS_CASE(5);
			default:
				return E_INVALIDARG;
			};
	
#undef NUM_BANDS_CASE

			VT_HR_END();
		}

	private:
		template <PixelPosition XPosition, PixelPosition YPosition>
		void ComputeMagAngSpan(unsigned char* mag, unsigned char* ang,
			const unsigned char* input, int length, int numInputStrideBytes)
			const
		{
			DiffHelper<XPosition, YPosition> dh(input, numInputStrideBytes);

			for (int i = 0; i < length; ++i)
			{
				typedef Int32Clamper<true, 
					StaticParams::numMagLUTInputBits, 8> clamper;

				const int dx = clamper::Clamp(dh.GetXDiffInt8());
				const int dy = clamper::Clamp(dh.GetYDiffInt8());

				mag[i] = magLUT.LookUp(dx, dy);
				ang[i] = angLUT.LookUp(dx, dy);

				dh.Advance();
			}
		}

		template <PixelPosition YPosition>
		void ComputeMagAngRow(unsigned char* mag, unsigned char* ang, 
			const unsigned char* input, int width, int numInputStrideBytes)
			const
		{
			ComputeMagAngSpan<First, YPosition>(&mag[0], &ang[0], &input[0], 
				1, numInputStrideBytes);
			ComputeMagAngSpan<Inner, YPosition>(&mag[1], &ang[1], &input[1], 
				width - 2, numInputStrideBytes);
			ComputeMagAngSpan<Last, YPosition>(&mag[width - 1], 
				&ang[width - 1], &input[width - 1], 1, numInputStrideBytes);
		}

		template <unsigned int NumBands> 
		HRESULT CombineMagAngImagesInternal(
			vt::CLumaByteImg& mag, vt::CLumaByteImg& ang, 
			const vt::CPlanarByteImg& allMag, const vt::CPlanarByteImg& allAng)
		{
			for (int y = 0; y < mag.Height(); ++y)
			{
				const unsigned char* srcMagPtr[NumBands];
				for (unsigned int b = 0; b < NumBands; ++b)
					srcMagPtr[b] = allMag.Band(b).Ptr(y);
				
				const unsigned char* srcAngPtr[NumBands];
				for (unsigned int b = 0; b < NumBands; ++b)
					srcAngPtr[b] = allAng.Band(b).Ptr(y);

				unsigned char* dstMagPtr = mag.Ptr(y);
				unsigned char* dstAngPtr = ang.Ptr(y);

				for (int x = 0; x < mag.Width(); ++x)
				{
					unsigned char maxVal = 0;
					unsigned int maxIdx = 0;
					for (unsigned int b = 0; b < NumBands; ++b)
					{
						if (srcMagPtr[b][x] > maxVal)
						{
							maxVal = srcMagPtr[b][x];
							maxIdx = b;
						}
					}

					dstMagPtr[x] = maxVal;
					dstAngPtr[x] = srcAngPtr[maxIdx][x];
				}
			}

			return S_OK;
		}

		template <> 
		HRESULT CombineMagAngImagesInternal<1>(
			vt::CLumaByteImg& mag, vt::CLumaByteImg& ang, 
			const vt::CPlanarByteImg& allMag, const vt::CPlanarByteImg& allAng)
		{
			VT_HR_BEGIN();

			VT_HR_EXIT(allMag.Band(0).CopyTo(mag));

			VT_HR_EXIT(allAng.Band(0).CopyTo(ang));

			VT_HR_END();
		}

	private:
		bool isInitialized;

		C2DSquareLookUp<true, 
			StaticParams::numMagLUTInputBits, 
			StaticParams::log2NumMagLUTBins, 
			typename StaticParams::MagLUTBinType> magLUT;

		C2DSquareLookUp<true, 
			StaticParams::numAngLUTInputBits, 
			StaticParams::log2NumAngLUTBins, 
			typename StaticParams::MagLUTBinType> angLUT;
	};	
}