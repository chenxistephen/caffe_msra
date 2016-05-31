#pragma once

namespace vt_internal
{
	///////////////////////////////////////////////////////////////////////////
	// integer clamping

	template <bool IsSigned>
	struct Int32Type
	{
		typedef typename 
			std::conditional<IsSigned, int, unsigned int>::type type;
	};

	template <bool IsSigned, unsigned int NumBits>
	struct MinInt32Value;

	template <>
	struct MinInt32Value<true, 0>
	{
		static const int value = -1;
	};

	template <unsigned int NumBits>
	struct MinInt32Value<true, NumBits>
	{
		static const int value = -(1 << (NumBits - 1));
	};

	template <unsigned int NumBits>
	struct MinInt32Value<false, NumBits>
	{
		static const unsigned int value = 0;
	};

	template <bool IsSigned, unsigned int NumBits>
	struct MaxInt32Value;

	template <unsigned int NumBits>
	struct MaxInt32Value<true, NumBits>
	{
		static const int value = (1i64 << (NumBits - 1)) - 1;
	};

	template <>
	struct MaxInt32Value<true, 0>
	{
		static const int value = 0;
	};

	template <unsigned int NumBits>
	struct MaxInt32Value<false, NumBits>
	{
		static const unsigned int value = (1i64 << NumBits) - 1;
	};

	template <bool IsSigned, unsigned int NumBits>
	inline bool UsesAtMostNumBits(typename Int32Type<IsSigned>::type v)
	{
		return 	v >= MinInt32Value<IsSigned, NumBits>::value &&
			v <= MaxInt32Value<IsSigned, NumBits>::value;
	}

	template <bool IsSigned,
		unsigned int NumOutputBits,
		unsigned int NumInputBits = 32>
	struct Int32Clamper
	{
	public:
		typedef typename Int32Type<IsSigned>::type IntType;

		static IntType Clamp(IntType input)
		{	
			return ClampInternal(input,
				std::conditional<needsToClamp, std::true_type, 
				std::false_type>::type());
		}

		template <unsigned int NumShiftBits>
		static IntType RightShiftAndClamp(IntType input)
		{	
			return RightShiftAndClampInternal<NumShiftBits>(input,
				std::conditional<needsToClamp, std::true_type, 
				std::false_type>::type());
		}

	private:
		static const bool needsToClamp = 
			NumInputBits > NumOutputBits;
		static const IntType minOutputValue = 
			MinInt32Value<IsSigned, NumOutputBits>::value;
		static const IntType maxOutputValue = 
			MaxInt32Value<IsSigned, NumOutputBits>::value;

		static IntType ClampInternal(IntType input, const std::true_type&)
		{
			VT_ASSERT((UsesAtMostNumBits<IsSigned, 
				NumInputBits>(input)));

			// Note: don't use the usual nested std::min/max here: in case of 
			// IntType = unsigned int, the compiler will generate an 
			// unnecessary test against > 0 .. haven't figured out why.
			const IntType tmp = input < minOutputValue ? 
				minOutputValue : input;
			return tmp > maxOutputValue ? maxOutputValue : tmp;
		}

		static IntType ClampInternal(IntType input, const std::false_type&)
		{
			VT_ASSERT((UsesAtMostNumBits<IsSigned, 
				NumInputBits>(input)));

			return input;
		}

		template <unsigned int NumShiftBits, typename ClampNeededType>
		static IntType RightShiftAndClampInternal(IntType input, 
			const ClampNeededType&)
		{
			return ClampInternal(input >> NumShiftBits, ClampNeededType());
		}
	};

	///////////////////////////////////////////////////////////////////////////
	// lookups

	template <bool IsIntSigned, unsigned int NumInputBits, 
		unsigned int Log2NumBins>
	struct IntBinConverter
	{
		static_assert(Log2NumBins <= NumInputBits, "too many bins");

		typedef typename Int32Type<IsIntSigned>::type IntType;

		static const unsigned int offset = 
			IsIntSigned ? ((1 << Log2NumBins) / 2) : 0;
		static const unsigned int shift = 
			NumInputBits - Log2NumBins;

		static unsigned int IntToBin(IntType v)
		{
			VT_ASSERT((UsesAtMostNumBits<IsIntSigned, NumInputBits>(v)));
			return (v >> shift) + offset;
		}

		static IntType BinToInt(unsigned int bin)
		{
			VT_ASSERT((UsesAtMostNumBits<false, Log2NumBins>(bin)));
			return (bin - offset) << shift;
		}
	};

	template <bool IsInputSignedX, bool IsInputSignedY, 
		unsigned int NumInputBitsX, unsigned int NumInputBitsY, 
		unsigned int Log2NumBinsX, unsigned int Log2NumBinsY,
		typename BinType>
	struct C2DLookUp
	{
	public:
		static const unsigned int numInputBitsX = NumInputBitsX;
		static const unsigned int numInputBitsY = NumInputBitsY;
		static const unsigned int numBinsX = (1 << Log2NumBinsX);
		static const unsigned int numBinsY = (1 << Log2NumBinsY);

	private:
		typedef typename Int32Type<IsInputSignedX>::type InputTypeX;
		typedef IntBinConverter<IsInputSignedX, NumInputBitsX, 
			Log2NumBinsX> I2BConverterX;

		typedef typename Int32Type<IsInputSignedY>::type InputTypeY;
		typedef IntBinConverter<IsInputSignedY, NumInputBitsY, 
			Log2NumBinsY> I2BConverterY;

	public:
		C2DLookUp()
		{
#ifdef _DEBUG
			lutAsCImg.Create(reinterpret_cast<vt::Byte*>(lut), numBinsX, 
				numBinsY, numBinsX * sizeof(BinType));
			hitCount.Create(numBinsX, numBinsY);
			hitCount.Fill(-1.f);
#ifndef _M_ARM
			VtImageDebuggerSetImageName(lutAsCImg, "LUT");
			VtImageDebuggerSetImageName(hitCount, "LUTHitCount");
#endif
#endif
		}

		template <typename BinaryFuncType>
		void Initialize(BinaryFuncType func)
		{
			for (unsigned int xBin = 0; xBin < numBinsX; ++xBin)
			{
				for (unsigned int yBin = 0; yBin < numBinsY; ++yBin)
				{
					InitializeEntry(xBin, yBin, 
						func(I2BConverterX::BinToInt(xBin), 
						I2BConverterY::BinToInt(yBin)));
				}
			}
		}

		const BinType LookUp(InputTypeX x, InputTypeY y) const
		{
			VT_ASSERT((UsesAtMostNumBits<IsInputSignedX, 
				NumInputBitsX>(x)));
			VT_ASSERT((UsesAtMostNumBits<IsInputSignedY, 
				NumInputBitsY>(y)));

			const auto xBin = I2BConverterX::IntToBin(x);
			const auto yBin = I2BConverterY::IntToBin(y);

#ifdef _DEBUG
			hitCount.Pix(xBin, yBin) += 1;
#endif

			return LookUpInternal(xBin, yBin);
		}

		const BinType LookUp(InputTypeX x) const
		{
			return LookUp(x, 0);
		}

		const BinType* GetTablePtr() const
		{
			return &lut[0];
		}		

	protected:
		template <typename T>
		void InitializeEntry(unsigned int xBin, unsigned int yBin, T v)
		{
			const auto cast = static_cast<BinType>(v);
			VT_ASSERT(cast == v);

			LookUpInternal(xBin, yBin) = cast;

#ifdef _DEBUG
			hitCount.Pix(xBin, yBin) = 0;
#endif
		}

	private:
		const BinType& LookUpInternal(unsigned int xBin, 
			unsigned int yBin) const
		{
			VT_ASSERT(xBin < numBinsX && yBin < numBinsY);
			const unsigned int index = (yBin << Log2NumBinsX) + xBin;

			VT_ASSERT(index < sizeof(lut)/sizeof(lut[0]));
			return lut[index];
		}

		BinType& LookUpInternal(unsigned int xBin, unsigned int yBin)
		{
			typedef C2DLookUp<IsInputSignedX, IsInputSignedY, NumInputBitsX, 
				NumInputBitsY, Log2NumBinsX, Log2NumBinsY, BinType> MyType;

			return const_cast<BinType&>(static_cast<const MyType*>(this)
				->LookUpInternal(xBin, yBin));
		}

	private:		
		BinType lut[numBinsX * numBinsY];

#ifdef _DEBUG
		vt::CCompositeImg<vt::LumaType<BinType>> lutAsCImg;
		mutable vt::CCompositeImg<vt::LumaType<float>, float> hitCount;
#endif
	};

	template <bool IsInputSigned, unsigned int NumInputBits, 
		unsigned int Log2NumBins, typename BinType>
	struct C2DSquareLookUp 
		: public C2DLookUp<IsInputSigned, IsInputSigned,
		NumInputBits, NumInputBits, Log2NumBins, Log2NumBins, BinType>
	{
	};

	template <bool IsInputSigned, unsigned int NumInputBits, 
		unsigned int Log2NumBins, typename BinType>
	struct C1DLookUp 
		: public C2DLookUp<IsInputSigned, IsInputSigned,
		NumInputBits, NumInputBits, Log2NumBins, 0, BinType>
	{
		static const unsigned int numBins = numBinsX;

		template <typename UnaryFuncType>
		void Initialize(UnaryFuncType func)
		{
			for (unsigned int xBin = 0; xBin < numBinsX; ++xBin)
			{
				InitializeEntry(xBin, 0, 
					func(I2BConverterX::BinToInt(xBin)));				
			}
		}
	};

#ifdef _M_ARM
	class NEON1dUInt8LUT32
	{	
	public:
		NEON1dUInt8LUT32()
		{
		}

		NEON1dUInt8LUT32(const unsigned char *lutData)
		{
			Initialize(lutData);
		}

		void Initialize(const unsigned char* lutData)
		{
			VT_ASSERT(lutData != nullptr);

			lut.val[0] = vld1_u8(lutData);
			lut.val[1] = vld1_u8(lutData + 8);			
			lut.val[2] = vld1_u8(lutData + 16);
			lut.val[3] = vld1_u8(lutData + 24);			
		}

		uint8x8_t LookUp(const uint8x8_t& input) const
		{
			return vtbl4_u8(lut, input);
		}

		uint8x16_t LookUp(const uint8x16_t& input) const
		{
			const auto fac_0 = vtbl4_u8(lut, vget_low_u8(input));
			const auto fac_1 = vtbl4_u8(lut, vget_high_u8(input));
			return vcombine_u8(fac_0, fac_1);
		}

	private:		
		uint8x8x4_t lut;
	};
	
	class NEON1dUInt8LUT64
	{	
	public:
		NEON1dUInt8LUT64()
		{
		}

		NEON1dUInt8LUT64(const unsigned char *lutData)
		{
			Initialize(lutData);
		}

		void Initialize(const unsigned char* lutData)
		{
			VT_ASSERT(lutData != nullptr);

			lutLo.val[0] = vld1_u8(lutData);
			lutLo.val[1] = vld1_u8(lutData + 8);			
			lutLo.val[2] = vld1_u8(lutData + 16);
			lutLo.val[3] = vld1_u8(lutData + 24);
			lutHi.val[0] = vld1_u8(lutData + 32);
			lutHi.val[1] = vld1_u8(lutData + 40);
			lutHi.val[2] = vld1_u8(lutData + 48);
			lutHi.val[3] = vld1_u8(lutData + 56);
			hiMask = vdupq_n_u8(0xE0);
		}

		uint8x8_t LookUp(const uint8x8_t& input) const
		{
			const auto inIsHi = vtst_u8(input, vget_low_u8(hiMask));
			uint8x8_t inIdx; inIdx.n64_u64[0] = 0;
			inIdx = vbif_u8(inIdx, input, vget_low_u8(hiMask));
			const auto facLo = vtbl4_u8(lutLo, inIdx);
			const auto facHi = vtbl4_u8(lutHi, inIdx);
			return vbsl_u8(inIsHi, facHi, facLo);
		}

		uint8x16_t LookUp(const uint8x16_t& input) const
		{
			const auto inIsHi = vtstq_u8(input, hiMask);
			uint8x16_t inIdx; inIdx.n128_u64[0] = 0; inIdx.n128_u64[1] = 0;
			inIdx = vbifq_u8(inIdx, input, hiMask);
			const auto facLo_l = vtbl4_u8(lutLo, vget_low_u8(inIdx));
			const auto facHi_l = vtbl4_u8(lutHi, vget_low_u8(inIdx));
			const auto facLo_h = vtbl4_u8(lutLo, vget_high_u8(inIdx));
			const auto facHi_h = vtbl4_u8(lutHi, vget_high_u8(inIdx));

			return vbslq_u8(inIsHi, vcombine_u8(facHi_l, facHi_h), 
				vcombine_u8(facLo_l, facLo_h));
		}

	private:
		uint8x16_t hiMask;
		uint8x8x4_t lutLo, lutHi;
	};
#endif

	///////////////////////////////////////////////////////////////////////////
	// other helpers

	template <typename T, typename U>
	inline bool AreSameSize(const T& a, const U& b)
	{
		if (a.Width() != b.Width())
			return false;
		if (a.Height() != b.Height())
			return false;

		return true;
	}

	HRESULT BlurImage(vt::CLumaByteImg& dst, const vt::CLumaByteImg& src, 
		vt::CLumaByteImg& buf, int numIterations);

	HRESULT BilinearResizeImage(vt::CImg& imgDst, const vt::CRect& rctDst, 
		const vt::CImg& imgSrc);
	
	HRESULT ScaleByteImage(vt::CLumaByteImg& dst, const vt::CLumaByteImg& src, 
		float factor);
}