#include "stdafx.h"

#include "dhelpers.h"
#include <algorithm>

using namespace vt;

namespace vt_internal
{
	HRESULT BlurImage(CLumaByteImg& dst, const CLumaByteImg& src, 
		CLumaByteImg& buf, int numIterations)
	{
		if (numIterations == 0)
			return src.CopyTo(dst);
		
		return VtSeparableFilter121(dst, dst.Rect(), src, CPoint(0, 0), true,
			numIterations, &buf);
	}

	HRESULT BilinearResizeImage(CImg& imgDst, const CRect& rctDst, 
		const CImg& imgSrc)
	{
		auto sampler = eSamplerKernelBilinear;

#ifdef _M_ARM
		// At the time of this writing, VtWarpImage() is implemented in 
		// NEON already, but VtResizeImage() is not

		CMtx3x3f scale;
		const float sx = float(imgSrc.Width()) / imgDst.Width();
		const float sy = float(imgSrc.Height()) / imgDst.Height();
		scale.MakeScale(sx, sy);
		CMtx3x3f shift;
		shift.MakeTranslate(.5f, .5f);	
		auto xfrm = shift.Inv()*scale*shift;

		return VtWarpImage(imgDst, rctDst, imgSrc, xfrm, sampler);
#else
		return VtResizeImage(imgDst, rctDst, imgSrc, sampler, 
			IMAGE_EXTEND(Extend), eResizeSamplingSchemeDual);
#endif
	}

	template<typename SrcElT, typename DstElT>
	struct ByteScaleOp
	{
	public:
		typedef unsigned char ImplSrcElType;
		typedef unsigned char ImplDstElType;

		static const unsigned int scaleShiftBits = 6;

		struct ParamType
		{
			unsigned char shiftedScale;			
		};

		struct TmpType
		{
#ifdef _M_ARM
			uint8x8_t factor;
#endif
		};

		static const int NumSrcBands = 0;
		static const int NumDstBands = 0;
		static const int NumSrcElPerGroup = 1;
		static const int NumDstElPerGroup = 1;
		static const int NumGroupsPerOpGeneric = 1;
		static const int NumGroupsPerOpNEON = 16;

	public:
		void EvalGeneric(const ImplSrcElType* pSrc, ImplDstElType* pDst, 
			ParamType* param, TmpType&)
		{
			const auto tmp = pSrc[0] * (*param).shiftedScale;

			pDst[0] = static_cast<unsigned char>(
				vt_internal::Int32Clamper<false, 8>
				::RightShiftAndClamp<scaleShiftBits>(tmp));
		}

#ifdef _M_ARM
		void InitNEON(ParamType* param, TmpType& tmp)
		{
			tmp.factor = vdup_n_u8((*param).shiftedScale);			
		}

		void EvalNEON(const ImplSrcElType* pSrc, ImplDstElType* pDst, 
			ParamType*, const TmpType& tmp)
		{
			auto u8_0 = vld1_u8(pSrc);
			auto u8_1 = vld1_u8(&pSrc[8]);
			auto prod_0 = vmull_u8(u8_0, tmp.factor);
			auto prod_1 = vmull_u8(u8_1, tmp.factor);
			auto res_0 = vqshrn_n_u16(prod_0, scaleShiftBits);
			auto res_1 = vqshrn_n_u16(prod_1, scaleShiftBits);
			vst1_u8(pDst, res_0);
			vst1_u8(&pDst[8], res_1);
		}
#endif
	};

	HRESULT ScaleByteImage(CLumaByteImg& dst, const CLumaByteImg& src, 
		float factor)
	{
		VT_HR_BEGIN();

		if (!src.IsValid())
			return E_INVALIDSRC;
		if (!dst.IsValid())
			return E_INVALIDDST;
		if (!AreSameSize(dst, src))
			return E_INVALIDARG;

#ifdef _M_ARM
		typedef ByteScaleOp<unsigned char, unsigned char> OT;
		OT::ParamType p;
		auto iscale = F2I(factor * (1 << OT::scaleShiftBits));
		p.shiftedScale = static_cast<unsigned char>(iscale);

		if (p.shiftedScale != iscale)
			return E_INVALIDARG;

		VT_HR_EXIT((
			UnaryImgOpSS<ByteScaleOp, unsigned char, unsigned char>(
			src, dst, &p)));
#else
		VT_HR_EXIT(VtScaleImage(dst, src, factor));
#endif		

		VT_HR_END();
	}	
}
