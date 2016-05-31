#include "stdafx.h"

#include "dhelpers.h"
#include "vt_rgb2luv.h"

using namespace vt;

namespace vt_internal
{
	struct LUVConversionConstants
	{
		LUVConversionConstants() 
			: unorm(0.197833f),
			vnorm(0.468331f),
			y0(0.00885645f),
			a(9.03296e+02f),
			umin(-0.325925f),
			vmin(-0.496296f),
			M(0.178325f, 0.341550f, 0.430574f, 
			0.071330f, 0.706655f, 0.222015f,
			0.939180f, 0.129553f, 0.020183f)
		{
			size_t i = 0;
			for (; i <= 1024; ++i)
			{
				const float y = i / 1024.f;
				const float f = 0.00370370f;
				if (y > y0)
					yToL[i] = f * (116 * powf(y, .333333f) - 16);
				else
					yToL[i] = f * a * y;
			}		
			for (; i < 1088; ++i)
				yToL[i] = yToL[i-1];
		}

		float unorm;
		float vnorm;
		float y0;
		float a;
		float umin;
		float vmin;

		CMtx3x3f M;	
		float yToL[1088];
	};

	inline CVec3f BGRToLUV(const CVec3f& bgr, const LUVConversionConstants& lc)
	{
		CVec3f xyz = lc.M * bgr;
		xyz.z = 1.f / (xyz.x + 15*xyz.y + 3*xyz.z + 1e-12f);

		const float l = lc.yToL[static_cast<int>(xyz.y*1024)];
		
		return CVec3f(l, 
			l * (52*xyz.x*xyz.z - 13*lc.unorm) - lc.umin,
			l * (117*xyz.y*xyz.z - 13*lc.vnorm) - lc.vmin);
	}

	template<typename SrcElT, typename DstElT>
	struct RGBToLUVPiotrOp
	{
	public:
		typedef float ImplSrcElType;
		typedef float ImplDstElType;
		typedef LUVConversionConstants ParamType;
		typedef int TmpType;

		static const int NumSrcBands = 3;
		static const int NumDstBands = 3;
		static const int NumSrcElPerGroup = 3;
		static const int NumDstElPerGroup = 3;
		static const int NumGroupsPerOpGeneric = 1;	

	public:
		void EvalGeneric(const ImplSrcElType* pSrc, ImplDstElType* pDst, 
			ParamType* p, TmpType&) 
		{
			const CVec3f luv = BGRToLUV(CVec3f(pSrc[0], pSrc[1], pSrc[2]), *p);
			
			pDst[0] = luv.x; 
			pDst[1] = luv.y;
			pDst[2] = luv.z;
		}
	};

	/*template <typename OT, typename IT>
	HRESULT ConvertToLUVPiotr(CPlanarTypedImg<OT>& dst, 
		const CPlanarTypedImg<IT>& src)
	{
		VT_HR_BEGIN();

		VT_HR_EXIT(src.IsValid() ? S_OK : E_INVALIDSRC);
		VT_HR_EXIT(src.Bands() == 3 ? S_OK : E_INVALIDSRC);

		if (dst.Width() != src.Width() || dst.Height() != src.Height()
			|| dst.Bands() != 3)
		{
			VT_HR_EXIT(dst.Create(src.Width(), src.Height(), 3));
		}

		LUVConversionConstants lc;

		for (int y = 0; y < src.Height(); ++y)
		{
			const IT* b = src.Band(0).Ptr(y); 
			const IT* g = src.Band(1).Ptr(y); 
			const IT* r = src.Band(2).Ptr(y); 
			OT* l = dst.Band(0).Ptr(y);
			OT* u = dst.Band(1).Ptr(y);
			OT* v = dst.Band(2).Ptr(y);

			for (int x = 0; x < src.Width(); ++x)
			{
				float bf, gf, rf;
				VtConv(&bf, b[x]);
				VtConv(&gf, g[x]);
				VtConv(&rf, r[x]);

				auto luv = BGRToLUV(CVec3f(bf, gf, rf), lc);
				
				VtConv(&l[x], luv.x);
				VtConv(&u[x], luv.y);
				VtConv(&v[x], luv.z);
			}
		}

		VT_HR_END();
	}*/
}

using namespace vt_internal;

namespace vt
{
	HRESULT VtConvertImageRGBToLUVPiotr(CByteImg& dst, const CRGBByteImg& src)
	{
		VT_HR_BEGIN();

		if (!src.IsValid())
			VT_HR_EXIT(E_INVALIDSRC);

		if (src.IsSharingMemory(dst))
			VT_HR_EXIT(E_INVALIDDST);

		if (dst.Width() != src.Width() || dst.Height() != src.Height()
			|| dst.Bands() != 3)
		{
			VT_HR_EXIT(InitDst(dst, src.Width(), src.Height(),
				VT_IMG_MAKE_TYPE(EL_FORMAT(src.GetType()), 3)));
		}
		
		LUVConversionConstants param;
		VT_HR_EXIT(UnaryImgOpDD<RGBToLUVPiotrOp>(src, dst, &param));

		VT_HR_END();
	}
}