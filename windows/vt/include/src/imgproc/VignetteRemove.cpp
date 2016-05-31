#include "StdAfx.h"

#include "vt_VignetteRemove.h"

using namespace vt;

#define F_INV_MAX						1.f // max scale of 4 at image corners
#define N_SAMPLES						63
#define N_HIST_LOG						500
#define ORIG_RANGE						256
#define PENALIZE_START					(ORIG_RANGE-10)

#define IS_SATURATED(x)                 (x >= 1.f)

#if (defined(_M_IX86) || defined(_M_AMD64))
int CorrectSpanSSE_B2F(float *q, Byte *p, float c2, float r2, int w, float f_inv)
{
    int c;

    __m128i x0 = _mm_setzero_si128();
    __m128 xf1 = _mm_set1_ps(1.f);
    __m128 xf2 = _mm_set1_ps(255.f);
    __m128 xff = _mm_set1_ps(f_inv);
    
    __m128 xfr2 = _mm_set1_ps(r2);
    xfr2 = _mm_mul_ps(xfr2, xfr2);

    for (c = 0; c < w - 3; c2 += 4, c += 4)
	{
        // load 4 1-byte pixels
        __m128i xp = _mm_loadu_si128((__m128i *) p);

        // convert 8bit to 16bit
        xp = _mm_unpacklo_epi8(xp, x0);    

        // convert 16bit to 32bit
        xp = _mm_unpacklo_epi16(xp, x0);

        // convert 32bit to float
        __m128 xfp = _mm_cvtepi32_ps(xp);

        // If origninal channel is too high or low don't trust it;
        // set output to 0 so it won't be counted in the histogram.
        __m128 xfl = _mm_cmple_ps(xfp, xf1);
        __m128 xfc = _mm_cmpge_ps(xfp, xf2);
        xfl = _mm_or_ps(xfl, xfc);
        xfp = _mm_andnot_ps(xfl, xfp);
        xfl = _mm_andnot_ps(xfl, xf1);

        // find radii = sqrt(c2 * c2 + r2 * r2)
        xfc = _mm_set_ps(c2, c2 + 1.f, c2 + 2.f, c2 + 3.f);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xfr2);
        xfc = _mm_sqrt_ps(xfc);

        // find attentuation factors = (1 + (f_inv * radius) ^ 2) ^ 2
        xfc = _mm_mul_ps(xfc, xff);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xf1);
        xfc = _mm_mul_ps(xfc, xfc);

        // apply attenuation factors and add 1 for log
        xfp = _mm_mul_ps(xfp, xfc);
        xfp = _mm_add_ps(xfp, xfl);

        // store 4 1-float pixels
        _mm_storeu_ps(q, xfp);

        p += 4;
        q += 4;
	}

    return c;
}
#endif

void CorrectSpan_B2F(float *q, Byte *p, int c2, int r2, int w, float f_inv)
{
	for (int c = 0; c < w; c2++ , c++)
	{
        // If origninal channel is too high or low don't trust it;
        // set output to 0 so it won't be counted in the histogram.
        if (p[c] <= 1 || p[c] >= 255)
            q[c] = 0.f;
        else
        {
    		float radius = _hypotf((float) r2, (float) c2);
            float factor = f_inv * radius;
            factor = 1.f + factor * factor;
            q[c] = p[c] * factor * factor + 1.f; // attenuation factor
        }
	}
}

void CorrectImage_B2F(CFloatImg& dst, const CRect& rct, CByteImg &src,
                      int w, int h, float f_inv)
{
    int c = 0, r, r2;

    VT_ASSERT( dst.Bands() == 1 && src.Bands() == 1 );

    for (r2 = rct.top - h/2, r = 0; r < rct.Height(); r2++, r++)
	{
		Byte *p = src.Ptr(r);
        float *q = dst.Ptr(r);

#if (defined(_M_IX86) || defined(_M_AMD64))
        if (g_SupportSSE2())
            c = CorrectSpanSSE_B2F(q, p, (float) (rct.left - w/2), (float) r2, rct.Width(), f_inv);
#endif
        if (c < rct.Width())
            CorrectSpan_B2F(q + c, p + c, rct.left - w/2 + c, r2, rct.Width() - c, f_inv);
	}
}

void CorrectSpanY(float *f, int c2, int r2, int w,
                  float f_inv, float global_intensity_scale)
{
	for (int c = 0; c < w; c2++ , c++)
	{
		float radius = _hypotf((float) r2, (float) c2);
		float factor = f_inv * radius;
        factor = 1.f + factor * factor;
        f[c] = factor * factor * global_intensity_scale; // attenuation factor
    }
}

#if (defined(_M_IX86) || defined(_M_AMD64))
int CorrectSpanY_SSE(float *f, float c2, float r2, int w, int bands,
                     float f_inv, float global_intensity_scale)
{
    int c;

    //__m128i x0 = _mm_setzero_si128();
    __m128 xf1 = _mm_set1_ps(1.f);
    __m128 xff = _mm_set1_ps(f_inv);
    __m128 xfg = _mm_set_ps(global_intensity_scale, global_intensity_scale, global_intensity_scale, 1.f);
   
    __m128 xfr2 = _mm_set1_ps(r2);
    xfr2 = _mm_mul_ps(xfr2, xfr2);

    for (c = 0; c < w - 3; c2 += 4, c += 4)
	{
        // find radii = sqrt(c2 * c2 + r2 * r2)
        __m128 xfc = _mm_set_ps(c2, c2 + 1.f, c2 + 2.f, c2 + 3.f);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xfr2);
        xfc = _mm_sqrt_ps(xfc);

        // find attentuation factors = (1 + (f_inv * radius) ^ 2) ^ 2
        xfc = _mm_mul_ps(xfc, xff);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xf1);
        xfc = _mm_mul_ps(xfc, xfc);

        // apply global intensity scale
        xfc = _mm_mul_ps(xfc, xfg);

        // unpack attentuation factors into separate registers
        __m128 xfa = _mm_unpackhi_ps(xfc, xfc);
        __m128 xfb = _mm_unpacklo_ps(xfa, xfa);
        xfa = _mm_unpackhi_ps(xfa, xfa);
        xfc = _mm_unpacklo_ps(xfc, xfc);
        __m128 xfd = _mm_unpacklo_ps(xfc, xfc); 
        xfc = _mm_unpackhi_ps(xfc, xfc);

        // store 4 RGB(A) factors
        _mm_storeu_ps(f, xfa);
        _mm_storeu_ps(f + bands, xfb);
        _mm_storeu_ps(f + bands * 2, xfc);
        _mm_storeu_ps(f + bands * 3, xfd);

        f += bands * 4;
	}

    return c;
}

int CorrectSpanRGB_SSE(float *f, float c2, float r2, int w, int bands,
                       const float f_inv[4], float global_intensity_scale)
{
    int c;

    //__m128i x0 = _mm_setzero_si128();
    __m128 xf1 = _mm_set1_ps(1.f);
    __m128 xff = _mm_loadu_ps(f_inv);
    __m128 xfg = _mm_set_ps(global_intensity_scale, global_intensity_scale, global_intensity_scale, 1.f);
   
    __m128 xfr2 = _mm_set1_ps(r2);
    xfr2 = _mm_mul_ps(xfr2, xfr2);

    for (c = 0; c < w - 3; c2 += 4, c += 4)
	{
        // find radii = sqrt(c2 * c2 + r2 * r2)
        __m128 xfc = _mm_set_ps(c2, c2 + 1.f, c2 + 2.f, c2 + 3.f);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xfr2);
        xfc = _mm_sqrt_ps(xfc);

        // unpack radii into separate registers
        __m128 xfa = _mm_unpackhi_ps(xfc, xfc);
        __m128 xfb = _mm_unpacklo_ps(xfa, xfa);
        xfa = _mm_unpackhi_ps(xfa, xfa);
        xfc = _mm_unpacklo_ps(xfc, xfc);
        __m128 xfd = _mm_unpacklo_ps(xfc, xfc); 
        xfc = _mm_unpackhi_ps(xfc, xfc);

        // find attentuation factors = (1 + (f_inv * radius) ^ 2) ^ 2
        xfa = _mm_mul_ps(xfa, xff);
        xfa = _mm_mul_ps(xfa, xfa);
        xfa = _mm_add_ps(xfa, xf1);
        xfa = _mm_mul_ps(xfa, xfa);
        xfb = _mm_mul_ps(xfb, xff);
        xfb = _mm_mul_ps(xfb, xfb);
        xfb = _mm_add_ps(xfb, xf1);
        xfb = _mm_mul_ps(xfb, xfb);
        xfc = _mm_mul_ps(xfc, xff);
        xfc = _mm_mul_ps(xfc, xfc);
        xfc = _mm_add_ps(xfc, xf1);
        xfc = _mm_mul_ps(xfc, xfc);
        xfd = _mm_mul_ps(xfd, xff);
        xfd = _mm_mul_ps(xfd, xfd);
        xfd = _mm_add_ps(xfd, xf1);
        xfd = _mm_mul_ps(xfd, xfd);

        // apply global intensity scale
        xfa = _mm_mul_ps(xfa, xfg);
        xfb = _mm_mul_ps(xfb, xfg);
        xfc = _mm_mul_ps(xfc, xfg);
        xfd = _mm_mul_ps(xfd, xfg);

        // store 4 RGB(A) factors
        _mm_storeu_ps(f, xfa);
        _mm_storeu_ps(f + bands, xfb);
        _mm_storeu_ps(f + bands * 2, xfc);
        _mm_storeu_ps(f + bands * 3, xfd);

        f += bands * 4;
	}

    return c;
}
#endif

void CorrectSpan(float *f, int c2, int r2, int w, int bands,
                 const float f_inv[4], float global_intensity_scale)
{
	for (int c = 0; c < w; c2++ , c++)
	{
		float radius = _hypotf((float) r2, (float) c2);
		float factor[3];
        for (int i = 0; i < 3; i++)
        {
            factor[i] = f_inv[i] * radius;
            factor[i] = 1.f + factor[i] * factor[i];
            f[i] = factor[i] * factor[i] * global_intensity_scale; // attenuation factor
        }

        // preserve alpha
        if (bands > 3)
            f[3] = 1.f;

        f += bands;
    }
}

void CorrectImage(CImg& dst, const CRect& rct, const CImg &src,
                  int w, int h, const float f_inv[4], float global_intensity_scale)
{
    int r, r2;

    VT_ASSERT( dst.Bands() == src.Bands() );

    // determine source and dest types
    int iSrcType  = EL_FORMAT(src.GetType());
    int iDstType  = EL_FORMAT(dst.GetType());
    int iSrcBands = src.Bands();

    // pixel buffer
    const int bufsize = 256;
    float tmpline[256 * 4];
    float facline[256 * 4];

    int iSpan = rct.Width();
    for (r2 = rct.top - h/2, r = 0; r < rct.Height(); r2++, r++)
	{
		const Byte *p = src.BytePtr(r);
        Byte *q = dst.BytePtr(r);

        int c = 0;
        for (int i = 0; i < iSpan; i += bufsize)
        {
            int iCurSpan = VtMin(bufsize, iSpan - i);
             
            // setup temp buffer as necessary
            float *pf, *qf;
            if (iSrcType != EL_FORMAT_FLOAT)
            {
                pf = tmpline;
                VtConvertSpan(pf, VT_IMG_MAKE_TYPE(EL_FORMAT_FLOAT, iSrcBands),
					          p, VT_IMG_MAKE_TYPE(iSrcType, iSrcBands), 
							  iCurSpan * iSrcBands);
            }
            else
            {
                pf = (float *) p;
            }
            qf = (iDstType != EL_FORMAT_FLOAT) ? tmpline : (float *) q;

            if (dst.Bands() == 1)
                CorrectSpanY(facline, rct.left - w/2 + c, r2, iCurSpan,
                             f_inv[0], global_intensity_scale);
            else
            {
                int n = 0;
#if (defined(_M_IX86) || defined(_M_AMD64))
                if (g_SupportSSE2())
                {
                    if (f_inv[0] == f_inv[1] && f_inv[1] == f_inv[2])
                        n = CorrectSpanY_SSE(facline, (float) (rct.left - w/2 + c), (float) r2,
                                             iCurSpan, dst.Bands(), f_inv[0], global_intensity_scale);
                    else
                        n = CorrectSpanRGB_SSE(facline, (float) (rct.left - w/2 + c), (float) r2,
                                               iCurSpan, dst.Bands(), f_inv, global_intensity_scale);
                }
#endif
                if (n < iCurSpan)
                {
                    CorrectSpan(facline + n * iSrcBands, rct.left - w/2 + c + n, r2,
                                iCurSpan - n, dst.Bands(), f_inv, global_intensity_scale);
                }
            }

            // apply attenuation factor and global intensity scale
            VtMulSpan(qf, facline, pf, iCurSpan * iSrcBands);
    
            // if necessary convert to dest type  
            if (iDstType != EL_FORMAT_FLOAT)
            {
                VtConvertSpan(q, VT_IMG_MAKE_TYPE(iDstType, iSrcBands), 
					          qf, VT_IMG_MAKE_TYPE(EL_FORMAT_FLOAT, iSrcBands),
							  iCurSpan * iSrcBands);
            }

            q += iCurSpan * src.PixSize();
            p += iCurSpan * dst.PixSize();
            c += iCurSpan;
        }
    }
}

float EstimateGlobalIntensityScale(CFloatImg& img, float* f_inv)
{
    int r, c, n_sat0[3] = { 0, 0, 0 };
    int width = img.Width(), height = img.Height();
    int bands = img.Bands() >= 3 ? 3 : 1;
	int max_delta = F2I( (float) width * (float) height * 0.05f );
    float excess_sat_total[3] = { 0, 0, 0 };
    float excess_orig_total[3] = { 0, 0, 0 };
    int n_excess_sat[3] = { 0, 0, 0 }; // number of excess saturated pixels (in corrected image)
	CFloatImg correctedImg(width, height, img.Bands());

	float global_intensity_scale = 1.0f;
	CorrectImage(correctedImg, img.Rect(), img, width, height, f_inv, 1.f);

	// now find the global intensity scale <= 1.0 that minimizes the change
	// in the number of saturated pixels (maximum of 5% of all pixels)
	// find the number of saturated pixels in the corrected image that
	// isn't already saturated in the original (possibly downsampled)
	// image
	for (r = 0; r < height; r++)
	{
		float *p = img.Ptr(r);
		float *q = correctedImg.Ptr(r);

		for (c = 0; c < width; c++)
		{
            for (int i = 0; i < bands; i++)
            {
			    if (IS_SATURATED(p[i])) 
				    n_sat0[i]++;
			    else if (IS_SATURATED(q[i]))
			    {
				    excess_orig_total[i] += p[i];
				    excess_sat_total[i] += q[i];
				    n_excess_sat[i]++;
			    }
            }

            p += img.Bands();
            q += correctedImg.Bands();
		}
	}

    for (int i = 0; i < bands; i++)
    {
	    if (n_excess_sat[i] > max_delta)
	    {
		    global_intensity_scale = VtMin(global_intensity_scale,
                (float) n_excess_sat[i] * 1.f / excess_sat_total[i]);
	    }

        VT_DEBUG_LOG( "IP_LOG_VIGNETTE: EstimateGlobalIntensityScale: No. of saturated pixels in original channel = %d\n", n_sat0[i]);
	    VT_DEBUG_LOG( "IP_LOG_VIGNETTE: EstimateGlobalIntensityScale: No. of excess saturated pixels in corrected channel = %d\n", n_excess_sat[i]);
    }

    VT_DEBUG_LOG( "IP_LOG_VIGNETTE: EstimateGlobalIntensityScale: global_intensity_scale = %f\n", global_intensity_scale);

    return global_intensity_scale;
}

inline bool IsTooDark(int iHist[256], int limit)
{
    int count = 0;
    for (int i = 0; i <= 25; i++)
        count += iHist[i];
    return count > limit;
}

inline bool IsTooBright(int iHist[256], int limit)
{
    int count = 0;
    for (int i = 250; i <= 255; i++)
        count += iHist[i];
    return count > limit;
}

HRESULT CVignetteAnalyze::AnalyzeImage(IImageReader* pReader, UINT uSourceId)
{
    VT_HR_BEGIN()
    
    if (pReader == NULL)
        VT_HR_EXIT( E_INVALIDARG );

    // Make sure all the sizes are the same.
    CImgInfo info = pReader->GetImgInfo();
    if (m_info.type != OBJ_UNDEFINED &&
        !((m_info.width == info.width && m_info.height == info.height) ||
          (m_info.width == info.height && m_info.height == info.width)))
        VT_HR_EXIT( E_INVALIDARG );
    if (info.Bands() != 1 && info.Bands() != 3 && info.Bands() != 4)
        VT_HR_EXIT( E_INVALIDARG );
    m_info = info;

    // Read RGB and Luminance images.
    int width = info.width, height = info.height, bands = info.Bands();
    CByteImg img, imgY;
    VT_HR_EXIT( pReader->ReadImg(img) );
    if (bands == 1)
        img.Share(imgY);
    else
    {
        VT_HR_EXIT( imgY.Create(width, height) );
        VT_HR_EXIT( VtConvertImage(imgY, img) );
    }

    // Build histograms for RGB and Luminance images.
    int iHist[4][256];
    VT_HR_EXIT( VtBuildHistogram(imgY, iHist[0]) );
    if (bands >= 3)
       VT_HR_EXIT( VtBuildHistogram(img, iHist[1], iHist[2], iHist[3]) );

    // Select candidate image with minimum entropy.
    int channel = 0;
    for (; channel < (bands == 1 ? 1 : 4); channel++)
    {
        int limit = (width * height) >> 4;  // dark/bright limit is 1/16 of pixels
        if (IsTooDark(iHist[channel], limit))
        {
        	VT_DEBUG_LOG( "IP_LOG_VIGNETTE: CVignetteAnalyze::AnalyzeImage: image %d channel %d too dark\n",
                uSourceId, channel );
            break;
        }
        else if (IsTooBright(iHist[channel], limit))
        {
        	VT_DEBUG_LOG( "IP_LOG_VIGNETTE: CVignetteAnalyze::AnalyzeImage: image %d channel %d too bright\n",
                uSourceId, channel );
            break;
        }
    }

    // If any of R,G,B too bright/dark force Luminance.
    if (channel > 0 && channel < 4)
        channel = 1;

    // Calculate discrete entropy.
    float ent[4] = { 0, 0, 0, 0 };
    for (int c = 0; c < channel; c++)
    {
        for (int h = 0; h < 256; h++)
            ent[c] -= iHist[c][h] <= 1 ? 0 : iHist[c][h] * log((float) iHist[c][h]);
    }

    // Pick Luminance or worst (highest entropy) of R,G,B for comparison.
    float fEntropyY = ent[0];
    float fEntropyRGB = VtMax(ent[1], VtMax(ent[2], ent[3]));

	VT_DEBUG_LOG( "IP_LOG_VIGNETTE: CVignetteAnalyze::AnalyzeImage: entropy = %f, %f\n",
        fEntropyY, fEntropyRGB );

    // Use RGB entropy if possible, else Luminance.
    if (m_fEntropyRGB > fEntropyRGB)
	{
		m_fEntropyRGB = fEntropyRGB;
        m_fEntropyY = 0;
		m_uSourceId = uSourceId;
    }
    if (m_fEntropyRGB == 0 && m_fEntropyY > fEntropyY)
    {
		m_fEntropyY = fEntropyY;
		m_uSourceId = uSourceId;
    }

	VT_HR_EXIT_LABEL()

    if (hr != S_OK)
        Reset();

    return hr;
}

float SSE(CFloatImg& correctedBW, CFloatLogHistogram& hist)
//float SSE(CFloatImg& correctedBW, CByteImg& imgBW)
{
	float sse = 0;
	float hist_log[N_HIST_LOG];

    hist.Clear();
    //ULONG uTotCnt = BuildLogHistogram(correctedBW, hist); // uTotCnt is not used
    BuildLogHistogram(correctedBW, hist);
#if 0
	int r, c;

    // I -> log (1 + I)/log 256

    float factor = ORIG_RANGE / log(256.f);

	for (int i = 0; i < N_HIST_LOG; i++)
		hist_log[i] = 0;
	for (r = 0; r < imgBW.Height(); r++)
	{
		Byte *p0 = imgBW.Ptr(r);
		float *p = correctedBW.Ptr(r);

		for (c = 0; c < imgBW.Width(); c++)
		{
			if (p0[c] <= 1 || p0[c] >= 255)
                continue;

            if (p[c] <= 0) 
			{
				continue;//hist_log[0] += 1;
			}
			else
			{
				int index = F2I( factor * log(p[c]) );

				if (index >= N_HIST_LOG)
                    index = N_HIST_LOG - 1;
				hist_log[index] += 1;
				if (index>0) hist_log[index-1] += 0.5;
				if (index>1) hist_log[index-2] += 0.25;
				//if (index>2) hist_log[index-3] += 0.125;
				if (index<N_HIST_LOG-1) hist_log[index+1] += 0.5;
				if (index<N_HIST_LOG-2) hist_log[index+2] += 0.25;
				//if (index<N_HIST_LOG-3) hist_log[index+3] += 0.125;
			}
		}
	}
#endif

    SmoothHistogram(N_HIST_LOG, hist.GetBuffer(), hist_log);

    for (UINT i = 0; i < PENALIZE_START; i++)
	{
        float n = hist_log[i];
        if (n <= 1)
            continue;
        sse -= n * log(n);
	}

    // de-emphasize rescaled intensities exceeding the dynamic range
	for (UINT i = PENALIZE_START; i < N_HIST_LOG; i++)
	{
        float n = hist_log[i];
        if (n <= 1)
            continue;
		sse -= (1.f - (float) (i - PENALIZE_START) /
            (float) (N_HIST_LOG - PENALIZE_START - 1)) *
            n * log(n);
	}

	return sse;
}

HRESULT vt::ComputeVignetteParams(IImageReader* pReader, eVignetteRemove eVR,
                                  VIGNETTE_PARAMS& params,
                                  CTaskProgress* pProgress)
{
    VT_HR_BEGIN()

    if (pReader == NULL)
        return E_INVALIDARG;

    if (eVR == eVRNone)
    {
        params = VIGNETTE_PARAMS();
        return S_OK;
    }

    CImgInfo info = pReader->GetImgInfo();

    // Find f_inv which minimizes SSE for candidate image.
    int width = info.width, height = info.height, bands = info.Bands();
    CByteImg img(width, height, eVR == eVRLum ? 1 : bands);
    VT_HR_EXIT( pReader->ReadImg(img) );

    CFloatImg correctedBW;
    VT_HR_EXIT( correctedBW.Create(width, height) );
    CByteImg imgBW;
    if (img.Bands() == 1)
        img.Share(imgBW);
    else
        VT_HR_EXIT( imgBW.Create(width, height) );

    // Check horizontal and vertical quarter-images also.
    CRect rctH4 = CRect(0, 3 * imgBW.Height() / 8, imgBW.Width(), 5 * imgBW.Height() / 8);
    CRect rctV4 = CRect(3 * imgBW.Width() / 8, 0, 5 * imgBW.Width() / 8, imgBW.Height());
    CFloatImg imgH4, imgV4;
    VT_HR_EXIT( correctedBW.Share(imgH4, rctH4) );
    VT_HR_EXIT( correctedBW.Share(imgV4, rctV4) );
	float sse_minH[4] = { 0, 0, 0, 0 };
	float sse_minV[4] = { 0, 0, 0, 0 };

	float sse_min[4] = { 0, 0, 0, 0 };
    float f_inv[4] = { 0, 0, 0, 0 };
    float f_inv_delta = F_INV_MAX / (N_SAMPLES-1);
    int n_channels = img.Bands() == 1 ? 1 : (eVR == eVRRGB ? 3 : 4); // 1 for Lum, 3 for RGB, 4 for Auto

    // Create a log histogram with 512 entries, hist_log(256) = 256.
    CFloatLogHistogram hist_log;
    VT_HR_EXIT( hist_log.Create(9, 0, 15) );

    float fPctPerIter = 100.f / (n_channels * N_SAMPLES);
    
	// scale is the scale of src and dst compared with the original image
	float scale = 0.5f * (float) width; // assumes aspect ratios are the same

	for (int channel = 0; channel < n_channels; channel++)
	{
        if (img.Bands() > 1)
        {
            if (channel < 3)
            {
                // note: modify channel because the order is actually BGRA
	            BandIndexType band = (BandIndexType) (2 - channel);
                VT_HR_EXIT( VtConvertBands(imgBW, img, 1, &band) );
            }
            else
            {
                // last "channel" is luminance for comparison
                VT_HR_EXIT( VtConvertImage(imgBW, img) );
            }
        }

		for (int i = 0; i < N_SAMPLES; i++)
		{
			float f_try = i * f_inv_delta / scale;

            CorrectImage_B2F(correctedBW, imgBW.Rect(), imgBW, width, height, f_try);

            // Horizontal and vertical quarter-images must also have improved entropy.
		    float sseH = SSE(imgH4, hist_log);
		    float sseV = SSE(imgV4, hist_log);

			float sse = SSE(correctedBW, hist_log);
            //float sse = SSE(correctedBW, imgBW);
			VT_DEBUG_LOG( "ComputeVignetteFactor: channel %d f_inv->sse %f->%f %f %f\n",
                channel, f_try, sse, sseH, sseV );

			if ((sse_min[channel] > sse) && 
				(sse_minH[channel]-sseH) > -0.005f*fabs(sse_minH[channel]) && 
				(sse_minV[channel]-sseV) > -0.005f*fabs(sse_minV[channel]) )
			{
				f_inv[channel] = f_try;
			}
            sse_min[channel] = VtMin(sse_min[channel], sse);
            sse_minH[channel] = VtMin(sse_minH[channel], sseH);
            sse_minV[channel] = VtMin(sse_minV[channel], sseV);

		    // Report progress.
			VT_HR_EXIT( CheckTaskProgressCancel(
				pProgress, float(channel * N_SAMPLES + i + 1) * fPctPerIter) );
		}
	}

	SetProgress100Percent(pProgress);

    float factor[4] = { 0, 0, 0, 0 };
	for (int channel = 0; channel < n_channels; channel++)
	{
        factor[channel] = 1.f + f_inv[channel] * f_inv[channel] * scale * scale;
        factor[channel] *= factor[channel];

        VT_DEBUG_LOG( "IP_LOG_VIGNETTE: ComputeVignetteFactor: channel %d: SSE min = %f, f_inv = %f, factor = %f\n",
			channel, sse_min[channel], f_inv[channel], factor[channel] );
	}

    if (n_channels == 1)
        f_inv[2] = f_inv[1] = f_inv[0];
    else if (n_channels == 4)
    {
        // Check if per channel factors are close enough to luminance.
        if (abs(factor[0] - factor[3]) / factor[3] > 0.1f ||
            abs(factor[1] - factor[3]) / factor[3] > 0.1f ||
            abs(factor[2] - factor[3]) / factor[3] > 0.1f)
        {
            VT_DEBUG_LOG( "IP_LOG_VIGNETTE :ComputeVignetteFactor: using luminance f_inv %f\n", f_inv[3] );
            f_inv[0] = f_inv[1] = f_inv[2] = f_inv[3];
        }
    }

    params.fFactorR = f_inv[0] * scale;
    params.fFactorG = f_inv[1] * scale;
    params.fFactorB = f_inv[2] * scale;

	// estimate a reasonable global intensity scale 
    CFloatImg imgFlt;
    VT_HR_EXIT( imgFlt.Create(width, height, img.Bands()) );
    VT_HR_EXIT( pReader->ReadImg(imgFlt) );
    params.fGlobalIntensityScale = EstimateGlobalIntensityScale(imgFlt, f_inv);
    
    VT_HR_END()
}

//+-----------------------------------------------------------------------------
//
// Class: CVignetteRemoveTransform
// 
// Synposis: IImageTransform implementation of CVignetteRemove::CorrectImage()
// 
//------------------------------------------------------------------------------
HRESULT
CVignetteRemoveTransform::Transform(OUT CImg* pimgDst, 
							        IN  const CRect& rctLayerDst,
							        IN  const CImg& imgSrc)
{
    VT_HR_BEGIN()

    VT_PTR_EXIT( pimgDst );

    CorrectImage(*pimgDst, rctLayerDst, imgSrc,
                 m_w, m_h, m_f_inv, m_params.fGlobalIntensityScale);

    VT_HR_END()
}
