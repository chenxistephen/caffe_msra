//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Separable filter specialization for 14641 kernel
//
//------------------------------------------------------------------------

#include "stdafx.h"

#include "../common/vt_intrinsics.h"
#include "vt_image.h"
#include "vt_utils.h"
#include "vt_separablefilter.h"

using namespace vt;

//-----------------------------------------------------------------------------
// routines to apply 14641 upsample filter to byte data; uses 16 bit accumulation
// so full precision of intermediate results are retained
//-----------------------------------------------------------------------------

// produce 8.3 16bpp result for vertical (first) pass
inline uint16_t FilterX4X4XP1(uint8_t a, uint8_t b)
{
    uint16_t a16 = (uint16_t)a;
    uint16_t b16 = (uint16_t)b;
    return (a16+b16)<<2;
}
inline uint16_t Filter1X6X1P1(uint8_t l, uint8_t c, uint8_t r)
{
    uint16_t l16 = (uint16_t)l;
    uint16_t c16 = (uint16_t)c;
    uint16_t r16 = (uint16_t)r;
    return ( (c16<<2) + (c16<<1) + l16 + r16 );
}

// produce rounded 8.0 result from 8.3 inputs for horizontal (second) pass
inline uint8_t FilterX4X4XP2(uint16_t a, uint16_t b)
{
    uint16_t acc = a+b;
    return (Byte)((acc+(1<<3))>>4);
}
inline uint8_t Filter1X6X1P2(uint16_t l, uint16_t c, uint16_t r)
{
    uint16_t acc = (c<<2)+(c<<1)+(l)+(r);
    return (Byte)((acc+(1<<5))>>6);
}

//-----------------------------------------------------------------------------
// routines to apply 14641 upsample filter to float data;
//-----------------------------------------------------------------------------

// produce scaled by 8.f results for vertical (first) pass
inline float FilterX4X4XP1(float a, float b)
{
    return (a+b)*4.f;
}
inline float Filter1X6X1P1(float l, float c, float r)
{
    return (l+(6.f*c)+r);
}

// produce result from 8.f scaled inputs for horizontal (second) pass
inline float FilterX4X4XP2(float a, float b)
{
    return (a+b)*.0625f;
}
inline float Filter1X6X1P2(float l, float c, float r)
{
    return (l+(6.f*c)+r)*(.015625f);
}

//-----------------------------------------------------------------------------
// apply 14641 upsample to 4 channels
//-----------------------------------------------------------------------------
template <typename Tdst, typename Tsrc>
void Filter4byX4X4XP1(Tdst* dst, const Tsrc*srca, const Tsrc* srcb)
{
    *(dst+0) = FilterX4X4XP1(*(srca+0), *(srcb+0));
    *(dst+1) = FilterX4X4XP1(*(srca+1), *(srcb+1));
    *(dst+2) = FilterX4X4XP1(*(srca+2), *(srcb+2));
    *(dst+3) = FilterX4X4XP1(*(srca+3), *(srcb+3));
}

template <typename Tdst, typename Tsrc>
void Filter4by1X6X1P1(Tdst* dst,const  Tsrc*srca, const Tsrc* srcb, const Tsrc* srcc)
{
    *(dst+0) = Filter1X6X1P1(*(srca+0), *(srcb+0), *(srcc+0));
    *(dst+1) = Filter1X6X1P1(*(srca+1), *(srcb+1), *(srcc+1));
    *(dst+2) = Filter1X6X1P1(*(srca+2), *(srcb+2), *(srcc+2));
    *(dst+3) = Filter1X6X1P1(*(srca+3), *(srcb+3), *(srcc+3));
}

template <typename Tdst, typename Tsrc>
void Filter4byX4X4XP2(Tdst* dst, const Tsrc*srca, const Tsrc* srcb)
{
    *(dst+0) = FilterX4X4XP2(*(srca+0), *(srcb+0));
    *(dst+1) = FilterX4X4XP2(*(srca+1), *(srcb+1));
    *(dst+2) = FilterX4X4XP2(*(srca+2), *(srcb+2));
    *(dst+3) = FilterX4X4XP2(*(srca+3), *(srcb+3));
}

template <typename Tdst, typename Tsrc>
void Filter4by1X6X1P2(Tdst* dst,const  Tsrc*srca, const Tsrc* srcb, const Tsrc* srcc)
{
    *(dst+0) = Filter1X6X1P2(*(srca+0), *(srcb+0), *(srcc+0));
    *(dst+1) = Filter1X6X1P2(*(srca+1), *(srcb+1), *(srcc+1));
    *(dst+2) = Filter1X6X1P2(*(srca+2), *(srcb+2), *(srcc+2));
    *(dst+3) = Filter1X6X1P2(*(srca+3), *(srcb+3), *(srcc+3));
}

//-----------------------------------------------------------------------------
// SIMD routines for SSE and Neon
//-----------------------------------------------------------------------------
#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))

// vector routines do multiples of 4 destination RGBA pixels at a time
static const int Upsample2to1by14641SIMDBatch = 4;

#if (defined(_M_IX86) || defined(_M_AMD64))
#define VT_UPSAMPLE14641_SIMD_SUPPORTED (g_SupportSSSE3())
#elif defined(_M_ARM)
#define VT_UPSAMPLE14641_SIMD_SUPPORTED true
#endif

// 'between' source pixels vertically; float I/O
template<void (*pfnStore)(__m128i*, const __m128i&)>
int Upsample2to1by14641ARowSIMD(float* pd, 
                               const float* pst, const float* psb,
                               int wd, const float* lvfpix, const float* rvfpix)
{
    int ret = wd&(~(Upsample2to1by14641SIMDBatch-1)); // number of dest pixels completed
    const float* pstlv = pst + (ret/2)*4;

    __m128 x0,x1;
    __m128 L,CL,CR,R;
    __m128 const4 = _mm_set1_ps(4.f);
    __m128 const6 = _mm_set1_ps(6.f);
    __m128 constNormX4X4X = _mm_set1_ps(.0625f);
    __m128 constNorm1X6X1 = _mm_set1_ps(.015625f);

    // TODO: aligned loads on sources

    // load L(eft) pix from top of left pix buffer
    L = _mm_load_ps(lvfpix+4);

    // load and vertically filter one RBBA pixels for C(enter)L(eft)
    x0 = _mm_loadu_ps(pst); pst += 4;
    x1 = _mm_loadu_ps(psb); psb += 4;
    x0 = _mm_add_ps(x0,x1); 
    CL = _mm_mul_ps(x0,const4); // vertically filtered and scaled by 8.f

    while ( wd >= Upsample2to1by14641SIMDBatch )
    {
        // load and vertically filter one RBBA pixels for C(enter)R(ight)
        x0 = _mm_loadu_ps(pst); pst += 4;
        x1 = _mm_loadu_ps(psb); psb += 4;
        x0 = _mm_add_ps(x0,x1); 
        CR = _mm_mul_ps(x0,const4); // vertically filtered and scaled by 8.f

        // load and vertically filter one pixel for R(ight)
        if (pst == pstlv)
        {
            R = _mm_load_ps(rvfpix);
        }
        else
        {
            x0 = _mm_loadu_ps(pst); pst += 4;
            x1 = _mm_loadu_ps(psb); psb += 4;
            x0 = _mm_add_ps(x0,x1);
            R = _mm_mul_ps(x0,const4); // vertically filtered and scaled by 8.f
        }

        // L,CL,CR,R
        // R0 = 1X6X1(L,CL,CR)
        // R1 = X4X4X(CL,CR)
        // R2 = 1X6X1(CL,CR,R)
        // R3 = X4X4(CR,R)

        // first result pixel is 1X6X1 filtered
        x0 = _mm_add_ps(L,CR);
        x1 = _mm_mul_ps(CL,const6);
        x0 = _mm_add_ps(x0,x1); // horizontally filtered by 1x6x1
        x0 = _mm_mul_ps(x0,constNorm1X6X1); // normalized
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // second result pixel is X4X4X filtered
        x0 = _mm_add_ps(CL,CR);
        x0 = _mm_mul_ps(x0,constNormX4X4X);
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // third result pixel is 1X6X1 filtered
        x0 = _mm_add_ps(CL,R);
        x1 = _mm_mul_ps(CR,const6);
        x0 = _mm_add_ps(x0,x1); // horizontally filtered by 1x6x1
        x0 = _mm_mul_ps(x0,constNorm1X6X1); // normalized
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // second result pixel is X4X4X filtered
        x0 = _mm_add_ps(CR,R);
        x0 = _mm_mul_ps(x0,constNormX4X4X);
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        L = CR;
        CL = R;
        wd -= Upsample2to1by14641SIMDBatch;
    }
    return ret;
}

// 'between' source pixels vertically; Byte I/O
template<void (*pfnStore)(__m128i*, const __m128i&)>
int Upsample2to1by14641ARowSIMD(Byte* pd, 
                               const Byte* pst, const Byte* psb,
                               int wd, const uint16_t* lvfpix, const uint16_t* rvfpix)
{
    int ret = wd&(~(Upsample2to1by14641SIMDBatch-1)); // number of dest pixels completed
    const Byte* pstlv = pst + (ret/2)*4;

#if (defined(_M_IX86) || defined(_M_AMD64))
    __m128i zero = _mm_setzero_si128();
#endif
    __m128i x0,x1,x2,x3;
    __m128i L,C,R;

    // load L(eft) pix from left pix buffer
    L = _mm_load_si128((const __m128i*)lvfpix);

    // load and vertically-filter two RGBA pixels for C(enter)
#if (defined(_M_IX86) || defined(_M_AMD64))
    x0 = _mm_loadl_epi64((const __m128i*)pst); pst += (4*2);
    x1 = _mm_loadl_epi64((const __m128i*)psb); psb += (4*2);
    x0 = _mm_unpacklo_epi8(x0,zero);
    x1 = _mm_unpacklo_epi8(x1,zero);
    x0 = _mm_adds_epu16(x0,x1);
#elif (defined(_M_ARM))
    {
        __n64 x0n = vld1_u8(pst); pst += (4*2);
        __n64 x1n = vld1_u8(psb); psb += (4*2);
        x0 = vaddl_u8(x0n,x1n);
    }
#endif
    C = _mm_slli_epi16(x0,2); // 8.3 vertically filtered

    while ( wd >= Upsample2to1by14641SIMDBatch )
    {
        if (pst == pstlv)
        {
            R = _mm_load_si128((const __m128i*)rvfpix);
        }
        else
        {
            // load and vertically-filter two RGBA pixels for R(ight)
#if (defined(_M_IX86) || defined(_M_AMD64))
            x0 = _mm_loadl_epi64((const __m128i*)pst); pst += (4*2); 
            x1 = _mm_loadl_epi64((const __m128i*)psb); psb += (4*2); 
            x0 = _mm_unpacklo_epi8(x0,zero);
            x1 = _mm_unpacklo_epi8(x1,zero);
            x0 = _mm_adds_epu16(x0,x1);
#elif (defined(_M_ARM))
            {
                __n64 x0n = vld1_u8(pst); pst += (4*2);
                __n64 x1n = vld1_u8(psb); psb += (4*2);
                x0 = vaddl_u8(x0n,x1n);
            }
#endif
            R = _mm_slli_epi16(x0,2); // 8.3 vertically filtered
        }

        // horizontal filtering
        //
        // for 4 results R0,R1,R2,R3:
        // R0,R2 are 1x6x1 filtered, and R1,R3 are x4x4x filtered
        // L,C,R each have two pixels L0,L1, etc.
        // R0 = 1x6x1(L1,C0,C1)
        // R1 = x4x4x(C0,C1)
        // R2 = 1x6x1(C0,C1,R0)
        // R3 = x4x4x(C1,R0)
        //
        // so can combine computation for R0,R2 and R1,R3:
        //  R0R2 = 1x6x1(L1C0,C0C1,C1R0)
        //  R1R3 = x4x4x(C0C1,C1R0)
        // where:
        //  L1C0 = (L<<64)|(C>>64)
        //  C0C1 = C
        //  C1R0 = (C<<64)|(R>>64)
        // then rejigger R0R2,R1R3 into R0R1R2R3 for store
        //
#if (defined(_M_IX86) || defined(_M_AMD64))
        x3 = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(L),_mm_castsi128_pd(C),1)); // L1C0
        x2 = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(C),_mm_castsi128_pd(R),1)); // C1R0
#elif (defined(_M_ARM))
        x3 = vextq_u16(L,C,4);
        x2 = vextq_u16(C,R,4);
#endif
        x1 = _mm_adds_epu16(x2,x3); // L1C0+C1R0
        x0 = _mm_slli_epi16(C,2);
        x1 = _mm_adds_epu16(x1,x0); // += (C0C1)<<2
        x0 = _mm_slli_epi16(C,1);
        x0 = _mm_adds_epu16(x1,x0); // += (C0C1)<<1 is unnormalized R0R2 in 8.6
        x1 = _mm_adds_epu16(C,x2);  // is unnormalized R1R3 in 8.4
#if (defined(_M_IX86) || defined(_M_AMD64))
        x0 = _mm_adds_epu16(x0,_mm_set1_epi16(1<<5));
        x0 = _mm_srli_epi16(x0,6); // normalized R0R2
        x1 = _mm_adds_epu16(x1,_mm_set1_epi16(1<<3));
        x1 = _mm_srli_epi16(x1,4); // normalized R1R3
        x0 = _mm_packus_epi16(x0,x1); // Byte R0R2R1R3
        x0 = _mm_shuffle_epi32(x0,_MM_SHUFFLE(3,1,2,0)); // R0R1R2R3
#elif (defined(_M_ARM))
        {
            uint8x8x2_t x0t;
            x0t.val[0] = vrshrn_n_u16(x0,6); // normalized and packed R0R2
            x0t.val[1] = vrshrn_n_u16(x1,4); // normalized and packed R1R3
            x0t = vzip_u32(x0t.val[0],x0t.val[1]); // R0R1, R2R3
            x0 = vcombine_u64(x0t.val[0],x0t.val[1]); // R0R1R2R3
        }
#endif
        pfnStore((__m128i*)pd,x0); pd += (4*4);

        L = C;
        C = R;
        wd -= Upsample2to1by14641SIMDBatch;
    }

    return ret; 
}

// 'on' source pixels vertically; float I/O
template<void (*pfnStore)(__m128i*, const __m128i&)>
int Upsample2to1by14641BRowSIMD(float* pd, 
                               const float* pst, const float* psc, const float* psb,
                               int wd, const float* lvfpix, const float* rvfpix)
{
    int ret = wd&(~(Upsample2to1by14641SIMDBatch-1)); // number of dest pixels completed
    const float* pstlv = pst + (ret/2)*4;

    __m128 x0,x1,x2;
    __m128 L,CL,CR,R;
    __m128 const6 = _mm_set1_ps(6.f);
    __m128 constNormX4X4X = _mm_set1_ps(.0625f);
    __m128 constNorm1X6X1 = _mm_set1_ps(.015625f);

    // TODO: aligned loads on sources

    // load L(eft) pix from top of left pix buffer
    L = _mm_load_ps(lvfpix+4);

    // load and vertically filter one RBBA pixels for C(enter)L(eft)
    x0 = _mm_loadu_ps(pst); pst += 4;
    x1 = _mm_loadu_ps(psb); psb += 4;
    x2 = _mm_loadu_ps(psc); psc += 4;
    x0 = _mm_add_ps(x0,x1);
    x2 = _mm_mul_ps(x2,const6);
    CL = _mm_add_ps(x0,x2); // horizontally filtered by 1x6x1

    while ( wd >= Upsample2to1by14641SIMDBatch )
    {
        // load and vertically filter one RBBA pixels for C(enter)R(ight)
        x0 = _mm_loadu_ps(pst); pst += 4;
        x1 = _mm_loadu_ps(psb); psb += 4;
        x2 = _mm_loadu_ps(psc); psc += 4;
        x0 = _mm_add_ps(x0,x1);
        x2 = _mm_mul_ps(x2,const6);
        CR = _mm_add_ps(x0,x2); // horizontally filtered by 1x6x1

        // load and vertically filter one pixel for R(ight)
        if (pst == pstlv)
        {
            R = _mm_load_ps(rvfpix);
        }
        else
        {
            x0 = _mm_loadu_ps(pst); pst += 4;
            x1 = _mm_loadu_ps(psb); psb += 4;
            x2 = _mm_loadu_ps(psc); psc += 4;
            x0 = _mm_add_ps(x0,x1);
            x2 = _mm_mul_ps(x2,const6);
            R = _mm_add_ps(x0,x2); // horizontally filtered by 1x6x1
        }

        // L,CL,CR,R
        // R0 = 1X6X1(L,CL,CR)
        // R1 = X4X4X(CL,CR)
        // R2 = 1X6X1(CL,CR,R)
        // R3 = X4X4(CR,R)

        // first result pixel is 1X6X1 filtered
        x0 = _mm_add_ps(L,CR);
        x1 = _mm_mul_ps(CL,const6);
        x0 = _mm_add_ps(x0,x1); // horizontally filtered by 1x6x1
        x0 = _mm_mul_ps(x0,constNorm1X6X1); // normalized
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // second result pixel is X4X4X filtered
        x0 = _mm_add_ps(CL,CR);
        x0 = _mm_mul_ps(x0,constNormX4X4X);
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // third result pixel is 1X6X1 filtered
        x0 = _mm_add_ps(CL,R);
        x1 = _mm_mul_ps(CR,const6);
        x0 = _mm_add_ps(x0,x1); // horizontally filtered by 1x6x1
        x0 = _mm_mul_ps(x0,constNorm1X6X1); // normalized
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        // second result pixel is X4X4X filtered
        x0 = _mm_add_ps(CR,R);
        x0 = _mm_mul_ps(x0,constNormX4X4X);
        pfnStore((__m128i*)pd,_mm_castps_si128(x0)); pd += 4;

        L = CR;
        CL = R;
        wd -= Upsample2to1by14641SIMDBatch;
    }
    return ret;
}

// 'on' source pixel vertically; Byte I/O
template<void (*pfnStore)(__m128i*, const __m128i&)>
int Upsample2to1by14641BRowSIMD(Byte* pd, 
                               const Byte* pst, const Byte* psc, const Byte* psb,
                               int wd, const uint16_t* lvfpix, const uint16_t* rvfpix)
{
    int ret = wd&(~(Upsample2to1by14641SIMDBatch-1)); // number of dest pixels completed
    const Byte* pstlv = pst + (ret/2)*4;

#if (defined(_M_IX86) || defined(_M_AMD64))
    __m128i zero = _mm_setzero_si128();
#endif
    __m128i x0,x1,x2,x3;
    __m128i L,C,R;

    // load L(eft) pix from left pix buffer
    L = _mm_load_si128((const __m128i*)lvfpix);

    // load and vertically-filter two RGBA pixels for C(enter)
#if (defined(_M_IX86) || defined(_M_AMD64))
    x0 = _mm_loadl_epi64((const __m128i*)pst); pst += (4*2);
    x1 = _mm_loadl_epi64((const __m128i*)psc); psc += (4*2);
    x2 = _mm_loadl_epi64((const __m128i*)psb); psb += (4*2);
    x0 = _mm_unpacklo_epi8(x0,zero);
    x1 = _mm_unpacklo_epi8(x1,zero);
    x2 = _mm_unpacklo_epi8(x2,zero);
    x0 = _mm_adds_epu16(x0,x2); // L+R
    x1 = _mm_slli_epi16(x1,1); 
    x0 = _mm_adds_epu16(x0,x1); // += (C<<1)
    x1 = _mm_slli_epi16(x1,1); 
    C = _mm_adds_epu16(x0,x1); // += (C<<2) is 8.3 vertically filtered
#elif (defined(_M_ARM))
    {
        __n64 x0n = vld1_u8(pst); pst += (4*2);
        __n64 x1n = vld1_u8(psc); psc += (4*2);
        __n64 x2n = vld1_u8(psb); psb += (4*2);
        x0 = vaddl_u8(x0n,x2n); // L+R
        x1 = vshll_n_u8(x1n,1); // 16bpp C << 1
    }
    x0 = vqaddq_u16(x0,x1);
    x1 = vshlq_n_u16(x1,1); // C << (1+1)
    C = vqaddq_u16(x0,x1);
#endif

    while ( wd >= Upsample2to1by14641SIMDBatch )
    {
        if (pst == pstlv)
        {
            R = _mm_load_si128((const __m128i*)rvfpix);
        }
        else
        {
            // load and vertically-filter two RGBA pixels for R(ight)
#if (defined(_M_IX86) || defined(_M_AMD64))
            x0 = _mm_loadl_epi64((const __m128i*)pst); pst += (4*2); 
            x1 = _mm_loadl_epi64((const __m128i*)psc); psc += (4*2); 
            x2 = _mm_loadl_epi64((const __m128i*)psb); psb += (4*2); 
            x0 = _mm_unpacklo_epi8(x0,zero);
            x1 = _mm_unpacklo_epi8(x1,zero);
            x2 = _mm_unpacklo_epi8(x2,zero);
            x0 = _mm_adds_epu16(x0,x2); // L+R
            x1 = _mm_slli_epi16(x1,1); 
            x0 = _mm_adds_epu16(x0,x1); // += (C<<1)
            x1 = _mm_slli_epi16(x1,1); 
            R = _mm_adds_epu16(x0,x1); // += (C<<2) is 8.3 vertically filtered
#elif (defined(_M_ARM))
            {
                __n64 x0n = vld1_u8(pst); pst += (4*2);
                __n64 x1n = vld1_u8(psc); psc += (4*2);
                __n64 x2n = vld1_u8(psb); psb += (4*2);
                x0 = vaddl_u8(x0n,x2n); // L+R
                x1 = vshll_n_u8(x1n,1); // 16bpp C << 1
            }
            x0 = vqaddq_u16(x0,x1);
            x1 = vshlq_n_u16(x1,1); // C << (1+1)
            R = vqaddq_u16(x0,x1);
#endif
        }

        // same horizontal filtering as above
#if (defined(_M_IX86) || defined(_M_AMD64))
        x3 = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(L),_mm_castsi128_pd(C),1)); // L1C0
        x2 = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(C),_mm_castsi128_pd(R),1)); // C1R0
#elif (defined(_M_ARM))
        x3 = vextq_u16(L,C,4);
        x2 = vextq_u16(C,R,4);
#endif
        x1 = _mm_adds_epu16(x2,x3); // L1C0+C1R0
        x0 = _mm_slli_epi16(C,2);
        x1 = _mm_adds_epu16(x1,x0); // += (C0C1)<<2
        x0 = _mm_slli_epi16(C,1);
        x0 = _mm_adds_epu16(x1,x0); // += (C0C1)<<1 is unnormalized R0R2 in 8.6
        x1 = _mm_adds_epu16(C,x2);  // is unnormalized R1R3 in 8.4
#if (defined(_M_IX86) || defined(_M_AMD64))
        x0 = _mm_adds_epu16(x0,_mm_set1_epi16(1<<5));
        x0 = _mm_srli_epi16(x0,6); // normalized R0R2
        x1 = _mm_adds_epu16(x1,_mm_set1_epi16(1<<3));
        x1 = _mm_srli_epi16(x1,4); // normalized R1R3
        x0 = _mm_packus_epi16(x0,x1); // Byte R0R2R1R3
        x0 = _mm_shuffle_epi32(x0,_MM_SHUFFLE(3,1,2,0)); // R0R1R2R3
#elif (defined(_M_ARM))
        {
            uint8x8x2_t x0t;
            x0t.val[0] = vrshrn_n_u16(x0,6); // normalized and packed R0R2
            x0t.val[1] = vrshrn_n_u16(x1,4); // normalized and packed R1R3
            x0t = vzip_u32(x0t.val[0],x0t.val[1]); // R0R1, R2R3
            x0 = vcombine_u64(x0t.val[0],x0t.val[1]); // R0R1R2R3
        }
#endif
        pfnStore((__m128i*)pd,x0); pd += (4*4);

        L = C;
        C = R;
        wd -= Upsample2to1by14641SIMDBatch;
    }

    return ret; 
}
#endif

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
static void GetKernelSetPosition(OUT UInt32& uIndexInCycle, OUT int& iSrcStart, 
						  int iDstStart, const C1dKernelSet& ks)
{
	// this function divides the supplied dst coord by the cycle count and
	// rounds to neg-infinity to find the first coordinate in the cycle
	int iCycleLen   = (int)ks.GetCycle();
	int iCycle      = iDstStart / iCycleLen;
	int iCycleStart = iCycle * iCycleLen;
	int iCycleRem   = iDstStart - iCycleStart;
	if( iCycleRem < 0 )
	{
		uIndexInCycle = (UInt32)(iCycleRem + iCycleLen);
		iSrcStart     = (iCycle-1) * ks.GetCoordShiftPerCycle();
	}
	else
	{
		uIndexInCycle = (UInt32)iCycleRem;
		iSrcStart     = iCycle * ks.GetCoordShiftPerCycle();
	}
}

template <typename Tio, typename Ttmp>
static HRESULT 
SeparableFilter14641Upsample2to1(CImg& imgDst, const CRect& rctDst,
                    const CImg& imgSrc, CPoint ptSrcOrigin,
                    C1dKernelSet& m_ksUp)
{
    VT_HR_BEGIN();

    // ASSERT that enough source data was provided
    CRect rctSrc    = vt::CRect(imgSrc.Rect()) + ptSrcOrigin;
    CRect rctReqSrc = GetRequiredSrcRect(rctDst,m_ksUp,m_ksUp);
    VT_ASSERT( rctSrc.RectInRect(&rctReqSrc) );

    UInt32 uKernIndex;
    int  iSrcOffsetV, iSrcOffsetH;
    GetKernelSetPosition(uKernIndex, iSrcOffsetH, rctDst.left, m_ksUp);
    GetKernelSetPosition(uKernIndex, iSrcOffsetV, rctDst.top, m_ksUp);

    // origin in imgSrc
    CPoint srcOrigin = CPoint( iSrcOffsetH - ptSrcOrigin.x, 
                               iSrcOffsetV - ptSrcOrigin.y);

    VT_HR_EXIT( (!imgSrc.IsValid()) ? E_INVALIDSRC : S_OK );
    // create the destination image if it needs to be
    VT_HR_EXIT( CreateImageForTransform(imgDst, rctDst.Width(), rctDst.Height(),
                                        imgSrc.GetType()) );

    int ys = srcOrigin.y;
    int ydr = rctDst.top;
    for (int yd=0; yd < rctDst.Height(); yd++,ydr++)
    {
        Ttmp pbufa[4], pbufb[4], pbufc[4];
        Tio* pd = (Tio*)imgDst.BytePtr(yd);
        int wd = rctDst.Width();
        int xdr = rctDst.left;
        if (ydr & 0x1)
        {
            // dest is 'between' src pixels vertically, so vfilter with 'x4x4x' (average)

            // ptrs to two scanline sources
            const Tio* pyt = (Tio*)imgSrc.BytePtr(srcOrigin.x,ys);
            const Tio* pyb = (Tio*)imgSrc.BytePtr(srcOrigin.x,ys+1);
            // advance source address when 'between' source only
            ys++;

            // vector implementation
#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))
            if ( (xdr & 0x1) && (wd >= (Upsample2to1by14641SIMDBatch+1)) )
            {
                // first pixel is 'between' in x so do scalar since vector code only does first pixel 'on'
                Filter4byX4X4XP1(pbufa,pyt,pyb);
                Filter4byX4X4XP1(pbufb,pyt+4,pyb+4);
                Filter4byX4X4XP2(pd,pbufa,pbufb); pd+=4;
                pyt+=4; pyb+=4;
                wd -= 1; xdr++;
            }
            if ( (VT_UPSAMPLE14641_SIMD_SUPPORTED) && (wd > Upsample2to1by14641SIMDBatch))
            // TODO: figure out what is failing when vector code does only 4 dest 
            {
                // compute left and right vertically filtered pixels for batch
                VT_DECLSPEC_ALIGN(16) Ttmp lpix[8];
                Filter4byX4X4XP1(&lpix[4],pyt-4,pyb-4);
                VT_DECLSPEC_ALIGN(16) Ttmp rpix[8];
                int lastROffsV = (wd&(~(Upsample2to1by14641SIMDBatch-1)))*(4/2);
                Filter4byX4X4XP1(&rpix[0],pyt+lastROffsV,pyb+lastROffsV);

                int wdv;
                SelectPSTFunc(IsAligned16(pd),wdv,Store4,
                    Upsample2to1by14641ARowSIMD,
                    pd,pyt,pyb,wd,lpix,rpix);

                wd -= wdv;
                pyt += ((wdv/2)*4);
                pyb += ((wdv/2)*4);
                pd += (wdv*4);
                xdr += wdv;
            }
#endif

            // scalar implementation/remainder
            for (int xd = rctDst.Width()-wd; xd < rctDst.Width(); xd++,xdr++)
            {
                if (xdr & 0x1)
                {
                    // 'between' src pixels horizontally and vertically, so
                    // filter two pixels vertically for horizontal inputs

                    // apply filter vertically for left
                    Filter4byX4X4XP1(pbufa,pyt,pyb);
                    Filter4byX4X4XP1(pbufb,pyt+4,pyb+4);

                    // apply horizontal filter and store
                    Filter4byX4X4XP2(pd,pbufa,pbufb); pd+=4;

                    // advance source ptrs for 'between' source only 
                    pyt+=4; pyb+=4;
                }
                else
                {
                    // 'on' src pixel horizontally and 'between' vertically,
                    // so filter three pixels vertically for horizontal inputs

                    // apply vertical filter
                    Filter4byX4X4XP1(pbufa,pyt-4,pyb-4);
                    Filter4byX4X4XP1(pbufb,pyt,pyb);
                    Filter4byX4X4XP1(pbufc,pyt+4,pyb+4);

                    // apply horizontal filter and store
                    Filter4by1X6X1P2(pd,pbufa,pbufb,pbufc); pd+=4;
                }
            }
        }
        else
        {
            // dest is 'on' source pixel vertically, so vfilter with '1x6x1'
            const Tio* pyt = (Tio*)imgSrc.BytePtr(srcOrigin.x,ys-1);
            const Tio* pyc = (Tio*)imgSrc.BytePtr(srcOrigin.x,ys);
            const Tio* pyb = (Tio*)imgSrc.BytePtr(srcOrigin.x,ys+1);

            // vector implementation
#if (defined(_M_IX86) || defined(_M_AMD64) || defined(_M_ARM))
            if ( (xdr & 0x1) && (wd >= (Upsample2to1by14641SIMDBatch+1)) )
            {
                // first pixel is 'between' in x so do scalar since vector code only does first pixel 'on'
                Filter4by1X6X1P1(pbufa,pyt,pyc,pyb);
                Filter4by1X6X1P1(pbufb,pyt+4,pyc+4,pyb+4);
                Filter4byX4X4XP2(pd,pbufa,pbufb); pd+=4;
                pyt+=4; pyc+=4; pyb+=4;
                wd -= 1; xdr++;
            }
            if ( (VT_UPSAMPLE14641_SIMD_SUPPORTED) && (wd > Upsample2to1by14641SIMDBatch))
            // TODO: figure out what is failing when vector code does only 4 dest 
            {
                // compute left and right vertically filtered pixels for batch
                VT_DECLSPEC_ALIGN(16) Ttmp lpix[8];
                Filter4by1X6X1P1(&lpix[4],pyt-4,pyc-4,pyb-4);
                VT_DECLSPEC_ALIGN(16) Ttmp rpix[8];
                int lastROffsV = (wd&(~(Upsample2to1by14641SIMDBatch-1)))*(4/2);
                Filter4by1X6X1P1(&rpix[0],pyt+lastROffsV,pyc+lastROffsV,pyb+lastROffsV);

                int wdv;
                SelectPSTFunc(IsAligned16(pd),wdv,Store4,
                    Upsample2to1by14641BRowSIMD,
                    pd,pyt,pyc,pyb,wd,lpix,rpix);

                wd -= wdv;
                pyt += ((wdv/2)*4);
                pyc += ((wdv/2)*4);
                pyb += ((wdv/2)*4);
                pd += (wdv*4);
                xdr += wdv;
            }
#endif

            // scalar implementation/remainder
            for (int xd = rctDst.Width()-wd; xd < rctDst.Width(); xd++,xdr++)
            {
                if (xdr & 0x1)
                {
                    // 'between' source pixels horizontally and 'on' vertically,
                    // so filter two pixels vertically for horizontal inputs
                    Filter4by1X6X1P1(pbufa,pyt,pyc,pyb);
                    Filter4by1X6X1P1(pbufb,pyt+4,pyc+4,pyb+4);

                    // apply horizontal filter and store
                    Filter4byX4X4XP2(pd,pbufa,pbufb); pd+=4;

                    // advance source ptrs for 'between' source only 
                    pyt+=4; pyc+=4; pyb+=4;
                }
                else
                {
                    // 'on' source pixel horizontally and 'on' vertically,
                    // so filter three pixels vertically for horizontal inputs
                    Filter4by1X6X1P1(pbufa,pyt-4,pyc-4,pyb-4);
                    Filter4by1X6X1P1(pbufb,pyt,pyc,pyb);
                    Filter4by1X6X1P1(pbufc,pyt+4,pyc+4,pyb+4);

                    // apply horizontal filter and store
                    Filter4by1X6X1P2(pd,pbufa,pbufb,pbufc); pd+=4;
                }
            }
        }
    }
   VT_HR_END();
}

//-----------------------------------------------------------------------------
// CSeparableFilter14641Transform routines
//+----------------------------------------------------------------------------

#ifndef VT_NO_XFORMS

HRESULT CSeparableFilter14641Transform::Clone(ITaskState **ppState)
{
	return CloneTaskState<CSeparableFilter14641Transform>(ppState, 
		[this](CSeparableFilter14641Transform* pN)
	{ return pN->InitializeUpsample2to1(m_dstType); });
}

HRESULT CSeparableFilter14641Transform::Transform(
    CImg* pimgDst, IN  const CRect& rctDst,
	const CImg& imgSrc, const CPoint& ptSrc)

{
    if (pimgDst->GetType() == OBJ_RGBAIMG)
    {
        return SeparableFilter14641Upsample2to1<Byte,uint16_t>(*pimgDst, rctDst,
            imgSrc, ptSrc, m_ksUp);
    }
    else if (pimgDst->GetType() == OBJ_RGBAFLOATIMG)
    {
        HRESULT hr =  SeparableFilter14641Upsample2to1<float,float>(*pimgDst, rctDst,
            imgSrc, ptSrc, m_ksUp);
        return hr;
    }
    else
    {
        return E_INVALIDARG;
    }
}

HRESULT CSeparableFilter14641Transform::InitializeUpsample2to1(int dstType)
{
    VT_HR_BEGIN();

    if ( (dstType != VT_IMG_FIXED(OBJ_RGBAIMG)) && (dstType != VT_IMG_FIXED(OBJ_RGBAFLOATIMG)) )
    {
        VT_HR_EXIT(E_INVALIDARG);
    }
    m_dstType = dstType;

    float arfKrnlUp0[] = {1.f/8.f, 6.f/8.f, 1.f/8.f};
	float arfKrnlUp1[] = {4.f/8.f, 4.f/8.f};
	VT_HR_EXIT( m_ksUp.Create(2, 1) );
	VT_HR_EXIT( m_ksUp.Set(0, -1, 3, arfKrnlUp0) );
	VT_HR_EXIT( m_ksUp.Set(1, 0, 2, arfKrnlUp1) );

    VT_HR_END();
}

#endif

//-----------------------------------------------------------------------------
// end