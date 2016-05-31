//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      specialized implementation of Harris feature detector
//
//      SSE specific code
//
//------------------------------------------------------------------------------
#pragma once

//------------------------------------------------------------------------------
inline void HarrisDoDrvLineIntSSE(int& w, int16_t*& res, size_t& resStride,
    const uint8_t*& srcP, const uint8_t*& srcC, const uint8_t*& srcN)
{
    if ( g_SupportSSE2() )
    {
        __m128i zeroReg = _mm_setzero_si128();
        while (w > 8)
        {
            // load x,y for previous and next
            __m128i pix8xp = _mm_loadl_epi64((__m128i*)(srcC-2));
            __m128i pix8xn = _mm_loadl_epi64((__m128i*)(srcC)); srcC+=8;
            __m128i pix8yp = _mm_loadl_epi64((__m128i*)(srcP)); srcP+=8;
            __m128i pix8yn = _mm_loadl_epi64((__m128i*)(srcN)); srcN+=8;
           // unpack 8 pixels to 16bpp in registers
           pix8xp = _mm_unpacklo_epi8(pix8xp,zeroReg);
           pix8xn = _mm_unpacklo_epi8(pix8xn,zeroReg);
           pix8yp = _mm_unpacklo_epi8(pix8yp,zeroReg);
           pix8yn = _mm_unpacklo_epi8(pix8yn,zeroReg);
           // differences
           __m128i dx8 = _mm_sub_epi16(pix8xp,pix8xn);
           __m128i dy8 = _mm_sub_epi16(pix8yn,pix8yp);
           // squared (unsigned result) differences >> 1
           __m128i xx8 = _mm_mullo_epi16(dx8,dx8); 
           xx8 = _mm_srli_epi16(xx8,1);
           __m128i yy8 = _mm_mullo_epi16(dy8,dy8); 
           yy8 = _mm_srli_epi16(yy8,1);
           // shift one of the inputs before multiply for (signed result) xy 
           dx8 = _mm_srai_epi16(dx8,1);
           __m128i xy8 = _mm_mullo_epi16(dx8,dy8);
           // store results
           _mm_store_si128((__m128i*)(res+(0*resStride)),xx8);
           _mm_store_si128((__m128i*)(res+(1*resStride)),yy8);
           _mm_store_si128((__m128i*)(res+(2*resStride)),xy8);
           res += 8;
           w -= 8;
        }
    }
}

inline void HarrisDoDrvLineUIntSSE(int& w, uint16_t*& res, size_t& resStride,
    const uint8_t*& srcP, const uint8_t*& srcC, const uint8_t*& srcN)
{
    if ( g_SupportSSE2() )
    {
        __m128i zeroReg = _mm_setzero_si128();
        __m128i halfReg = _mm_set1_epi16(-32768); // 0x8000
        while (w > 8)
        {
            // load x,y for previous and next
            __m128i pix8xp = _mm_loadl_epi64((__m128i*)(srcC-2));
            __m128i pix8xn = _mm_loadl_epi64((__m128i*)(srcC)); srcC+=8;
            __m128i pix8yp = _mm_loadl_epi64((__m128i*)(srcP)); srcP+=8;
            __m128i pix8yn = _mm_loadl_epi64((__m128i*)(srcN)); srcN+=8;
           // unpack 8 pixels to 16bpp in registers
           pix8xp = _mm_unpacklo_epi8(pix8xp,zeroReg);
           pix8xn = _mm_unpacklo_epi8(pix8xn,zeroReg);
           pix8yp = _mm_unpacklo_epi8(pix8yp,zeroReg);
           pix8yn = _mm_unpacklo_epi8(pix8yn,zeroReg);
           // differences
           __m128i dx8 = _mm_sub_epi16(pix8xp,pix8xn);
           __m128i dy8 = _mm_sub_epi16(pix8yn,pix8yp);
           // squared (unsigned result) differences
           __m128i xx8 = _mm_mullo_epi16(dx8,dx8); 
           __m128i yy8 = _mm_mullo_epi16(dy8,dy8); 
           // shift one input down so signed multiply result fits in int16
           dx8 = _mm_srai_epi16(dx8,1);
           // multiply result is s.15, then bias by 1/2
           __m128i xy8 = _mm_mullo_epi16(dx8,dy8);
           xy8 = _mm_add_epi16(xy8,halfReg);
           // store results
           _mm_store_si128((__m128i*)(res+(0*resStride)),xx8);
           _mm_store_si128((__m128i*)(res+(1*resStride)),yy8);
           _mm_store_si128((__m128i*)(res+(2*resStride)),xy8);
           res += 8;
           w-= 8;
        }
    }
}

//------------------------------------------------------------------------------
// SSE for 5x5 Harris filtering
//------------------------------------------------------------------------------
inline void HarrisFilterVert5SSE(int& w, int srcStride,
    uint16_t*& dst0, uint16_t*& dst1, uint16_t*& dst2, 
    uint16_t*& src0, uint16_t*& src1, uint16_t*& src2, 
    uint16_t*& src3, uint16_t*& src4)
{
    if ( g_SupportSSE2() )
    {
        while (w > 8)
        {
            __m128i acc0;
            __m128i acc1;
            __m128i acc2;
            __m128i src;

            // process first sample and assign to accumulator registers
            src = _mm_load_si128((__m128i*)src0);
            acc0 = _mm_srli_epi16(src,4);
            src = _mm_load_si128((__m128i*)(src0+srcStride));
            acc1 = _mm_srli_epi16(src,4);
            src = _mm_load_si128((__m128i*)(src0+(2*srcStride))); src0 += 8;
            acc2 = _mm_srli_epi16(src,4);

#define ACC_SINGLE_SHIFT(_i_,_sh_,_d_) \
    src = _mm_load_si128((__m128i*)(src##_i_##)); \
    src = _mm_srli_epi16(src,_sh_); \
    acc0 = _mm_add_epi16(acc0,src); \
    src = _mm_load_si128((__m128i*)(src##_i_##+srcStride)); \
    src = _mm_srli_epi16(src,_sh_); \
    acc1 = _mm_add_epi16(acc1,src); \
    src = _mm_load_si128((__m128i*)(src##_i_##+(2*srcStride))); src##_i_## += _d_; \
    src = _mm_srli_epi16(src,_sh_); \
    acc2 = _mm_add_epi16(acc2,src); \

            ACC_SINGLE_SHIFT(1,2,8);
            ACC_SINGLE_SHIFT(2,3,0);
            ACC_SINGLE_SHIFT(2,2,8); // TODO: could avoid reloading...
            ACC_SINGLE_SHIFT(3,2,8);
            ACC_SINGLE_SHIFT(4,4,8);
#undef ACC_SINGLE_SHIFT

            _mm_store_si128((__m128i*)dst0, acc0); dst0 += 8;
            _mm_store_si128((__m128i*)dst1, acc1); dst1 += 8;
            _mm_store_si128((__m128i*)dst2, acc2); dst2 += 8;
            w-=8;
        }
    }
}

#define FILTER8_HORIZ_5(_src_,_res_) \
{ \
    int32_t v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,vA,vB; \
    v0 = *(_src_+0); v1 = *(_src_+1); v2 = *(_src_+2); v3 = *(_src_+3); v4 = *(_src_+4); \
    v5 = *(_src_+5); v6 = *(_src_+6); v7 = *(_src_+7); v8 = *(_src_+8); v9 = *(_src_+9); \
    vA = *(_src_+10); vB = *(_src_+11); \
    _res_[0] = (v0>>4)+(v1>>2)+(v2>>3)+(v2>>2)+(v3>>2)+(v4>>4); \
    _res_[1] = (v1>>4)+(v2>>2)+(v3>>3)+(v3>>2)+(v4>>2)+(v5>>4); \
    _res_[2] = (v2>>4)+(v3>>2)+(v4>>3)+(v4>>2)+(v5>>2)+(v6>>4); \
    _res_[3] = (v3>>4)+(v4>>2)+(v5>>3)+(v5>>2)+(v6>>2)+(v7>>4); \
    _res_[4] = (v4>>4)+(v5>>2)+(v6>>3)+(v6>>2)+(v7>>2)+(v8>>4); \
    _res_[5] = (v5>>4)+(v6>>2)+(v7>>3)+(v7>>2)+(v8>>2)+(v9>>4); \
    _res_[6] = (v6>>4)+(v7>>2)+(v8>>3)+(v8>>2)+(v9>>2)+(vA>>4); \
    _res_[7] = (v7>>4)+(v8>>2)+(v9>>3)+(v9>>2)+(vA>>2)+(vB>>4); \
}
inline void HarrisFilterHoriz5SSE(int& w, float*& dst, 
    uint16_t*& buf0, uint16_t*& buf1, uint16_t*& buf2)
{
    if ( g_SupportSSE2() )
    {
        // horizontal filter and compute 8 harris function pixels at a time to enable
        // reuse of read data on the vertical filtering and to use SSE
        __m128i zeroReg = _mm_setzero_si128();
        __m128i halfReg = _mm_set1_epi32(0x10000);
        while (w > 8)
        {
            int16_t xx[8],yy[8],xy[8];
            FILTER8_HORIZ_5(buf0,xx); buf0 += 8;
            FILTER8_HORIZ_5(buf1,yy); buf1 += 8;
            FILTER8_HORIZ_5(buf2,xy); buf2 += 8;

            __m128i ixx8 = _mm_loadu_si128((__m128i*)xx);
            __m128i iyy8 = _mm_loadu_si128((__m128i*)yy);
            __m128i ixy8 = _mm_loadu_si128((__m128i*)xy);

            __m128i ixx4 = _mm_unpacklo_epi16(ixx8,zeroReg); 
            __m128i iyy4 = _mm_unpacklo_epi16(iyy8,zeroReg);
            __m128i ixy4 = _mm_unpacklo_epi16(ixy8,zeroReg); 
            ixy4 = _mm_slli_epi32(ixy4,1); ixy4 = _mm_sub_epi32(ixy4,halfReg);
                        
            __m128 fxx4 = _mm_cvtepi32_ps(ixx4);
            __m128 fyy4 = _mm_cvtepi32_ps(iyy4);
            __m128 fxy4 = _mm_cvtepi32_ps(ixy4);
            __m128 fHarris4 = _mm_sub_ps(_mm_mul_ps(fxx4,fyy4),_mm_mul_ps(fxy4,fxy4));
            _mm_store_ps(dst,fHarris4); dst+=4;

            ixx4 = _mm_unpackhi_epi16(ixx8,zeroReg); 
            iyy4 = _mm_unpackhi_epi16(iyy8,zeroReg);
            ixy4 = _mm_unpackhi_epi16(ixy8,zeroReg); 
            ixy4 = _mm_slli_epi32(ixy4,1); ixy4 = _mm_sub_epi32(ixy4,halfReg);

                        
            fxx4 = _mm_cvtepi32_ps(ixx4);
            fyy4 = _mm_cvtepi32_ps(iyy4);
            fxy4 = _mm_cvtepi32_ps(ixy4);
            fHarris4 = _mm_sub_ps(_mm_mul_ps(fxx4,fyy4),_mm_mul_ps(fxy4,fxy4));
            _mm_store_ps(dst,fHarris4); dst+=4;

            w -= 8;
        }
    }
}
#undef FILTER8_HORIZ_5
//------------------------------------------------------------------------------
    
    
//------------------------------------------------------------------------------
// SSE for 7x7 Harris filtering
//------------------------------------------------------------------------------
inline void HarrisFilterVert7SSE(int& w, int srcStride,
    int32_t*& dst0, int32_t*& dst1, int32_t*& dst2, 
    int16_t*& src0, int16_t*& src1, int16_t*& src2, 
    int16_t*& src3, int16_t*& src4, int16_t*& src5, int16_t*& src6)
{
    if ( g_SupportSSE2() )
    {
        __m128i zeroReg = _mm_setzero_si128();
        // TODO: just do 4 at a time here?
        while (w > 8)
        {
            __m128i acc0l,acc0h;
            __m128i acc1l,acc1h;
            __m128i acc2l,acc2h;
            __m128i src;

            // process first sample and assign to accumulator registers
            src = _mm_load_si128((__m128i*)src0);
            acc0l = _mm_unpacklo_epi16(zeroReg,src); acc0l = _mm_srai_epi32(acc0l,16-0);
            acc0h = _mm_unpackhi_epi16(zeroReg,src); acc0h = _mm_srai_epi32(acc0h,16-0);
            src = _mm_load_si128((__m128i*)(src0+srcStride));
            acc1l = _mm_unpacklo_epi16(zeroReg,src); acc1l = _mm_srai_epi32(acc1l,16-0);
            acc1h = _mm_unpackhi_epi16(zeroReg,src); acc1h = _mm_srai_epi32(acc1h,16-0);
            src = _mm_load_si128((__m128i*)(src0+(2*srcStride))); src0 += 8;
            acc2l = _mm_unpacklo_epi16(zeroReg,src); acc2l = _mm_srai_epi32(acc2l,16-0);
            acc2h = _mm_unpackhi_epi16(zeroReg,src); acc2h = _mm_srai_epi32(acc2h,16-0);

            __m128i srcl, srch;
#define ACC_SINGLE_SHIFT(_i_,_sh_,_d_) \
    src = _mm_load_si128((__m128i*)(src##_i_##)); \
    srcl = _mm_unpacklo_epi16(zeroReg,src); srcl = _mm_srai_epi32(srcl, _sh_); \
    srch = _mm_unpackhi_epi16(zeroReg,src); srch = _mm_srai_epi32(srch, _sh_); \
    acc0l = _mm_add_epi32(acc0l,srcl); acc0h = _mm_add_epi32(acc0h,srch); \
    src = _mm_load_si128((__m128i*)(src##_i_##+srcStride)); \
    srcl = _mm_unpacklo_epi16(zeroReg,src); srcl = _mm_srai_epi32(srcl, _sh_); \
    srch = _mm_unpackhi_epi16(zeroReg,src); srch = _mm_srai_epi32(srch, _sh_); \
    acc1l = _mm_add_epi32(acc1l,srcl); acc1h = _mm_add_epi32(acc1h,srch); \
    src = _mm_load_si128((__m128i*)(src##_i_##+(2*srcStride))); src##_i_## += _d_; \
    srcl = _mm_unpacklo_epi16(zeroReg,src); srcl = _mm_srai_epi32(srcl, _sh_); \
    srch = _mm_unpackhi_epi16(zeroReg,src); srch = _mm_srai_epi32(srch, _sh_); \
    acc2l = _mm_add_epi32(acc2l,srcl); acc2h = _mm_add_epi32(acc2h,srch); 

            ACC_SINGLE_SHIFT(1,16-2,8);
            ACC_SINGLE_SHIFT(2,16-3,8);
            ACC_SINGLE_SHIFT(3,16-3,0);
            ACC_SINGLE_SHIFT(3,16-1,8); // TODO: could avoid reloading...
            ACC_SINGLE_SHIFT(4,16-3,8);
            ACC_SINGLE_SHIFT(5,16-2,8);
            ACC_SINGLE_SHIFT(6,16-0,8);
#undef ACC_SINGLE_SHIFT

            _mm_store_si128((__m128i*)dst0, acc0l); dst0 += 4;
            _mm_store_si128((__m128i*)dst0, acc0h); dst0 += 4;
            _mm_store_si128((__m128i*)dst1, acc1l); dst1 += 4;
            _mm_store_si128((__m128i*)dst1, acc1h); dst1 += 4;
            _mm_store_si128((__m128i*)dst2, acc2l); dst2 += 4;
            _mm_store_si128((__m128i*)dst2, acc2h); dst2 += 4;
            w-=8;
        }
    }
}

// |A B C D|E F G H|I J - - |
// result is (A)+(B<<2)+(C<<3)+(D<<3)+(D<<1)+(E<<3)+(F<<2)+(G)  13 ops
// = A+((((B+((C+D+E)<<1)+F)<<1)+D)<<1)+G   10 ops
// = ((B+((C+D+E)<<1)+F)<<2) + (A+(D<<1)+G) same 10 ops with fewer intermediate registers
// notation: A# is vector |A B C D|, B# is vector |B C D E|, etc.
#define FILTER4_HORIZ_7_SSE(_src_,_res_) \
{ \
    __m128i resr4 = _mm_load_si128((__m128i*)(_src_+0)); /* |A B C D| */ \
    __m128i midl4 = _mm_load_si128((__m128i*)(_src_+4)); /* |E F G H| */ \
    __m128i resl4, algn4, temp4; \
    algn4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(resr4), \
                _mm_castsi128_ps(midl4),_MM_SHUFFLE(1,0,3,2))); /* C# */ \
    resl4 = _mm_add_epi32(midl4,algn4); /* resl = E#+C# */ \
    algn4 = _mm_srli_si128(resr4,4*3);  /* |D 0 0 0| */ \
    temp4 = _mm_slli_si128(midl4,4*1);  /* |0 E F G| */ \
    temp4 = _mm_or_si128(algn4,temp4);  /* D# */ \
    resl4 = _mm_add_epi32(resl4,temp4); /* resl += D# */ \
    resl4 = _mm_slli_epi32(resl4,1);    /* resl = (C#+D#+E#)<<1 */ \
    temp4 = _mm_slli_epi32(temp4,1);    /* D#<<1 */ \
    algn4 = _mm_srli_si128(resr4,4*1);  /* |B C D 0| */ \
    resr4 = _mm_add_epi32(resr4,temp4); /* resr = A#+(D#<<1) */ \
    temp4 = _mm_slli_si128(midl4,4*3);  /* |0 0 0 E| */ \
    algn4 = _mm_or_si128(algn4,temp4);  /* B# */ \
    resl4 = _mm_add_epi32(resl4,algn4); /* resl += B# */ \
    temp4 = _mm_load_si128((__m128i*)(_src_+8)); /* |I J - -| */ \
    algn4 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(midl4), \
                _mm_castsi128_ps(temp4),_MM_SHUFFLE(1,0,3,2))); /* G# */ \
    resr4 = _mm_add_epi32(resr4,algn4); /* resr = A#+(D#<<1)+G# */ \
    algn4 = _mm_srli_si128(midl4,4*1);  /* |F G H 0| */ \
    temp4 = _mm_slli_si128(temp4,4*3);  /* |0 0 0 I| */ \
    algn4 = _mm_or_si128(algn4,temp4);  /* F# */ \
    resl4 = _mm_add_epi32(resl4,algn4); /* resl += F# */ \
    resl4 = _mm_slli_epi32(resl4,2);    /* (B#+((C#+D#+E#)<<1)+F#)<<2 */ \
    resl4 = _mm_add_epi32(resl4,resr4); /* result */ \
    _res_ = _mm_cvtepi32_ps(resl4); \
} // 22 ops

inline void HarrisFilterHoriz7SSE(int& w, float*& dst, 
    int32_t*& buf0, int32_t*& buf1, int32_t*& buf2,
    float  fscaleSq)
{
    if ( g_SupportSSE2() )
    {
        const __m128 fScaleSq4 = _mm_set_ps1(fscaleSq);
        while (w > 4)
        {
            __m128 fxx4,fyy4,fxy4;
            FILTER4_HORIZ_7_SSE(buf0,fxx4); buf0 += 4;
            FILTER4_HORIZ_7_SSE(buf1,fyy4); buf1 += 4;
            FILTER4_HORIZ_7_SSE(buf2,fxy4); buf2 += 4;
            __m128 fHarris4 = _mm_sub_ps(_mm_mul_ps(fxx4,fyy4),_mm_mul_ps(fxy4,fxy4));
            fHarris4 = _mm_mul_ps(fHarris4,fScaleSq4);
            _mm_store_ps(dst,fHarris4); dst+=4;
            w -= 4;
        }
    }
}
#undef FILTER4_HORIZ_7_SSE
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// end