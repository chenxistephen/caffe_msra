//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      specialized implementation of Harris feature detector - 
//
//      Arm/Neon specific code
//
//------------------------------------------------------------------------------
#pragma once

#if !defined(vreinterpretq_u16_i16)
#define vreinterpretq_u16_i16 vreinterpretq_u16_s16
#endif
#if !defined(vreinterpretq_i16_u16)
#define vreinterpretq_i16_u16 vreinterpretq_s16_u16
#endif
#if !defined(vreinterpretq_i32_u32)
#define vreinterpretq_i32_u32 vreinterpretq_s32_u32
#endif

//------------------------------------------------------------------------------
// Neon for Harris difference computations
//------------------------------------------------------------------------------
inline void HarrisDoDrvLineIntNeon(int& w, int16_t*& res, size_t& resStride,
    const uint8_t*& srcP, const uint8_t*& srcC, const uint8_t*& srcN)
{
    while (w > 8)
    {
        //
        int16x8_t pix8xp = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcC-2)));
        int16x8_t pix8xn = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcC))); srcC+=8;
        int16x8_t pix8yp = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcP))); srcP+=8;
        int16x8_t pix8yn = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcN))); srcN+=8;
        //
        int16x8_t dx8 = vsubq_s16(pix8xp,pix8xn);
        int16x8_t dy8 = vsubq_s16(pix8yn,pix8yp);
        //
        int16x8_t xx8 = vmulq_s16(dx8,dx8);
        xx8 = vreinterpretq_i16_u16(vshrq_n_u16(vreinterpretq_u16_i16(xx8),1));
        int16x8_t yy8 = vmulq_s16(dy8,dy8);
        yy8 = vreinterpretq_i16_u16(vshrq_n_u16(vreinterpretq_u16_i16(yy8),1));
        //
        dx8 = vshrq_n_s16(dx8,1);
        int16x8_t xy8 = vmulq_s16(dx8,dy8);
        //
        vst1q_s16(res+(0*resStride),xx8);
        vst1q_s16(res+(1*resStride),yy8);
        vst1q_s16(res+(2*resStride),xy8);
        res += 8;
        w -= 8;
    }
}
inline void HarrisDoDrvLineUIntNeon(int& w, uint16_t*& res, size_t& resStride,
    const uint8_t*& srcP, const uint8_t*& srcC, const uint8_t*& srcN)
{
    int16x8_t halfReg = vmovq_n_s16((short)-32768); // 0x8000
    while (w > 8)
    {
        //
        int16x8_t pix8xp = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcC-2)));
        int16x8_t pix8xn = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcC))); srcC+=8;
        int16x8_t pix8yp = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcP))); srcP+=8;
        int16x8_t pix8yn = vreinterpretq_i16_u16(vmovl_u8(vld1_u8(srcN))); srcN+=8;
        //
        int16x8_t dx8 = vsubq_s16(pix8xp,pix8xn);
        int16x8_t dy8 = vsubq_s16(pix8yn,pix8yp);
        //
        int16x8_t xx8 = vmulq_s16(dx8,dx8);
        int16x8_t yy8 = vmulq_s16(dy8,dy8);
        //
        dx8 = vshrq_n_s16(dx8,1);
        int16x8_t xy8 = vmulq_s16(dx8,dy8);
        xy8 = vaddq_s16(xy8,halfReg);
        //
        vst1q_s16((int16_t*)res+(0*resStride),xx8);
        vst1q_s16((int16_t*)res+(1*resStride),yy8);
        vst1q_s16((int16_t*)res+(2*resStride),xy8);
        res += 8;
        w -= 8;
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Neon for 5x5 Harris filtering
//------------------------------------------------------------------------------
inline void HarrisFilterVert5Neon(int& w, int srcStride,
    uint16_t*& dst0, uint16_t*& dst1, uint16_t*& dst2, 
    uint16_t*& src0, uint16_t*& src1, uint16_t*& src2, 
    uint16_t*& src3, uint16_t*& src4)
{
    while (w > 8)
    {
        uint16x8_t acc0,acc1,acc2;
        // process first sample and assign to accumulator registers
        acc0 = vld1q_u16(src0+(0*srcStride));
        acc1 = vld1q_u16(src0+(1*srcStride));
        acc2 = vld1q_u16(src0+(2*srcStride));
        acc0 = vshrq_n_u16(acc0,4);
        acc1 = vshrq_n_u16(acc1,4);
        acc2 = vshrq_n_u16(acc2,4);
        src0 += 8;

        uint16x8_t t0,t1,t2;
#define ACC_LOAD(_i_) \
        t0 = vld1q_u16(src##_i_##+(0*srcStride)); \
        t1 = vld1q_u16(src##_i_##+(1*srcStride)); \
        t2 = vld1q_u16(src##_i_##+(2*srcStride)); \
        src##_i_## += 8;
#define ACC_SHIFT(_sh_) \
        t0 = vshrq_n_u16(t0,_sh_); \
        t1 = vshrq_n_u16(t1,_sh_); \
        t2 = vshrq_n_u16(t2,_sh_);
#define ACC_ACCUM() \
        acc0 = vaddq_u16(acc0,t0); \
        acc1 = vaddq_u16(acc1,t1); \
        acc2 = vaddq_u16(acc2,t2);

        ACC_LOAD(1); ACC_SHIFT(2); ACC_ACCUM();
        ACC_LOAD(2); ACC_SHIFT(2); ACC_ACCUM()
                        ACC_SHIFT(1); ACC_ACCUM() // ((>>2)>>1) = >>3
        ACC_LOAD(3); ACC_SHIFT(2); ACC_ACCUM()
        ACC_LOAD(4); ACC_SHIFT(4); ACC_ACCUM()
#undef ACC_LOAD
#undef ACC_SHIFT
#undef ACC_ACCUM
        vst1q_u16(dst0, acc0); dst0 += 8;
        vst1q_u16(dst1, acc1); dst1 += 8;
        vst1q_u16(dst2, acc2); dst2 += 8;
        w-=8;
    }
}

#define FILTER8_HORIZ_5(_src_,_res_) \
{ \
    int32_t v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,vA,vB; \
    v0 = *(_src_+0); v1 = *(_src_+1); v2 = *(_src_+2); v3 = *(_src_+3); v4 = *(_src_+4); \
    v5 = *(_src_+5); v6 = *(_src_+6); v7 = *(_src_+7); v8 = *(_src_+8); v9 = *(_src_+9); \
    vA = *(_src_+10); vB = *(_src_+11); \
    _res_[0] = (int16_t)((v0>>4)+(v1>>2)+(v2>>3)+(v2>>2)+(v3>>2)+(v4>>4)); \
    _res_[1] = (int16_t)((v1>>4)+(v2>>2)+(v3>>3)+(v3>>2)+(v4>>2)+(v5>>4)); \
    _res_[2] = (int16_t)((v2>>4)+(v3>>2)+(v4>>3)+(v4>>2)+(v5>>2)+(v6>>4)); \
    _res_[3] = (int16_t)((v3>>4)+(v4>>2)+(v5>>3)+(v5>>2)+(v6>>2)+(v7>>4)); \
    _res_[4] = (int16_t)((v4>>4)+(v5>>2)+(v6>>3)+(v6>>2)+(v7>>2)+(v8>>4)); \
    _res_[5] = (int16_t)((v5>>4)+(v6>>2)+(v7>>3)+(v7>>2)+(v8>>2)+(v9>>4)); \
    _res_[6] = (int16_t)((v6>>4)+(v7>>2)+(v8>>3)+(v8>>2)+(v9>>2)+(vA>>4)); \
    _res_[7] = (int16_t)((v7>>4)+(v8>>2)+(v9>>3)+(v9>>2)+(vA>>2)+(vB>>4)); \
}
inline void HarrisFilterHoriz5Neon(int& w, float*& dst, 
    uint16_t*& buf0, uint16_t*& buf1, uint16_t*& buf2)
{
    const int32x4_t halfReg = vmovq_n_s32(0x10000);
    while (w > 8)
    {
        int16_t xx[8],yy[8],xy[8];
        FILTER8_HORIZ_5(buf0,xx); buf0 += 8;
        FILTER8_HORIZ_5(buf1,yy); buf1 += 8;
        FILTER8_HORIZ_5(buf2,xy); buf2 += 8;

        for (int i=0; i<8; i+=4)
        {
            int32x4_t ixx4 = vreinterpretq_i32_u32(vmovl_u16(vld1_u16((const uint16_t*)&xx[i])));
            int32x4_t iyy4 = vreinterpretq_i32_u32(vmovl_u16(vld1_u16((const uint16_t*)&yy[i])));
            int32x4_t ixy4 = vreinterpretq_i32_u32(vmovl_u16(vld1_u16((const uint16_t*)&xy[i])));
            ixy4 = vshlq_n_s32(ixy4,1);
            ixy4 = vsubq_s32(ixy4,halfReg);

            float32x4_t fxx4 = vcvtq_f32_s32(ixx4);
            float32x4_t fyy4 = vcvtq_f32_s32(iyy4);
            float32x4_t fxy4 = vcvtq_f32_s32(ixy4);

            float32x4_t fHarris4 = vmulq_f32(fxx4,fyy4); 
            fHarris4 = vmlsq_f32(fHarris4,fxy4,fxy4);
            vst1q_f32(dst, fHarris4); dst+=4;

        }
        w -= 8;
    }
}
#undef FILTER8_HORIZ_5
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Neon for 7x7 Harris filtering
//------------------------------------------------------------------------------
inline void HarrisFilterVert7Neon(int& w, int srcStride,
    int32_t*& dst0, int32_t*& dst1, int32_t*& dst2, 
    int16_t*& src0, int16_t*& src1, int16_t*& src2, 
    int16_t*& src3, int16_t*& src4, int16_t*& src5, int16_t*& src6)
{
    while (w > 4)
    {
        // this implementation has more ops (since the shifts are not
        // combined), but is faster due to the minimal use of registers
        // and minimal loads

        int32x4_t acc0,acc1,acc2;
        // process first sample and assign to accumulator registers
        acc0 = vmovl_s16(vld1_s16(src0+(0*srcStride)));
        acc1 = vmovl_s16(vld1_s16(src0+(1*srcStride)));
        acc2 = vmovl_s16(vld1_s16(src0+(2*srcStride)));
        src0 += 4;

        int32x4_t t0,t1,t2;
#define ACC_LOAD(_i_) \
        t0 = vmovl_s16(vld1_s16(src##_i_##+(0*srcStride))); \
        t1 = vmovl_s16(vld1_s16(src##_i_##+(1*srcStride))); \
        t2 = vmovl_s16(vld1_s16(src##_i_##+(2*srcStride))); \
        src##_i_## += 4;
#define ACC_SHIFT(_sh_) \
        t0 = vshlq_n_s32(t0,_sh_); \
        t1 = vshlq_n_s32(t1,_sh_); \
        t2 = vshlq_n_s32(t2,_sh_);
#define ACC_ACCUM() \
        acc0 = vaddq_s32(acc0,t0); \
        acc1 = vaddq_s32(acc1,t1); \
        acc2 = vaddq_s32(acc2,t2);

        ACC_LOAD(1); ACC_SHIFT(2); ACC_ACCUM();
        ACC_LOAD(2); ACC_SHIFT(3); ACC_ACCUM();
        ACC_LOAD(3); ACC_SHIFT(1); ACC_ACCUM();
                        ACC_SHIFT(2); ACC_ACCUM(); // ((<<1)<<2) = <<3
        ACC_LOAD(4); ACC_SHIFT(3); ACC_ACCUM();
        ACC_LOAD(5); ACC_SHIFT(2); ACC_ACCUM();
        ACC_LOAD(6);               ACC_ACCUM();
#undef ACC_LOAD
#undef ACC_SHIFT
#undef ACC_ACCUM
        vst1q_s32(dst0, acc0); dst0 += 4;
        vst1q_s32(dst1, acc1); dst1 += 4;
        vst1q_s32(dst2, acc2); dst2 += 4;
        w-=4;
    }
}


inline void HarrisFilterHoriz7Neon(int& w, float*& dst, 
    int32_t*& buf0, int32_t*& buf1, int32_t*& buf2,
    float  fscaleSq)
{
// non-Neon version - was faster on Mango with hybrid SDK, but not for Apollo
//
// result is (A)+(B<<2)+(C<<3)+(D<<3)+(D<<1)+(E<<3)+(F<<2)+(G)  13 ops
// = A+((((B+((C+D+E)<<1)+F)<<1)+D)<<1)+G   10 ops
// = ((B+((C+D+E)<<1)+F)<<2) + (A+(D<<1)+G) same 10 ops with fewer intermediate registers
#define FILTER4_HORIZ_7(_src_,_res_) \
    { \
        int32_t v0,v1,v2,v3,v4,v5,v6,v7,v8,v9; \
        v0 = *(_src_+0); v1 = *(_src_+1); v2 = *(_src_+2); v3 = *(_src_+3); v4 = *(_src_+4); \
        v5 = *(_src_+5); v6 = *(_src_+6); v7 = *(_src_+7); v8 = *(_src_+8); v9 = *(_src_+9); \
        _res_[0] = (v0)+(v3<<1)+((v1+((v2+v3+v4)<<1)+v5)<<2)+(v6); \
        _res_[1] = (v1)+(v4<<1)+((v2+((v3+v4+v5)<<1)+v6)<<2)+(v7); \
        _res_[2] = (v2)+(v5<<1)+((v3+((v4+v5+v6)<<1)+v7)<<2)+(v8); \
        _res_[3] = (v3)+(v6<<1)+((v4+((v5+v6+v7)<<1)+v8)<<2)+(v9); \
    }

// Neon version
//
// |A B C D|E F G H|I J - - |
// notation: A# is vector |A B C D|, B# is vector |B C D E|, etc.
// ((B#+((C#+D#+E#)<<1)+F#)<<2) + (A#+(D#<<1)+G#)
#define FILTER4_HORIZ_7_NEON(_src_,_res_) \
    { \
        uint32x4_t resr4 = vld1q_u32((const uint32_t*)(_src_+0)); /* A# */ \
        uint32x4_t midl4 = vld1q_u32((const uint32_t*)(_src_+4)); /* E# */ \
        uint32x4_t tmpa4, resl4, tmpb4; \
        resl4 = vextq_u32(resr4,midl4,2); /* C# */ \
        resl4 = vaddq_u32(resl4,midl4);   /* resl = C#+E# */ \
        tmpb4 = vextq_u32(resr4,midl4,3); /* D# */ \
        resl4 = vaddq_u32(resl4,tmpb4);   /* resl = C#+D#+E# */ \
        resl4 = vshlq_n_u32(resl4,1);     /* resl = (C#+D#+E#)<<1 */ \
        tmpa4 = vextq_u32(resr4,midl4,1); /* B# */ \
        resl4 = vaddq_u32(resl4,tmpa4);   /* resl B#+((C#+D#+E#)<<1) */ \
        tmpb4 = vshlq_n_u32(tmpb4,1);     /* D#<<1 */ \
        resr4 = vaddq_u32(resr4,tmpb4);   /* resr = A#+(D#<<1) */ \
        tmpb4 = vld1q_u32((const uint32_t*)(_src_+8)); /* |I J - -| */ \
        tmpa4 = vextq_u32(midl4,tmpb4,1); /* F# */ \
        resl4 = vaddq_u32(resl4,tmpa4);   /* resl = B#+((C#+D#+E#)<<1)+F# */ \
        resl4 = vshlq_n_u32(resl4,2);     /* resl = (B#+((C#+D#+E#)<<1)+F#)<<2 */ \
        tmpa4 = vextq_u32(midl4,tmpb4,2); /* G# */ \
        resr4 = vaddq_u32(resr4,tmpa4);   /* resr = A#+(D#<<1)+G# */ \
        resl4 = vaddq_u32(resl4,resr4);   /* result */ \
        vst1q_u32((uint32_t*)(_res_),resl4); \
    }

    const float32x4_t fScaleSq4 = vmovq_n_f32(fscaleSq);
    while (w > 4)
    {
        int32_t xx[4],yy[4],xy[4];
        FILTER4_HORIZ_7_NEON(buf0,xx); buf0 += 4;
        FILTER4_HORIZ_7_NEON(buf1,yy); buf1 += 4;
        FILTER4_HORIZ_7_NEON(buf2,xy); buf2 += 4;
        int32x4_t ixx4 = vld1q_s32(xx);
        int32x4_t iyy4 = vld1q_s32(yy);
        int32x4_t ixy4 = vld1q_s32(xy);
        float32x4_t fxx4 = vcvtq_f32_s32(ixx4);
        float32x4_t fyy4 = vcvtq_f32_s32(iyy4);
        float32x4_t fxy4 = vcvtq_f32_s32(ixy4);

        float32x4_t fHarris4 = vmulq_f32(fxx4,fyy4); 
        fHarris4 = vmlsq_f32(fHarris4,fxy4,fxy4);
        fHarris4 = vmulq_f32(fHarris4,fScaleSq4);
        vst1q_f32(dst, fHarris4); dst+=4;
        w -= 4;
    }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// end