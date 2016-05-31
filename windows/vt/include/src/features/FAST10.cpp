//
// Copyright (c) Microsoft Corporation. All rights reserved.
//
#include "stdafx.h"
using namespace vt;
#include "FAST10.h"

#if defined(_M_ARM)

#define vld4_u8(pcD)        ( __neon_D1Adr( 0xf420000f, (pcD)) )
#define vld4q_u8(pcQ)       ( __neon_Q1Adr( 0xf420010f, (pcQ)) )
#define vst4_u8(pD, D4)     ( __neon_AdrD1( 0xf400000f, (pD), (D4)) )
#define vst4q_u8(pD, Q4)    ( __neon_AdrQ1( 0xf400010f, (pD), (Q4)) )

#define neon_const_uint8x8_t(p0,p1,p2,p3,p4,p5,p6,p7) \
    ( ((uint64_t)(p0)<<56)| \
      ((uint64_t)(p1)<<48)| \
      ((uint64_t)(p2)<<40)| \
      ((uint64_t)(p3)<<32)| \
      ((uint64_t)(p4)<<24)| \
      ((uint64_t)(p5)<<16)| \
      ((uint64_t)(p6)<< 8)| \
      ((uint64_t)(p7)<< 0) )

#endif

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// FAST corner detect 10
//
// Implementation is based on "Machine learning for high-speed corner detection" published by
// Edward Rosten and Tom Drummond; Department of Engineering, Cambridge University, UK 2006
//
// function operates on 8 bits/pixel image data; examines 16 pixel bresenham circle ring around
// each pixel; if at least 10 contiguous of the 16 pixels are different from the center by a
// threshold (high or low), then that pixel is marked as a corner
//
// ring pixels are labeled as:
//  -----------------------------
//  |   |   | P | A | B |   |   | 
//  -----------------------------
//  |   | O |   |   |   | C |   | 
//  -----------------------------
//  | N |   |   |   |   |   | D | 
//  -----------------------------
//  | M |   |   | X |   |   | E | 
//  -----------------------------
//  | L |   |   |   |   |   | F | 
//  -----------------------------
//  |   | K |   |   |   | G |   | 
//  -----------------------------
//  |   |   | J | I | H |   |   | 
//  -----------------------------
//
// the goal is to quickly prove that pixel is not a corner for early exit, then confirm exactly that 
// pixels have at least 10 contiguous ring pixels that exceed the same (high or low) threshold
//
// test logic to prove that a pixel is not a FAST10 corner (and, by the process of elimination, that it
// *is* a FAST10 corner) is to compare each pixel to the 5 opposing pixels (such as A against G-K), and
// if neither that pixel nor the opposing 5 pixels are outside the same threshold, then the pixel cannot
// be a FAST10 corner; this can be done incrementally, so if it is shown that, for example, either A or
// I pass the FAST10 candidacy, then A can subsequently tested independently against G,H,J,K; also the
// opposing pixels can be grouped with & and tested against the single pixel (since (a|b)&(a|c) == (a|(b&c)))
//
// conceptually there are 16*5 tests comparing each pixel to the 5 opposite, but since tests are duplicated
// across neighboring pixels there are actually only 40 tests that have to be done (16 * (5/2))
//
// offsets for ring pixels and mask values
// (sadly these cannot be grouped into a structure since the compiler does not efficiently use members as symbolic constants)
#define RING_PIXEL_OFFSET(_Label_,_Xco_,_Yco_,_Xlo_,_Ylo_) \
static const int RP##_Label_##cxo = _Xco_; \
static const int RP##_Label_##cyo = _Yco_; \
static const int RP##_Label_##lxo = _Xlo_; \
static const int RP##_Label_##lyo = _Ylo_;
RING_PIXEL_OFFSET(A, 0,-3, 1, 0)
RING_PIXEL_OFFSET(B, 1,-3, 1, 0)
RING_PIXEL_OFFSET(C, 2,-2, 1, 1)
RING_PIXEL_OFFSET(D, 3,-1, 1, 1)
RING_PIXEL_OFFSET(E, 3, 0, 0, 1)
RING_PIXEL_OFFSET(F, 3, 1, 0, 1)
RING_PIXEL_OFFSET(G, 2, 2,-1, 1)
RING_PIXEL_OFFSET(H, 1, 3,-1, 1)
RING_PIXEL_OFFSET(I, 0, 3,-1, 0)
RING_PIXEL_OFFSET(J,-1, 3,-1, 0)
RING_PIXEL_OFFSET(K,-2, 2,-1,-1)
RING_PIXEL_OFFSET(L,-3, 1,-1,-1)
RING_PIXEL_OFFSET(M,-3, 0, 0,-1)
RING_PIXEL_OFFSET(N,-3,-1, 0,-1)
RING_PIXEL_OFFSET(O,-2,-2, 1,-1)
RING_PIXEL_OFFSET(P,-1,-3, 1,-1)

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// C version
//
inline void FASTCornerDetect10CPixel(int x, int y, FAST10Corners &corners, const unsigned char *src, int stride, int threshold)
{
    const unsigned char *center = src + (y*stride) + x;
    int pix = *center;
    int threshH = pix + threshold;
    int threshL = pix - threshold;
            
// load one of 16 ring pixels and test against high and low threshold values
#define LOAD_AND_TEST_PIX(_label_) \
unsigned int test##_label_ = 0x0; \
{ \
	int _pix_ = *(center+((RP##_label_##cyo)*stride)+(RP##_label_##cxo)); \
	if (_pix_ > threshH) { test##_label_ |= 0x1; } \
	if (_pix_ < threshL) { test##_label_ |= 0x2; } \
}
// test one pixel against one or more opposing positions in ring - if neither one nor all others are outside threshold, then pixel is not a FAST10 corner
#define TEST1_1(l1, l2)       { testAccum &= (test##l1 | (test##l2)); if (!(testAccum)) { return; } }
#define TEST1_2(l1, l2,l3)    { testAccum &= (test##l1 | (test##l2&test##l3)); if (!(testAccum)) { return; } }
#define TEST1_3(l1, l2,l3,l4) { testAccum &= (test##l1 | (test##l2&test##l3&test##l4)); if (!(testAccum)) { return; } }
            
    // load and test 'compass' and diagonal pixels
    unsigned int testAccum = 0x3;
    // incrementally load up to every other pixel in ring and do tests (this is fast and eliminates most non-corners) (12 tests)
    LOAD_AND_TEST_PIX(A) LOAD_AND_TEST_PIX(I) TEST1_1(A,I)
    LOAD_AND_TEST_PIX(E) LOAD_AND_TEST_PIX(M) TEST1_1(E,M)
    LOAD_AND_TEST_PIX(C) TEST1_2(C, I,M) 
    LOAD_AND_TEST_PIX(K) TEST1_3(K, A,C,E)
    LOAD_AND_TEST_PIX(G) TEST1_2(G, M,A)
    LOAD_AND_TEST_PIX(O) TEST1_3(O, E,G,I)
    // complete compass and diagonal tests first, then do remaining tests (28 total)
    LOAD_AND_TEST_PIX(L) LOAD_AND_TEST_PIX(N) TEST1_2(E, L,N) 
    LOAD_AND_TEST_PIX(P) LOAD_AND_TEST_PIX(B) TEST1_2(I, B,P) TEST1_2(G, N,P) TEST1_1(B, L)
    LOAD_AND_TEST_PIX(D) LOAD_AND_TEST_PIX(F) TEST1_2(M, D,F) TEST1_2(K, B,D) TEST1_2(D, L,N) TEST1_3(F, L,N,P)
    LOAD_AND_TEST_PIX(H) LOAD_AND_TEST_PIX(J) TEST1_2(A, H,J) TEST1_2(O, F,H) TEST1_2(C, J,L) TEST1_3(H, N,B,P) TEST1_3(J, P,B,D)
#undef LOAD_AND_TEST_PIX            
#undef TEST1_1  
#undef TEST1_2  
#undef TEST1_3  
	corners.Add((int)(center-src));
}
void FASTCornerDetect10C(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold)
{
	for (int y=3; y<(h-3); y++)
	{
		for (int x=3; x<(w-3); x++)
		{
            FASTCornerDetect10CPixel(x,y,corners,src,stride,threshold);
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(_M_ARM)
//
// single-bit masks for positions within AddCorner mask of each of the 8 pixels (interleaved by the vand/vpadd instructions) - hi/lo masks are ORd together
static const UINT32 positionMasks[8] =
    { (0x01<<0)|(0x80<<0),(0x02<<0)|(0x40<<0),(0x04<<8)|(0x20<<8),(0x08<<8)|(0x10<<8),(0x10<<16)|(0x08<<16),(0x20<<16)|(0x04<<16),(0x40<<24)|(0x02<<24),(0x80<<24)|(0x01<<24) };

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// NEON version
//
// operates on 8 pixels at a time (a single 'y' and 8 'x's); approach is virtually identical to the
// C implementation; 
//
// NEON vector booleans are (at least) 8 bits and are set to 0xFF or 0x00 by conditional instructions; code below
// combines 8 vector booleans together by ANDing out all but one bit each in a different position, then
// using the vpadd instruction to add adjacent vectors to combine them; the high and low masks are also
// in different positions so the high and low results can be ORd together for the final vector test bitmask
//
void FASTCornerDetect10Neon(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold)
{
    // preload first 6 scanlines
    int16x8_t vthreshold = vdupq_n_s16((int16_t)threshold);
	for (int y=3; y<(h-3); y++)
	{
        int x; // scope x outside loop for vector remainder processing
		for (x=3; x<(w-3-7); x=x+8)
		{
            const unsigned char *center = src + (y*stride)+x;
            int16x8_t threshH, threshL;
            {
                int16x8_t pixC = vreinterpretq_i16_u16(vmovl_u8(vld1_u8((const uint8x8_t*)center))); // load uint8 and convert to int16
                threshH = vaddq_s16(pixC,vthreshold);
                threshL = vsubq_s16(pixC,vthreshold);
            }
   
// load pixels and compute threshold booleans
#define LOAD_AND_TEST_PIX(_label_) \
register UINT32 test##_label_; \
{ \
int16x8_t pixval = vreinterpretq_i16_u16(vmovl_u8(vld1_u8((const uint8x8_t*)(center+((RP##_label_##cyo)*stride)+(RP##_label_##cxo))))); \
uint8x8_t tmp1 = vmovn_u16(vcgtq_s16(pixval,threshH)); tmp1 = vand_u8(tmp1,boolMergeMaskH); tmp1 = vpadd_u8(tmp1,tmp1); \
uint8x8_t tmp2 = vmovn_u16(vcltq_s16(pixval,threshL)); tmp2 = vand_u8(tmp2,boolMergeMaskL); tmp2 = vpadd_u8(tmp2,tmp2); \
test##_label_ = vget_lane_u32(vorr_u32(tmp1,tmp2),0); \
}
// test one pixel against one or more opposing positions in ring - if neither one nor all others are outside threshold, then pixel is not a FAST10 corner
#define TEST1_1(l1, l2)       { testAccum &= (test##l1 | (test##l2)); if (!(testAccum)) { continue; } }
#define TEST1_2(l1, l2,l3)    { testAccum &= (test##l1 | (test##l2&test##l3)); if (!(testAccum)) { continue; } }
#define TEST1_3(l1, l2,l3,l4) { testAccum &= (test##l1 | (test##l2&test##l3&test##l4)); if (!(testAccum)) { continue; } }
            
            uint8x8_t boolMergeMaskH, boolMergeMaskL; 
            boolMergeMaskH.n64_u64 = neon_const_uint8x8_t( 1<<7, 1<<6, 1<<5, 1<<4, 1<<3, 1<<2, 1<<1, 1<<0 );
            boolMergeMaskL.n64_u64 = neon_const_uint8x8_t( 1<<0, 1<<1, 1<<2, 1<<3, 1<<4, 1<<5, 1<<6, 1<<7 );
            register UINT32 testAccum = 0xffffffff;
            // incrementally load up to every other pixel in ring and do tests (this is fast and eliminates most non-corners) (12 tests)
            LOAD_AND_TEST_PIX(A) LOAD_AND_TEST_PIX(I) TEST1_1(A,I)
            LOAD_AND_TEST_PIX(E) LOAD_AND_TEST_PIX(M) TEST1_1(E,M)
            LOAD_AND_TEST_PIX(C) TEST1_2(C, I,M) 
            LOAD_AND_TEST_PIX(K) TEST1_3(K, A,C,E)
            LOAD_AND_TEST_PIX(G) TEST1_2(G, M,A)
            LOAD_AND_TEST_PIX(O) TEST1_3(O, E,G,I)
            // complete compass and diagonal tests first, then do remaining tests (28 total)
            LOAD_AND_TEST_PIX(L) LOAD_AND_TEST_PIX(N) TEST1_2(E, L,N) 
            LOAD_AND_TEST_PIX(P) LOAD_AND_TEST_PIX(B) TEST1_2(I, B,P) TEST1_2(G, N,P) TEST1_1(B, L)
            LOAD_AND_TEST_PIX(D) LOAD_AND_TEST_PIX(F) TEST1_2(M, D,F) TEST1_2(K, B,D) TEST1_2(D, L,N) TEST1_3(F, L,N,P)
            LOAD_AND_TEST_PIX(H) LOAD_AND_TEST_PIX(J) TEST1_2(A, H,J) TEST1_2(O, F,H) TEST1_2(C, J,L) TEST1_3(H, N,B,P) TEST1_3(J, P,B,D)
#undef LOAD_AND_TEST_PIX
#undef TEST1_1
#undef TEST1_2
#undef TEST1_3

#define TEST_ADD_CORNER(_pos_) \
if (testAccum & positionMasks[_pos_]) { corners.AddFast((int)(center+(_pos_)-src)); }
            corners.CheckAdd8();
            TEST_ADD_CORNER(0) TEST_ADD_CORNER(1) TEST_ADD_CORNER(2) TEST_ADD_CORNER(3)
            TEST_ADD_CORNER(4) TEST_ADD_CORNER(5) TEST_ADD_CORNER(6) TEST_ADD_CORNER(7)
#undef TEST_ADD_CORNER            
        }
        // do remainder (<8 pixels wide) of row with C version
        for (;x<(w-3);x++)
        {
            FASTCornerDetect10CPixel(x,y,corners,src,stride,threshold);
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
#endif

#if (defined(_M_IX86) || defined(_M_AMD64))
//////////////////////////////////////////////////////////////////////////////////////////////////
//
// SSE version
//
// operates on 16 pixels at a time (a single 'y' and 16 'x's); approach is virtually identical to the
// C implementation; 
//
// does 16 at a time in 8bpp registers by using saturated add/subtract
//
void FASTCornerDetect10SSE(FAST10Corners &corners, const unsigned char *src, int w, int h, int stride, int threshold)
{
    __m128i vthresh = _mm_set1_epi8((Byte)threshold);
    __m128i vzero = _mm_setzero_si128();
	for (int y=3; y<(h-3); y++)
	{
        int x; // scope x outside loop for vector remainder processing
		for (x=3; x<(w-3-15); x=x+16)
		{
            const unsigned char *center = src + (y*stride)+x;
            __m128i pixC = _mm_loadu_si128((__m128i*)center);
            __m128i pixCpTH = _mm_adds_epu8(pixC,vthresh);
            __m128i pixCmTH = _mm_subs_epu8(pixC,vthresh);
 
// load pixels and compute threshold booleans
#define LOAD_AND_TEST_PIX(_label_) \
register UInt32 test##_label_; \
{ \
__m128i pixT = _mm_loadu_si128((__m128i*)(center+((RP##_label_##cyo)*stride)+(RP##_label_##cxo))); \
__m128i tstLo = _mm_cmpeq_epi8(vzero, _mm_subs_epu8(pixT,pixCpTH) ); \
__m128i tstHi = _mm_cmpeq_epi8(vzero, _mm_subs_epu8(pixCmTH,pixT) ); \
test##_label_ = ~(_mm_movemask_epi8(tstLo) | ((_mm_movemask_epi8(tstHi))<<16)); \
}
// test one pixel against one or more opposing positions in ring - if neither one nor all others are outside threshold, then pixel is not a FAST10 corner
#define TEST1_1(l1, l2)       { testAccum &= (test##l1 | (test##l2)); if (!(testAccum)) { continue; } }
#define TEST1_2(l1, l2,l3)    { testAccum &= (test##l1 | (test##l2&test##l3)); if (!(testAccum)) { continue; } }
#define TEST1_3(l1, l2,l3,l4) { testAccum &= (test##l1 | (test##l2&test##l3&test##l4)); if (!(testAccum)) { continue; } }
            
            register UInt32 testAccum = 0xffffffff;
            // incrementally load up to every other pixel in ring and do tests (this is fast and eliminates most non-corners) (12 tests)
            LOAD_AND_TEST_PIX(A) LOAD_AND_TEST_PIX(I) TEST1_1(A,I)
            LOAD_AND_TEST_PIX(E) LOAD_AND_TEST_PIX(M) TEST1_1(E,M)
            LOAD_AND_TEST_PIX(C) TEST1_2(C, I,M) 
            LOAD_AND_TEST_PIX(K) TEST1_3(K, A,C,E)
            LOAD_AND_TEST_PIX(G) TEST1_2(G, M,A)
            LOAD_AND_TEST_PIX(O) TEST1_3(O, E,G,I)
            // complete compass and diagonal tests first, then do remaining tests (28 total)
            LOAD_AND_TEST_PIX(L) LOAD_AND_TEST_PIX(N) TEST1_2(E, L,N) 
            LOAD_AND_TEST_PIX(P) LOAD_AND_TEST_PIX(B) TEST1_2(I, B,P) TEST1_2(G, N,P) TEST1_1(B, L)
            LOAD_AND_TEST_PIX(D) LOAD_AND_TEST_PIX(F) TEST1_2(M, D,F) TEST1_2(K, B,D) TEST1_2(D, L,N) TEST1_3(F, L,N,P)
            LOAD_AND_TEST_PIX(H) LOAD_AND_TEST_PIX(J) TEST1_2(A, H,J) TEST1_2(O, F,H) TEST1_2(C, J,L) TEST1_3(H, N,B,P) TEST1_3(J, P,B,D)
#undef LOAD_AND_TEST_PIX
#undef TEST1_1
#undef TEST1_2
#undef TEST1_3

#define TEST_ADD_CORNER(_pos_) \
if (testAccum & ((1<<(_pos_))|((1<<(_pos_))<<16))) { corners.AddFast((int)(center+(_pos_)-src)); }
            corners.CheckAdd16();
            TEST_ADD_CORNER( 0) TEST_ADD_CORNER( 1) TEST_ADD_CORNER( 2) TEST_ADD_CORNER( 3)
            TEST_ADD_CORNER( 4) TEST_ADD_CORNER( 5) TEST_ADD_CORNER( 6) TEST_ADD_CORNER( 7)
            TEST_ADD_CORNER( 8) TEST_ADD_CORNER( 9) TEST_ADD_CORNER(10) TEST_ADD_CORNER(11)
            TEST_ADD_CORNER(12) TEST_ADD_CORNER(13) TEST_ADD_CORNER(14) TEST_ADD_CORNER(15)
#undef TEST_ADD_CORNER            
        }
        // do remainder (<16 pixels wide) of row with C version
        for (;x<(w-3);x++)
        {
            FASTCornerDetect10CPixel(x,y,corners,src,stride,threshold);
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Non-Maximum Suppression
//
// computes score for all corner pixels, then scans pixel list for the 8 immediate neighboring
// pixels, and marks the pixel as a 'max' if it is greater than or equal to all of it's neighbors
//
void FASTNonMaxSuppression(FAST10Corners &corners, const unsigned char *src, int stride, int threshold)
{
	//
	// compute score for every corner pixel
	//

    for (int i=0; i<corners.count; i++)
    {
        const unsigned char *addr = src + corners.offset[i];
        int threshH, threshL;
        {
            int pix = (int)(*addr);
            threshH = pix + threshold;
            threshL = pix - threshold;
        }

		// Center pixel is C and threshold is T
		// For each pixel on the ring with luminance L
        //       If it is brighter than C + T, add (L - (C + T)) to SumBright
        //       If it is darker than C-T, add (C - T - L ) to SumDark
		// Score is the max of (SumBright, SumDark)
        int sumBright = 0; 
		int sumDark = 0;
// load one of 16 ring pixels and accumulate, incrementally updating address - needs to be instanced in alphabetical order
#define LOAD_AND_SUM_PIX(_label_) \
{ \
    addr += ((RP##_label_##lxo)+(stride*(RP##_label_##lyo))); \
	int _pix_ = (int)(*addr); \
	if (_pix_ > threshH) { sumBright += (_pix_ - threshH); } \
	else if (_pix_ < threshL) { sumDark += (threshL - _pix_); } \
}
        addr += 3; // set address to pixel E
        LOAD_AND_SUM_PIX(F) LOAD_AND_SUM_PIX(G) LOAD_AND_SUM_PIX(H) LOAD_AND_SUM_PIX(I)
        LOAD_AND_SUM_PIX(J) LOAD_AND_SUM_PIX(K) LOAD_AND_SUM_PIX(L) LOAD_AND_SUM_PIX(M)
        LOAD_AND_SUM_PIX(N) LOAD_AND_SUM_PIX(O) LOAD_AND_SUM_PIX(P) LOAD_AND_SUM_PIX(A)
        LOAD_AND_SUM_PIX(B) LOAD_AND_SUM_PIX(C) LOAD_AND_SUM_PIX(D) LOAD_AND_SUM_PIX(E)
#undef LOAD_AND_SUM_PIX        
		corners.score[i] = max(sumBright,sumDark);
	}

    //
	// do non-max suppression - suppress any for which the score is not >= its 8 immediate neighbors
	//
    // loops over corner pixels finding neighbors for suppression comparisons; outside loop steps forward in list to 
    // test each 'i' pixel against inside loop 'j' pixels later in list, performing both forward and reverse comparisons
    // (i.e. both 'i to j' (forward) and 'j to i' (reverse)); reverse suppressions are noted by setting an upper bit in
    // the offset for that pixel (masked out when using offsets for indexing and comparison); forward suppressions are noted
    // and comparisons continue since potential reverse suppressions need to be found even if that 'i' pixel is already
    // known to be suppressed
    //
    // jOffs, once computed in the first 'j' iteration, holds the number of list corners to skip ahead when looking for
    // neighbors in the subsequent scanline; at the end of each 'j' iteration, jOffs is set to skip to the (x,y+1) neighbor
    // if it is present else to the next corner in the list after the position (x,y+1) would be in
    //
    int jOffs = 0; // count to skip ahead in list when checking for neighbors in next scanline
    int i = 0; // scoped outside loop so last corner can be tested for reverse suppression after loop
    for (; i<(corners.count-1); i++)
    {
        int thisOffs = corners.offset[i];
        int isReverseSuppressed = (thisOffs>>30);
        bool isForwardSuppressed = false;
        corners.offset[i] &= ~(1<<30); // clean up offset - only really needed for testing...
        thisOffs = corners.offset[i]; 
        // first see next corner in list is the (forward) neighbor on this scanline (x+1,y)
        if (((corners.offset[i+1])&(~(1<<30))) == (thisOffs+1))
        {
            if (corners.score[i] < corners.score[i+1]) { isForwardSuppressed = true; }
            else if (corners.score[i+1] < corners.score[i]) { corners.offset[i+1] |= (1<<30); }
        }
        // then look ahead for neighbors in next row
        int thisOffsMax = thisOffs+stride+1;
        for (int j=i+1+jOffs; j<corners.count; j++)
        {
            if (j==(i+1+jOffs)) { jOffs--; } // decrement start offset after it is applied (first iteration only)
            int testOffs = corners.offset[j]; testOffs &= ~(1<<30);
            if (testOffs > thisOffsMax) { break; } // already beyond (x+1,y+1)
            if (testOffs < (thisOffsMax-2)) { jOffs++; continue; } // not reached (x-1,y+1); start past this pixel for next iteration
            if (testOffs==thisOffsMax)
            {
                // found (x+1,y+1) immediate neighbor
                if (corners.score[i] < corners.score[j]) { goto nextCorner; }
                else if (corners.score[j] < corners.score[i]) { corners.offset[j] |= (1<<30); }
                break; // already reached most remote neighbor
            }
            if (testOffs==(thisOffsMax-1))
            {
                // found (x,y+1) immediate neighbor
                if (corners.score[i] < corners.score[j]) { goto nextCorner; }
                else if (corners.score[j] < corners.score[i]) { corners.offset[j] |= (1<<30); }
                // don't advance jOffs - need to test this pixel next iteration
            }
            else
            {
                // found (x-1,y+1) immediate neighbor
                if (corners.score[i] < corners.score[j]) { isForwardSuppressed = true; } // still need to test later reverse suppressions
                else if (corners.score[j] < corners.score[i]) { corners.offset[j] |= (1<<30); }
                jOffs++; // start past this pixel next iteration
            }
		}
        if (!(isReverseSuppressed||isForwardSuppressed)) { corners.AddMax(thisOffs); }
    nextCorner:;
	}
    if (i<corners.count)
    {
        if (!(corners.offset[i] & (1<<30))) { corners.AddMax(corners.offset[i]); } else { corners.offset[i] &= ~(1<<30); }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
//
HRESULT FAST10Detect(vt::vector<HARRIS_FEATURE_POINT>& pts, 
    const vt::CByteImg& src, const vt::CRect& srcRect,
    int border, float threshold)
{
    HRESULT hr = S_OK;
    int w = src.Width();
    int h = src.Height();
    FAST10Corners corners;
    hr = corners.GetHr();
    if (hr != S_OK) { return hr; }
#if (defined(_M_IX86) || defined(_M_AMD64))
    FASTCornerDetect10SSE(corners, src.Ptr(), w, h, src.StrideBytes(), (int)threshold);
#elif defined(_M_ARM)
    FASTCornerDetect10Neon(corners, src.Ptr(), w, h, src.StrideBytes(), (int)threshold);
#else
    FASTCornerDetect10C(corners, src.Ptr(), w, h, src.StrideBytes(), (int)threshold);
#endif
    hr = corners.GetHr();
    if (hr != S_OK) { return hr; }
    FASTNonMaxSuppression(corners, src.Ptr(), src.StrideBytes(), (int)threshold);
    for (int i=0; i<corners.GetMaxCount(); i++)
    {
        int x,y,score; corners.GetMaxCoordinate(i,x,y,score,src.StrideBytes());
        if ((x<border)||(y<border)||(x>(w-border))||(y>(h-border))) { continue; }
        HARRIS_FEATURE_POINT fp;
        fp.x = (float)x;
        fp.y = (float)y;
        fp.score = (float)score;
        pts.push_back(fp);
    }
    return hr;
}

// end
