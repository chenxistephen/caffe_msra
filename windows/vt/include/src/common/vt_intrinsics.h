#pragma once

#include "vt_basetypes.h"

// If you enable this warning, you'll get a "compile with /arch:AVX" warning.
// According to this article 
// http://stackoverflow.com/questions/7839925/using-avx-cpu-instructions
// AVX instructions are always generated where AVX intrinsics are used
// (I've double checked this by looking at the generated code of the 
// ConvertOp class)
// Disabling the warning tells the VC compiler that we know what we're doing 
// and not mix SSE and AVX instructions "improperly", which would result in 
// penalty cycles (that's why the warning).
#pragma warning (disable: 4752)

#ifdef VT_GCC
#if (defined(_M_IX86) || defined(_M_AMD64))
#include <x86intrin.h>
#endif
#define VT_INTERLOCKED_INCREMENT(a) (__sync_add_and_fetch(a, 1))
#define VT_INTERLOCKED_DECREMENT(a) (__sync_add_and_fetch(a, -1))
#else
#if (defined(_M_IX86) || defined(_M_AMD64))
#pragma warning (disable : 6540)
#include <intrin.h>
#pragma warning (default : 6540)
#elif defined(_M_ARM)
#include <arm_neon.h>
#endif
#define VT_INTERLOCKED_INCREMENT(a) InterlockedIncrement(a)
#define VT_INTERLOCKED_DECREMENT(a) InterlockedDecrement(a)
#endif

#ifdef VT_GCC
#define VT_DECLSPEC_ALIGN(nbytes) __attribute__((__aligned__(nbytes)))
#else
#define VT_DECLSPEC_ALIGN(nbytes) __declspec(align(nbytes))
#endif	

#ifdef VT_GCC
#define nullptr NULL 
#endif

#ifdef VT_GCC
#define VT_FORCEINLINE inline
#else
#define VT_FORCEINLINE __forceinline
#endif	


namespace vt {

#include "SSEonNeon.h"

#if (defined(_M_IX86) || defined(_M_AMD64))

// SSE2 emulation of corresponding intrinsics from SSE4 set
#ifndef VT_GCC
#define SSE2_mm_extract_epi8 (x, imm)	((((imm) & 0x1) == 0) ? _mm_extract_epi16((x), (imm) >> 1) & 0xff : _mm_extract_epi16(_mm_srli_epi16((x), 8), (imm) >> 1))
#else
inline int SSE2_mm_extract_epi8(__m128i x, int imm) { return ((imm & 0x1) == 0) ? _mm_extract_epi16(x, imm >> 1) & 0xff : _mm_extract_epi16(_mm_srli_epi16(x, 8), imm >> 1); }
#endif
#define SSE2_mm_extract_epi16(x, imm)	_mm_extract_epi16((x), (imm));
#define SSE2_mm_extract_epi32(x, imm)	_mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#ifdef VT_GCC
#define SSE2_mm_extract_epf32(x, imm)	_mm_cvtss_f32((__m128)_mm_srli_si128((__m128i)(x), 4 * (imm)))
#else
#define SSE2_mm_extract_epf32(x, imm)	_mm_cvtss_f32(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(x), 4 * (imm))))
#endif
#define SSE2_mm_extract_epi64(x, imm)	_mm_cvtsi128_si64(_mm_srli_si128((x), 8 * (imm)))

template<bool BypassCache>
struct StoreIntrinsicSSE;

template<>
struct StoreIntrinsicSSE<false>
{
	static void StoreAligned(Byte* dst, const __m128i& A)
	{ _mm_store_si128((__m128i*)dst, A); }
	
	static void StoreAligned(UInt16* dst, const __m128i& A)
	{ _mm_store_si128((__m128i*)dst, A); }

	static void StoreAligned(float* dst, const __m128& A)
	{ _mm_store_ps(dst, A); }
};

template<>
struct StoreIntrinsicSSE<true>
{
	static void StoreAligned(Byte* dst, const __m128i& A)
	{ _mm_stream_si128((__m128i*)dst, A); }
	
	static void StoreAligned(UInt16* dst, const __m128i& A)
	{ _mm_stream_si128((__m128i*)dst, A); }

	static void StoreAligned(float* dst, const __m128& A)
	{ _mm_stream_ps(dst, A); }
};

#ifndef VT_GCC
template<bool BypassCache>
struct StoreIntrinsicAVX;

template<>
struct StoreIntrinsicAVX<false>
{
	template<typename T>
	static void StoreAligned(T* dst, const __m256& A)
	{ _mm256_store_ps((float*)dst, A); }	
};

template<>
struct StoreIntrinsicAVX<true>
{
	template<typename T>
	static void StoreAligned(T* dst, const __m256& A)
	{ _mm256_stream_ps((float*)dst, A); }
};
#endif // VT_GCC

inline __m128 Load4AlignedSSE(const float* p)
{ return _mm_load_ps(p); }
inline __m128i Load4AlignedSSE(const __m128i* p)
{ return _mm_load_si128(p); }
inline __m128 Load4UnAlignedSSE(const float* p)
{ return _mm_loadu_ps(p); }
inline __m128i Load4UnAlignedSSE(const __m128i* p)
{ return _mm_loadu_si128(p); }
inline void Store4AlignedSSE(float* p, const __m128& v)
{ _mm_store_ps(p, v); }
inline void Store4AlignedSSE(__m128i* p, const __m128i& v)
{ _mm_store_si128(p, v); }
inline void Store4UnAlignedSSE(float* p, const __m128& v)
{ _mm_storeu_ps(p, v); }
inline void Store4UnAlignedSSE(__m128i* p, const __m128i& v)
{ _mm_storeu_si128(p, v); }

#ifndef VT_GCC
inline __m256 Load8AlignedSSE(const float* p)
{ return _mm256_load_ps(p); }
inline __m256i Load8AlignedSSE(const __m256i* p)
{ return _mm256_load_si256(p); }
inline __m256 Load8UnAlignedSSE(const float* p)
{ return _mm256_loadu_ps(p); }
inline __m256i Load8UnAlignedSSE(const __m256i* p)
{ return _mm256_loadu_si256(p); }
inline void Store8AlignedSSE(float* p, const __m256& v)
{ _mm256_store_ps(p, v); }
inline void Store8AlignedSSE(__m256i* p, const __m256i& v)
{ _mm256_store_si256(p, v); }
inline void Store8UnAlignedSSE(float* p, const __m256& v)
{ _mm256_storeu_ps(p, v); }
inline void Store8UnAlignedSSE(__m256i* p, const __m256i& v)
{ _mm256_storeu_si256(p, v); }
#endif

// macros to select Processor Specialization functions with templated store func
#define SelectPSTFunc( __cond, __ret, __f1, __fname, ... ) \
    if (__cond) { __ret = __fname<__f1##AlignedSSE>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedSSE>(__VA_ARGS__); }
#define SelectPSTFunc2( __cond, __ret, __f1, __f2, __fname, ... ) \
    if (__cond) { __ret = __fname<__f1##AlignedSSE,__f2##AlignedSSE>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedSSE,__f2##UnAlignedSSE>(__VA_ARGS__); }
#define SelectPSTFunc2c( __cond1, __cond2, __ret, __f1, __f2, __fname, ... ) \
    if ((__cond1)&&(__cond2)) { __ret = __fname<__f1##AlignedSSE,__f2##AlignedSSE>(__VA_ARGS__); } \
    else if (__cond1) { __ret = __fname<__f1##AlignedSSE,__f2##UnAlignedSSE>(__VA_ARGS__); } \
    else if (__cond2) { __ret = __fname<__f1##UnAlignedSSE,__f2##AlignedSSE>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedSSE,__f2##UnAlignedSSE>(__VA_ARGS__); }

#elif _M_ARM

inline __m128 Load4AlignedNEON(const float* p)
{ return vld1q_f32_ex(p,128); }
inline __m128 Load4AlignedNEON(const __m128i* p)
{ return vld1q_u8_ex((uint8_t*)p,128); }
inline __m128 Load4UnAlignedNEON(const float* p)
{ return vld1q_f32(p); }
inline __m128 Load4UnAlignedNEON(const __m128i* p)
{ return vld1q_u8((uint8_t*)p); }
inline void Store4AlignedNEON(float* p, const __m128& v)
{ vst1q_f32_ex(p, v, 128); }
inline void Store4AlignedNEON(__m128i* p, const __m128i& v)
{ vst1q_u8_ex((uint8_t*)p, v, 128); }
inline void Store4UnAlignedNEON(float* p, const __m128& v)
{ vst1q_f32(p, v); }
inline void Store4UnAlignedNEON(__m128i* p, const __m128i& v)
{ vst1q_u8((uint8_t*)p, v); }

// macros to select Processor Specialization functions with templated store func
#define SelectPSTFunc( __cond, __ret, __f1, __fname, ... ) \
    if (__cond) { __ret = __fname<__f1##AlignedNEON>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedNEON>(__VA_ARGS__); }
#define SelectPSTFunc2( __cond, __ret, __f1, __f2, __fname, ... ) \
    if (__cond) { __ret = __fname<__f1##AlignedNEON,__f2##AlignedNEON>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedNEON,__f2##UnAlignedNEON>(__VA_ARGS__); }
#define SelectPSTFunc2c( __cond1, __cond2, __ret, __f1, __f2, __fname, ... ) \
    if ((__cond1)&&(__cond2)) { __ret = __fname<__f1##AlignedNEON,__f2##AlignedNEON>(__VA_ARGS__); } \
    else if (__cond1) { __ret = __fname<__f1##AlignedNEON,__f2##UnAlignedNEON>(__VA_ARGS__); } \
    else if (__cond2) { __ret = __fname<__f1##UnAlignedNEON,__f2##AlignedNEON>(__VA_ARGS__); } \
    else { __ret = __fname<__f1##UnAlignedNEON,__f2##UnAlignedNEON>(__VA_ARGS__); }

#endif

};
