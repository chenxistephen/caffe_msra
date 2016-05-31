#ifndef FAST_OP_BASE_H
#define FAST_OP_BASE_H

#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#  ifndef NOMINMAX
#     define NOMINMAX
#  endif

#  pragma inline_recursion(on)
#  pragma inline_depth(255)

#  ifndef FORCE_INLINE
#     define FORCE_INLINE __forceinline
#  endif

#  ifndef NO_INLINE
#     define NO_INLINE __declspec(noinline)
#  endif

#ifndef HAS_LONG_LONG
#  define HAS_LONG_LONG
#endif

#  pragma warning(disable : 4996)
#  pragma warning(disable : 4714)

#else
#  define FORCE_INLINE inline
#endif

#ifndef MAX
#  define MAX(arg1,arg2) (((arg1)>(arg2))? (arg1) : (arg2))
#endif

#ifndef MIN
#  define MIN(arg1,arg2) (((arg1)<(arg2))? (arg1) : (arg2))
#endif

#ifndef ASSERT 
#  ifdef _DEBUG
#     ifdef _WIN64
#        define ASSERT( COND )	if(bool cond = !(COND)) __debugbreak();
#     else
#        define ASSERT( COND )	if(bool cond = !(COND)) __asm {int 3};
#     endif
#  else
#     define ASSERT( COND )
#  endif
#endif

#define ASSERT_ASSUME( COND )  {ASSERT(COND); __assume(COND);}

#ifdef _DEBUG
#  define DEBUG_ONLY(code) code
#  ifdef USE_MULTI_THREAD
#     undef USE_MULTI_THREAD
#  endif
#  define USE_MULTI_THREAD
#else
#  define DEBUG_ONLY(code)
#  define USE_MULTI_THREAD
#endif

#define DBL_DIG         15                      /* # of decimal digits of precision */
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#define DBL_MANT_DIG    53                      /* # of bits in mantissa */
#define DBL_MAX         1.7976931348623158e+308 /* max value */
#define DBL_MAX_10_EXP  308                     /* max decimal exponent */
#define DBL_MAX_EXP     1024                    /* max binary exponent */
#define DBL_MIN         2.2250738585072014e-308 /* min positive value */
#define DBL_MIN_10_EXP  (-307)                  /* min decimal exponent */
#define DBL_MIN_EXP     (-1021)                 /* min binary exponent */
#define _DBL_RADIX      2                       /* exponent radix */
#define _DBL_ROUNDS     1                       /* addition rounding: near */

#define FLT_DIG         6                       /* # of decimal digits of precision */
#define FLT_EPSILON     1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#define FLT_GUARD       0
#define FLT_MANT_DIG    24                      /* # of bits in mantissa */
#define FLT_MAX         3.402823466e+38F        /* max value */
#define FLT_MAX_10_EXP  38                      /* max decimal exponent */
#define FLT_MAX_EXP     128                     /* max binary exponent */
#define FLT_MIN         1.175494351e-38F        /* min positive value */
#define FLT_MIN_10_EXP  (-37)                   /* min decimal exponent */
#define FLT_MIN_EXP     (-125)                  /* min binary exponent */
#define FLT_NORMALIZE   0
#define FLT_RADIX       2                       /* exponent radix */
#define FLT_ROUNDS      1                       /* addition rounding: near */

#define LDBL_DIG        DBL_DIG                 /* # of decimal digits of precision */
#define LDBL_EPSILON    DBL_EPSILON             /* smallest such that 1.0+LDBL_EPSILON != 1.0 */
#define LDBL_MANT_DIG   DBL_MANT_DIG            /* # of bits in mantissa */
#define LDBL_MAX        DBL_MAX                 /* max value */
#define LDBL_MAX_10_EXP DBL_MAX_10_EXP          /* max decimal exponent */
#define LDBL_MAX_EXP    DBL_MAX_EXP             /* max binary exponent */
#define LDBL_MIN        DBL_MIN                 /* min positive value */
#define LDBL_MIN_10_EXP DBL_MIN_10_EXP          /* min decimal exponent */
#define LDBL_MIN_EXP    DBL_MIN_EXP             /* min binary exponent */
#define _LDBL_RADIX     DBL_RADIX               /* exponent radix */
#define _LDBL_ROUNDS    DBL_ROUNDS              /* addition rounding: near */

typedef unsigned int    u_int;
typedef unsigned long   u_long;
typedef unsigned char   u_char;
typedef unsigned short  u_short;

typedef unsigned short     u_int16; 
typedef unsigned int       u_int32; 
typedef unsigned long long u_int64; 

#ifdef HAS_LONG_LONG
   typedef long long t_longlong;
   typedef unsigned long long t_ulonglong;
#else
   typedef long t_longlong;
   typedef unsigned long t_ulonglong;
#endif

#define MACRO_CONCATENATE_(X,Y) X##Y
#define MACRO_CONCATENATE(X,Y) MACRO_CONCATENATE_(X,Y)

struct S_Nought{};
extern const S_Nought nought;

namespace ns_base
{
   template <bool x> struct STATIC_ASSERT_FAILED;
   template <> struct STATIC_ASSERT_FAILED<true> {};
}

#define STATIC_ASSERT(B)\
   enum { MACRO_CONCATENATE(SAWeirdName_,__LINE__) = sizeof(ns_base::STATIC_ASSERT_FAILED<bool(B)>) }

namespace ns_base
{
   template <typename P>
   void checked_delete(P *ptr) 
   {
      STATIC_ASSERT(sizeof(*ptr)); 
      delete ptr;
   }
   
   template <typename P>
   void checked_array_delete(P *ptr) 
   {
      STATIC_ASSERT(sizeof(*ptr)); 
      delete[] ptr;
   }
}

#include <limits>
#include <type_traits>

namespace ns_base
{
   void MLTKInit();
}

#endif
