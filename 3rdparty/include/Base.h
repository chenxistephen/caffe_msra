#ifndef BASE_H
#define BASE_H

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
