//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Traits classes
//
//      If you ever find yourself in a templated function writing T(1),
//      or even T(0), in order to get a "zero" or "one" of type T, you
//      should probably be using a traits template.  There are two main
//      problems with casting to get these.  First, there may be
//      sensible int constructors, e.g. for a vector.  Fine, you say,
//      let's construct from a double: T(0.0).  However, this does not
//      work for float, which will narrow, and does not express the
//      right idea.  Second, the "right" idea is to decide what exactly
//      it is that you want, and to ask for it.  With zero for example,
//      one generally wants an initializer for an accumulator, i.e. the
//      additive identity.  So to add all the elements in a vector, write
//
//         T total = numeric_traits<T>::zero;
//         for(...) total += v[i];
//
//      And to find the floating-point precision for a T, use
//         numeric_traits<T>::epsilon;
//      
//
//  History:
//      2008/3/21-awf
//          Created
//      2011/11/4-kopf
//          Added max, min traits, corrected epsilon values
//
//------------------------------------------------------------------------
#pragma once

namespace vt {

template <class T>
struct numeric_traits {
  static const T zero;    // additive identity
  static const T one;     // multiplicative identity
  static const T max;     // maximal value
  static const T min;     // minimal positive value
  typedef double abs_t;   // return type of the "abs" function
  static const abs_t epsilon; // smallest T such that (unit + epsilon) == unit, where unit is numeric_traits<abs_t>::one
};

// Specializations:
template <>
struct numeric_traits<float> {
  static const float zero;    // additive identity
  static const float one;     // multiplicative identity
  static const float max;     // maximal value
  static const float min;     // minimal positive value
  typedef float abs_t;   // return type of the "abs" function
  static const abs_t epsilon; // smallest T such that (unit + epsilon) == unit, where unit is numeric_traits<abs_t>::one
};

template <>
struct numeric_traits<double> {
  static const double zero;    // additive identity
  static const double one;     // multiplicative identity
  static const double max;     // maximal value
  static const double min;     // minimal positive value
  typedef double abs_t;   // return type of the "abs" function
  static const abs_t epsilon; // smallest T such that (unit + epsilon) == unit, where unit is numeric_traits<abs_t>::one
};

};