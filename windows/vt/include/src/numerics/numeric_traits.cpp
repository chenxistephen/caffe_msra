
#include "stdafx.h"

#include "vt_numeric_traits.h"

const float vt::numeric_traits<float>::zero = 0.0f;
const float vt::numeric_traits<float>::one = 1.0f;
const float vt::numeric_traits<float>::max = 3.402823466e+38f;
const float vt::numeric_traits<float>::min = 1.175494351e-38f;
const float vt::numeric_traits<float>::epsilon = 1.192092896e-07f;

const double vt::numeric_traits<double>::zero = 0.0;
const double vt::numeric_traits<double>::one = 1.0;
const double vt::numeric_traits<double>::max = 1.7976931348623158e+308;
const double vt::numeric_traits<double>::min = 2.2250738585072014e-308;
const double vt::numeric_traits<double>::epsilon = 2.2204460492503131e-016;

