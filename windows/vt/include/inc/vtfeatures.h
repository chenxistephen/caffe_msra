//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      VisionTools 'feature' library include file
//
//  History:
//      2011/10/14-v-mitoel
//          Created
//
//------------------------------------------------------------------------

#pragma once

#ifndef VT_GCC
#include "../src/features/vt_stabilize.h"
#include "../src/features/vt_harrisdetect.h"
#else
#include "vtcore.h"
#endif
#include "../src/features/rollingshutter.h"
