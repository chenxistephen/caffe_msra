//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      Main VisionTools include file
//
//  History:
//      2004/11/08-mattu/swinder
//          Created
//
//------------------------------------------------------------------------

#pragma once

#include "vtcommon.h"

#include "vtnumerics.h"

#include "../src/core/vt_atltypes.h" // for MFC CPoint CSize and CRect
#include "../src/core/vt_pad.h"
#include "../src/core/vt_layer.h"
#include "../src/core/vt_vimage.h"
#include "../src/core/vt_convert.h"
#include "../src/core/vt_convert.inl"
#include "../src/core/vt_bicubicbspline.h"
#include "../src/core/vt_imgmath.h"
#include "../src/core/vt_imgmath.inl"
#include "../src/core/vt_kernel.h"
#include "../src/core/vt_addressgen.h"
#include "../src/core/vt_pixeliterators.h"
#include "../src/core/vt_resample.h"
#include "../src/core/vt_warp.h"
#include "../src/core/vt_stats.h"
#include "../src/core/vt_separablefilter.h"
#include "../src/core/vt_steerablefilter.h"
#include "../src/core/vt_rotate.h"
#include "../src/core/vt_filter.h"
#include "../src/core/vt_morphology.h"

#include "../src/core/vt_pyramid.h"
#if !defined(VT_GCC)
#include "../src/core/vt_timer.h"
#endif

#if !defined(VT_GCC)
#include "../src/core/vt_draw.h"
#include "../src/core/vt_params.h"
#include "../src/core/vt_dib.h"
#include "../src/core/vt_archivereaderwriter.h"
#include "../src/core/vt_blobstore.h"
#include "../src/core/vt_fileblobstore.h"
#include "../src/core/vt_winrtblobstore.h"
#include "../src/core/vt_taskmanager.h"
#include "../src/core/vt_readerwriter.h"
#include "../src/core/vt_transform.h"
#include "../src/core/vt_imagecache.h"
#include "../src/core/vt_utils.h"
#include "../src/core/vt_fft.h"
#include "../src/core/vt_imgdist.h"
#include "../src/core/vt_memstream.h"
#include "../src/core/vt_exifmap.h"
#include "../src/core/vt_image.h"
#include "../src/core/vt_imgdbg.h"
#include "../src/core/vt_video.h"
#endif
