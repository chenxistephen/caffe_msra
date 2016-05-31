
#pragma once

#include "vtcommon.h"
#include "vtcore.h"

namespace vt {

enum FusionMethod { FuseAverage, FuseFeather, FusePoisson, FusePoissonTiled,
    FuseLaplacian, FuseLaplacianStack, FuseBilateral };

HRESULT VtCompleteImage(IImageReaderWriter * pPyramid,
    vt::CTaskProgress* pProgress = NULL,
    FusionMethod fusionMethod = FuseAverage,
    int numRepetitions = 30, int numRepetitionsDec = 3,
    float discontinuousPenalty = 0.5f,
    float windowEnergyThresh = 4.f,
    const wchar_t * pwszTimingFile = NULL,
    CImg * pDebugUV = NULL);

};
