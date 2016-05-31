#pragma once

#include "vtcore.h"

#include "vt_detector.h"

namespace vt
{
	enum NMSMethod
	{
		NMSMInvalid = -1,
		NMSMExhaustive,
		NMSMGreedy
	};

	enum NMSDenominator
	{
		NMSDInvalid = -1,
		NMSDUnion,
		NMSDMin
	};

	HRESULT VtNonMaximumSuppression(CDetection* inOutDetections, 
		unsigned int numDetectionsIn, unsigned int* numDetectionsOut, 
		int* outDetectionIndex, 
		float overlapThreshold, NMSMethod method = NMSMGreedy, 
		NMSDenominator denominator = NMSDMin);
}
