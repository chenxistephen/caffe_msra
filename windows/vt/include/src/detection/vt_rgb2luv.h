#pragma once

#include "vtcore.h"

namespace vt
{
	HRESULT VtConvertImageRGBToLUVPiotr(CByteImg& dst, 
		const CRGBByteImg& src);
}