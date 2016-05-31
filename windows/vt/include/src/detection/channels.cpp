#include "stdafx.h"

using namespace vt;

namespace vt_internal 
{	
	HRESULT FillColorChannels(CLumaByteImg& channel, 
		const CLumaByteImg& input, int shrinkFactor)
	{
		VT_HR_BEGIN();

		if (!input.IsValid())
			return E_INVALIDSRC;
		if (!channel.IsValid())
			return E_INVALIDDST;
		if (channel.Width() * shrinkFactor != input.Width() ||
			channel.Height() * shrinkFactor != input.Height())
			return E_INVALIDARG;
		if (shrinkFactor != 4)
			return E_INVALIDARG;

		VT_HR_EXIT(VtSeparableFilterBoxDecimate4to1(channel, 
			channel.Rect(), input, CPoint(0, 0)));

		VT_HR_END();
	}		
}