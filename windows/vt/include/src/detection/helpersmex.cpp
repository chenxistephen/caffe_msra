#include <vtcore.h>

#include "vt_rgb2luv.h"
#include "vt_planarimage.h"

#include "vtmatlab.h"
#include "mex.h"

using namespace vt;

void PrintUsage(bool isError = false)
{
	const char* msg = "usage: \n"
		"luv = vthelpers('rgb2luv', rgb)";

	if (isError)
	{
		mexErrMsgTxt(msg);
	}
	else
	{
		mexPrintf(msg);
		mexPrintf("\n");
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || !mxIsChar(prhs[0]))
	{
		PrintUsage();
		return;
	}

	auto cmd = VtMatlabToVtString(prhs[0]);
	if (cmd == "rgb2luv")
	{
		if (nrhs != 2)		
		{
			PrintUsage(true);
			return;
		}

		CPlanarByteImg rgbPlanar;
		VtWrapMatlabArrayInImg(rgbPlanar, prhs[1]);
		
		if (rgbPlanar.Bands() != 3)
			mexErrMsgTxt("input image must have 3 channels");

		CRGBByteImg rgb;
		VtHRToMatlabError(VtConvertPlanarToPackedImage(rgb, rgbPlanar));
		
		CRGBByteImg bgr;		
		VtHRToMatlabError(VtRGBColorSwapImage(bgr, rgb));

		CByteImg luvPacked;
		VtHRToMatlabError(VtConvertImageRGBToLUVPiotr(luvPacked, bgr));

		CPlanarByteImg luv;
		plhs[0] = VtCreateMatlabArrayAsImg(luv, luvPacked.Width(), 
			luvPacked.Height(), luvPacked.Bands());

		VtHRToMatlabError(VtConvertPackedToPlanarImage(luv, luvPacked));
	}
	else
	{
		PrintUsage(true);
		return;
	}
}
