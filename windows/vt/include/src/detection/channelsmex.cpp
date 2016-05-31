#include <vtcore.h>
#include <vtdetection.h>

#include "vtmatlab.h"
#include "mex.h"

using namespace vt;

int g_maxInputWidth = 2048;
int g_maxInputHeight = 2048;

typedef CChannelsDefaultStaticParams StaticParams;
CChannelsEngine<StaticParams> g_eng;
CChannelsParams g_chparams;
CDetectionPyramidParams g_dpyrparams;
CChannelsPyramidParams g_chpyrparams;

void PrintUsage(bool isError = false)
{
	const char* msg = "usage: \n"
		"[chparams, pyrparams] = vtchannels('init', [chparams], [pyrparams])\n"
		"channels = vtchannels('run', img)";

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

void ChannelsParamsFromMatlab(CChannelsParams& p, const mxArray* arr)
{
	VT_MATLAB_INIT_GET_FIELD("chparams")

	VT_MATLAB_GET_FIELD(p, numInputChannels, arr);
	VT_MATLAB_GET_FIELD(p, shrinkFactor, arr);
	VT_MATLAB_GET_FIELD(p, numInputBlurIterations, arr);
	VT_MATLAB_GET_FIELD(p, numInputNormBlurIterations, arr);
	VT_MATLAB_GET_FIELD(p, inputNormRegularization, arr);
	VT_MATLAB_GET_FIELD(p, inputNormScaling, arr);
	VT_MATLAB_GET_FIELD(p, numGradNormBlurIterations, arr);
	VT_MATLAB_GET_FIELD(p, gradNormRegularization, arr);
	VT_MATLAB_GET_FIELD(p, gradNormScaling, arr);
	VT_MATLAB_GET_FIELD(p, numOrientBins, arr);
}

void ChannelsParamsToMatlab(mxArray*& arr, const CChannelsParams& p)
{
	const char* fieldNames[] = 
	{
		"numInputChannels",
		"shrinkFactor",
		"numInputBlurIterations",
		"numInputNormBlurIterations",
		"inputNormRegularization",
		"inputNormScaling",
		"numGradNormBlurIterations",
		"gradNormRegularization",
		"gradNormScaling",
		"numOrientBins"
	};

	arr = mxCreateStructMatrix(1, 1, _countof(fieldNames), fieldNames);

	VT_MATLAB_INIT_SET_FIELD("chparams");

	VT_MATLAB_SET_FIELD(arr, numInputChannels, p);
	VT_MATLAB_SET_FIELD(arr, shrinkFactor, p);
	VT_MATLAB_SET_FIELD(arr, numInputBlurIterations, p);
	VT_MATLAB_SET_FIELD(arr, numInputNormBlurIterations, p);
	VT_MATLAB_SET_FIELD(arr, inputNormRegularization, p);
	VT_MATLAB_SET_FIELD(arr, inputNormScaling, p);
	VT_MATLAB_SET_FIELD(arr, numGradNormBlurIterations, p);
	VT_MATLAB_SET_FIELD(arr, gradNormRegularization, p);
	VT_MATLAB_SET_FIELD(arr, gradNormScaling, p);
	VT_MATLAB_SET_FIELD(arr, numOrientBins, p);
}

void PyramidParamsFromMatlab(CDetectionPyramidParams& dp, 
	CChannelsPyramidParams& cp, const mxArray* arr)
{
	VT_MATLAB_INIT_GET_FIELD("pyrparams");

	VT_MATLAB_GET_FIELD(dp, minLevelWidth, arr);
	VT_MATLAB_GET_FIELD(dp, minLevelHeight, arr);
	VT_MATLAB_GET_FIELD(dp, maxNumLevels, arr);
	VT_MATLAB_GET_FIELD(dp, maxNumOctaves, arr);
	VT_MATLAB_GET_FIELD(dp, numLevelsPerOctave, arr);
	VT_MATLAB_GET_FIELD(cp, approxSubOctaveGradients, arr);
	VT_MATLAB_GET_FIELD(cp, approxSubOctaveGradientLambda, arr);
	VT_MATLAB_GET_FIELD(cp, numPostSmoothIterations, arr);
}

void PyramidParamsToMatlab(mxArray*& arr, const CDetectionPyramidParams& dp, 
	const CChannelsPyramidParams& cp)
{
	const char* fieldNames[] = 
	{
		"minLevelWidth",
		"minLevelHeight",
		"maxNumLevels",
		"maxNumOctaves",		
		"numLevelsPerOctave",
		"approxSubOctaveGradients",
		"approxSubOctaveGradientLambda",
		"numPostSmoothIterations"
	};

	arr = mxCreateStructMatrix(1, 1, _countof(fieldNames), fieldNames);

	VT_MATLAB_INIT_SET_FIELD("pyrparams");

	VT_MATLAB_SET_FIELD(arr, minLevelWidth, dp);
	VT_MATLAB_SET_FIELD(arr, minLevelHeight, dp);
	VT_MATLAB_SET_FIELD(arr, maxNumLevels, dp);
	VT_MATLAB_SET_FIELD(arr, maxNumOctaves, dp);
	VT_MATLAB_SET_FIELD(arr, numLevelsPerOctave, dp);
	VT_MATLAB_SET_FIELD(arr, approxSubOctaveGradients, cp);
	VT_MATLAB_SET_FIELD(arr, approxSubOctaveGradientLambda, cp);
	VT_MATLAB_SET_FIELD(arr, numPostSmoothIterations, cp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || !mxIsChar(prhs[0]))
	{
		PrintUsage();
		return;
	}

	auto cmd = VtMatlabToVtString(prhs[0]);
	if (cmd == "init")
	{
		g_chparams = CChannelsParams();
		g_dpyrparams = CDetectionPyramidParams();
		g_chpyrparams = CChannelsPyramidParams();
		
		if (nrhs > 1)		
			ChannelsParamsFromMatlab(g_chparams, prhs[1]);
		if (nrhs > 2)
			PyramidParamsFromMatlab(g_dpyrparams, g_chpyrparams, prhs[2]);
		
		// make sure "transposed" flag is set since we're in Matlab
		g_chparams.inputIsTransposed = true;
		VtHRToMatlabError(g_eng.Initialize(g_maxInputWidth, g_maxInputHeight, 
			g_chparams));

		ChannelsParamsToMatlab(plhs[0], g_chparams);
		PyramidParamsToMatlab(plhs[1], g_dpyrparams, g_chpyrparams);		
	}
	else if (cmd == "run")
	{
		if (nrhs != 2)
		{
			PrintUsage(true);
			return;
		}

		if (!g_eng.IsInitialized())
			mexErrMsgTxt("not initialized");

		CPlanarByteImg input;
		VtWrapMatlabArrayInImg(input, prhs[1]);

		if (input.Width() > g_maxInputWidth || 
			input.Height() > g_maxInputHeight)
		{
			mexErrMsgTxt("input image too large");
		}

		CChannelsPyramidEngine<StaticParams> peng;
		VtHRToMatlabError(peng.Initialize(input.Width(), input.Height(), 
			g_eng, g_chpyrparams));

		CSize cSize;
		VtHRToMatlabError(g_eng.GetChannelSize(cSize, input.Width(), 
			input.Height()));
		unsigned int numChannels;
		VtHRToMatlabError(g_eng.GetNumChannels(numChannels));

		CDetectionPyramid channels;
		VtHRToMatlabError(channels.Initialize(cSize.cx, cSize.cy, numChannels,
			g_dpyrparams));		
		CByteImg packedInput;
		VtHRToMatlabError(VtConvertPlanarToPackedImage(packedInput, input));
		VtHRToMatlabError(peng.ComputeChannelsPyramid(channels, packedInput));		

		mwSize ndim = 2;
		mwSize dims[2] = { channels.NumLevels(), 1 };

		mxArray* res = mxCreateCellArray(ndim, dims);

		for (unsigned int i = 0; i < channels.NumLevels(); ++i)
		{
			CPlanarByteImg cl;
			VtHRToMatlabError(channels.ShareLevel(cl, i));

			CPlanarByteImg dstTmp;
			mxArray* level = VtCreateMatlabArrayAsImg(dstTmp, 
				cl.Width(), cl.Height(), cl.Bands());
			mxSetCell(res, i, level);

			for (int j = 0; j < cl.Bands(); ++j)
				VtHRToMatlabError(cl.Band(j).CopyTo(dstTmp.Band(j)));
		}

		plhs[0] = res;		
	}
	else
	{
		PrintUsage(true);
		return;
	}
}
