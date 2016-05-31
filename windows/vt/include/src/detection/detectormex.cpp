#include <vector>
#include <memory>

#include <vtcore.h>

#include "vt_planarimage.h"
#include "vt_detector.h"

#include "vtmatlab.h"
#include "mex.h"

using namespace vt;

std::vector<std::shared_ptr<CDetector>> g_detector;

void PrintUsage(bool isError = false)
{
	const char* msg = "usage: \nvtdetector('init', filename(s))\n"
		"heatmap = vtdetector('run', channels)";

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

void AddDetector(const char* filename)
{
	wstring fname;
	VtMultiByteToWideChar(fname, filename);
	
	g_detector.push_back(std::make_shared<CDetector>());

	CDetectorParams params;
	params.inputIsTransposed = true;
	VtHRToMatlabError(g_detector.back()
		->InitializeFromFile(fname.get_constbuffer(), params));
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
		if (nrhs != 2)
		{
			PrintUsage(true);
			return;
		}

		g_detector.clear();
				
		if (!mxIsCell(prhs[1]))
		{
			AddDetector(VtMatlabToVtString(prhs[1]));
		}
		else
		{
			for (mwSize i = 0; i < mxGetNumberOfElements(prhs[1]); ++i)
				AddDetector(VtMatlabToVtString(mxGetCell(prhs[1], i)));			
		}
	}
	else if (cmd == "run")
	{
		if (nrhs != 2)
		{
			PrintUsage(true);
			return;
		}
		
		if (g_detector.empty())
			mexErrMsgTxt("not initialized");

		CPlanarByteImg channelsTmp;
		VtWrapMatlabArrayInImg(channelsTmp, prhs[1]);

		// make copy so that rows are aligned
		// (VT default is 64 bytes, detector API requires multiple of 16)
		CPlanarByteImg channels;
		VtHRToMatlabError(channels.Create(channelsTmp.Width(), 
			channelsTmp.Height(), channelsTmp.Bands()));
		for (int i = 0; i < channelsTmp.Bands(); ++i)
			VtHRToMatlabError(channelsTmp.Band(i).CopyTo(channels.Band(i)));

		const mwSize dims[2] = { 1, (mwSize)g_detector.size() };
		mxArray* heatmaps = mxCreateCellArray(2, dims);

		for (int i = 0; i < (int)g_detector.size(); ++i)
		{
			CPlanarByteImg output;
			mxSetCell(heatmaps, i, 
				VtCreateMatlabArrayAsImg(output, channels.Width(),
				channels.Height(), 1));

			if (!g_detector[i]->IsInputSizeRegistered(channels.StrideBytes(),
				channels.Height(), channels.Bands()))
			{
#if 0
				mexPrintf("detector: registering new size: %d %d %d\n", 
					channels.StrideBytes(), channels.Height(), 
					channels.Bands());
#endif

				VtHRToMatlabError(g_detector[i]->RegisterInputSize(channels));
			}

			VtHRToMatlabError(
				g_detector[i]->Evaluate(output.Band(0), channels));
		}
		
		plhs[0] = heatmaps;		
	}
	else
	{
		PrintUsage(true);
		return;
	}
}
