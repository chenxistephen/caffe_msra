#pragma once

#include <vtcore.h>
#include "vt_planarimage.h"
#include "vtmatlab.h"
#include "mex.h"

struct MatlabToVTTypeMap
{
	mxClassID matlabType;
	int vtType;
};

const MatlabToVTTypeMap m2vtMaps[] = { 
	{ mxUINT8_CLASS, EL_FORMAT_BYTE },
	{ mxUINT16_CLASS, EL_FORMAT_SHORT },
	{ mxINT32_CLASS, EL_FORMAT_INT },
	{ mxSINGLE_CLASS, EL_FORMAT_FLOAT },
	{ mxDOUBLE_CLASS, EL_FORMAT_DOUBLE }		
};

namespace vt
{
	string VtMatlabToVtString(const mxArray* str)
	{
		if (!mxIsChar(str))
			mexErrMsgTxt("not a string");

		vector<char> buf;
		VtHRToMatlabError(buf.resize(1024));

		if (mxGetString(str, &buf[0], (int)buf.size()) != 0)
			mexErrMsgTxt("mxGetString() failed");

		string res;
		res.format_with_resize("%s", &buf[0]);

		return res;
	}

	void VtHRToMatlabError(HRESULT hr)
	{
		if (SUCCEEDED(hr))
			return;

		wchar_t errbuf[256];
		VtErrorToString(hr, errbuf, _countof(errbuf));
		string errstr;
		VtWideCharToMultiByte(errstr, errbuf);
		mexErrMsgTxt(errstr);
	}

	mxClassID VtImgToMatlabElType(int type)
	{
		for (size_t i = 0; i < _countof(m2vtMaps); ++i)
		{
			if (m2vtMaps[i].vtType == type)
				return m2vtMaps[i].matlabType;
		}

		mexErrMsgTxt("unsupported format");
		return mxUNKNOWN_CLASS;
	}

	int VtMatlabToImgElType(mxClassID type)
	{
		for (size_t i = 0; i < _countof(m2vtMaps); ++i)
		{
			if (m2vtMaps[i].matlabType == type)
				return m2vtMaps[i].vtType;
		}

		mexErrMsgTxt("unsupported format");
		return -1;
	}

	template <>
	bool VtGetMatlabData(bool& val, const mxArray* arr)
	{
		if (arr == nullptr)
			return false;
		val = mxGetScalar(arr) != 0;
		return true;
	}
}