#pragma once

#include "mex.h"
#include <vtcore.h>

namespace vt
{
	void VtHRToMatlabError(HRESULT hr);

	string VtMatlabToVtString(const mxArray* str);

	mxClassID VtImgToMatlabElType(int type);

	int VtMatlabToImgElType(mxClassID type);

	template <typename T>
	bool VtGetMatlabData(T& val, const mxArray* arr)
	{
		if (arr == nullptr)
			return false;
		val = static_cast<T>(mxGetScalar(arr));
		return true;
	}

	// specialization to prevent warning C4800 (narrowing to bool)
	template <>
	bool VtGetMatlabData(bool& val, const mxArray* arr);

	template <unsigned int N>
	bool VtGetMatlabData(char (&str)[N], const mxArray* arr)
	{
		if (arr == nullptr)
			return false;

		return mxGetString(arr, str, N) == 0;
	}

	template <typename T>
	struct MatlabClassID;

#define MAKE_MATLAB_CLASS_ID_MAP(cppType, mTypeName) \
	template <> \
	struct MatlabClassID<cppType> \
	{ \
		static const mxClassID value = mTypeName; \
	};

	MAKE_MATLAB_CLASS_ID_MAP(bool, mxLOGICAL_CLASS);
	MAKE_MATLAB_CLASS_ID_MAP(int, mxINT32_CLASS);
	MAKE_MATLAB_CLASS_ID_MAP(unsigned int, mxUINT32_CLASS);
	MAKE_MATLAB_CLASS_ID_MAP(float, mxSINGLE_CLASS);

#undef MAKE_MATLAB_CLASS_ID_MAP

	template <typename T>
	bool VtMakeMatlabData(mxArray*& arr, const T& val)
	{
		arr = mxCreateNumericMatrix(1, 1, MatlabClassID<T>::value, mxREAL);
		if (arr == nullptr)
			return false;
		*reinterpret_cast<T*>(mxGetData(arr)) = val;
		return true;
	}

#define VT_MATLAB_INIT_GET_FIELD(structName) \
	const char *errFormat1 = "Missing field '%s' in " ##structName; \
	const char *errFormat2 = "Error parsing field '%s' in " ##structName; \
	mxArray *tmp; char errMsg[1024];

#define VT_MATLAB_GET_FIELD(structname, field, arrname) \
	tmp=mxGetField(arrname,0,#field); \
	if(!tmp) { sprintf_s(errMsg,errFormat1,#field); mexErrMsgTxt(errMsg); } \
	if (!VtGetMatlabData(structname.field, tmp)) \
	{ sprintf_s(errMsg,errFormat2,#field); mexErrMsgTxt(errMsg); }

#define VT_MATLAB_INIT_SET_FIELD(structName) \
	const char *errFormat = "Error setting field '%s' in " ##structName; \
	mxArray *tmp; char errMsg[1024];

#define VT_MATLAB_SET_FIELD(arrname, field, structname) \
	tmp = nullptr; \
	if(!VtMakeMatlabData(tmp, structname.field) || \
	mxGetFieldNumber(arrname, #field) == -1) { \
	sprintf_s(errMsg,errFormat,#field); mexErrMsgTxt(errMsg); } \
	mxSetField(arrname, 0, #field, tmp);

	template <typename T>
	void VtWrapMatlabArrayInImg(CPlanarTypedImg<T>& dst, const mxArray* src)
	{
		const mwSize* dims = mxGetDimensions(src);
		const mwSize ndim = mxGetNumberOfDimensions(src);
		const mwSize sz[3] = { dims[0], ndim > 1 ? dims[1] : 1, 
			ndim > 2 ? dims[2] : 1 };

		if (ElTraits<T>::ElFormat() != VtMatlabToImgElType(mxGetClassID(src)))
			mexErrMsgTxt("invalid image type");

		VtHRToMatlabError(dst.Create((T*)mxGetData(src), (int)sz[0], 
			(int)sz[1], (int)sz[2], (int)sz[0] * sizeof(T)));
	}

	template <typename T>
	mxArray* VtCreateMatlabArrayAsImg(CPlanarTypedImg<T>& dst, int width, 
		int height, int numChannels)
	{
		const mwSize chsz[3] = { width, height, numChannels };
		mxArray* ch = mxCreateNumericArray(3, chsz, 
			VtImgToMatlabElType(ElTraits<T>::ElFormat()), mxREAL);

		if (ch == nullptr)
			mexErrMsgTxt("cannot allocate image");

		VtWrapMatlabArrayInImg(dst, ch);

		return ch;
	}
}