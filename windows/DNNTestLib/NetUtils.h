
#pragma once

#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <stdlib.h>

#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif

namespace DNNTestLib
{
	enum ConvErrorCode : int
	{
		Success = 0,
		ConfigFileError,
		WeightFileError,
		MeanFileError,
		MemoryAllocationError, // you'd better exit the program when receiving this error
		ParameterError,
		ModelPackageFileError,
		BufferError,
		InitError,
		FATAL,
	};

	struct ConvException
	{
		ConvErrorCode code;
		const char *message;

		ConvException(ConvErrorCode err_code, const char *msg)
			: code(err_code), message(msg)
		{}
	};

	inline bool ThrowException(ConvErrorCode err_code, 
		const char *err_msg, bool assert_condition = false)
	{
		if (!assert_condition)
		{
			std::cout << "Error code " << int(err_code) << ": " << err_msg << std::endl;
			//throw ConvException(err_code, err_msg);
			return true;
		}
		return false;
	}

	template <class T>
	inline void SaveRelease(T *&ptr)
	{
		if (ptr == NULL)
			return;
		delete ptr;
		ptr = NULL;
	}

	template <class T>
	inline void SaveAlignedRelease(T *&arr)
	{
		if (arr == NULL)
			return;
		_aligned_free(arr);
		arr = NULL;
	}

	template <class T>
	inline T ReadData(std::istream &stream)
	{
		T r;
		stream.read((char *)&r, sizeof(T));
		ThrowException(BufferError, "Read data file!", !stream.fail());
		return r;
	}

	template <class T>
	inline void WriteData(std::ostream &stream, const T &t)
	{
		stream.write((const char *)&t, sizeof(T));
	}

	template <class T>
	inline int ReadArray(std::istream &stream, T *arr)
	{
		int n = ReadData<int>(stream);
		stream.read((char *)arr, n * sizeof(T));
		return n;
	}

	template <class T>
	inline void WriteArray(std::ostream &stream, const T *arr, int n)
	{
		WriteData(stream, n);
		stream.write((const char *)arr, n * sizeof(T));
	}

	inline std::string ReadString(std::istream &stream)
	{
		char buf[200];
		int n = ReadArray(stream, buf);
		assert(n < 199);
		buf[n] = '\0';

		return std::string(buf);
	}

	inline void WriteString(std::ostream &stream, const std::string &str)
	{
		WriteArray(stream, str.c_str(), (int)str.size());
	}
}