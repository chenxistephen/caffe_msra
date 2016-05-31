#ifndef _READ_CONFIGURE_H_
#define _READ_CONFIGURE_H_
#include <fstream>
#include <stdio.h>
#include <vector>
#include <tchar.h>
#include <windows.h>

struct ReadConfigureInfoStruct
{
	std::string strVariable;
	int iFlag;
	ReadConfigureInfoStruct()
	{
		iFlag = 0;
	}
};

#define MaxCommentLength 256

#define BeginReadConfigure(File) \
{\
	std::ifstream __fin(File);\
	if (__fin.fail())\
	{\
		printf("Can't open configure file: %s\n", File);\
		throw ("Can't open configure file!\n");\
	}\
	std::string __str;\
	std::string __strEqu;\
	char __pBuff[MaxCommentLength];\
	std::vector<ReadConfigureInfoStruct> __vInfo;\
	int __iTotalVariableNumber = 0;\
	int __iCurrentVariable;\
	bool __bFlag;\
	__fin >> __str;\
	while(!__fin.eof() &&!__fin.fail())\
	{\
		__iCurrentVariable = 0;\
		__bFlag = false;


#define BeginReadConfigureDefault(DefaultFileName) \
	const char *__pDefaultFileName = NULL;\
	if (argc == 1)\
	{\
		__pDefaultFileName = DefaultFileName;\
	}\
	else\
	{\
		__pDefaultFileName = argv[1];\
		for (int __i = 2; __i < argc; __i++)\
		{\
			::ShellExecuteA(NULL, "open", argv[0], argv[__i], NULL, SW_SHOWNORMAL);\
		}\
	}\
	BeginReadConfigure(__pDefaultFileName)

#define ReadVariable(A) \
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = #A;\
		}\
		if (strcmp(__str.c_str(), #A) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", #A);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			__fin >> A;\
			if (__fin.fail()) \
				printf("Error: Read \"%s\" failed!\n", #A); \
			std::cout << #A << " = " << A << ";\n";\
			__fin.getline(__pBuff, MaxCommentLength);\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


#define ReadCharString(A) \
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = #A;\
		}\
		if (strcmp(__str.c_str(), #A) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", #A);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			std::string __tmpStr;\
			__fin >> __tmpStr;\
			##A = new char[__tmpStr.size()+1];\
			strcpy_s(A, __tmpStr.size()+1, __tmpStr.c_str());\
			std::cout << #A << " = " << A << ";\n";\
			__fin.getline(__pBuff, MaxCommentLength);\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


#define ReadWcharString(A) \
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = #A;\
		}\
		if (strcmp(__str.c_str(), #A) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", #A);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			std::string __tmpStr;\
			__fin >> __tmpStr;\
			##A = new WCHAR[__tmpStr.size()+1];\
			MultiByteToWideChar(0, 0, __tmpStr.c_str(), (int)__tmpStr.size()+1, ##A, (int)__tmpStr.size()+1);\
			std::wcout << L#A << L" = " << A << L";\n";\
			__fin.getline(__pBuff, MaxCommentLength);\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


#define ReadArray(A, C) \
		sprintf_s(__pBuff, MaxCommentLength, "%s[%d]", #A, C);\
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = __pBuff;\
		}\
		if (strcmp(__str.c_str(), __pBuff) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", __pBuff);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			for (int __i = 0; __i < C; __i++)\
			{\
				__fin >> ##A[__i];\
				std::cout << #A << "[" << __i << "] = " << ##A[__i] <<";\n";\
				__fin.getline(__pBuff, MaxCommentLength);\
			}\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


#define ReadVector(A) \
		sprintf_s(__pBuff, MaxCommentLength, "%s[", #A);\
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = __pBuff;\
		}\
		if (strncmp(__str.c_str(), __pBuff, strlen(__pBuff)) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", __pBuff);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			int __C = atoi(__str.c_str()+strlen(__pBuff));\
			__fin >> __strEqu;\
			##A.resize(__C);\
			for (int __i = 0; __i < __C; __i++)\
			{\
				__fin >> ##A[__i];\
				std::cout << #A << "[" << __i << "] = " << ##A[__i] <<";\n";\
				__fin.getline(__pBuff, MaxCommentLength);\
			}\
			__bFlag = true;\
		}\
		__iCurrentVariable++;

// char **str;
// ...
// ReadCharStringArray(str, 10);
// ...
// DestroyCharStringArray(str, 10);
#define ReadCharStringArray(A, C) \
		sprintf_s(__pBuff, MaxCommentLength, "%s[%d]", #A, C);\
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = __pBuff;\
		}\
		if (strcmp(__str.c_str(), __pBuff) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", __pBuff);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			##A = new char*[C];\
			if (##A == NULL)\
			{\
				printf("Out of memory!\n");\
				throw ("Out of memory!\n");\
			}\
			for (int __i = 0; __i < (int)C; __i++)\
			{\
				std::string __strTmp;\
				__fin >> __strTmp;\
				##A[__i] = new char[__strTmp.size()+1];\
				strcpy_s(##A[__i], __strTmp.size()+1, __strTmp.c_str());\
				std::cout << #A << "[" << __i << "] = " << ##A[__i] <<";\n";\
				__fin.getline(__pBuff, MaxCommentLength);\
			}\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


// char **str;
// ...
// ReadCharStringArray(str, 10);
// ...
// DestroyCharStringArray(str, 10);
#define ReadWcharStringArray(A, C) \
		sprintf_s(__pBuff, MaxCommentLength, "%s[%d]", #A, C);\
		if ((int)__vInfo.size() > __iCurrentVariable && __vInfo[__iCurrentVariable].iFlag == 0)\
		{\
			__vInfo[__iCurrentVariable].iFlag = 1;\
			__vInfo[__iCurrentVariable].strVariable = __pBuff;\
		}\
		if (strcmp(__str.c_str(), __pBuff) == 0)\
		{\
			if ((int)__vInfo.size() <= __iCurrentVariable)\
			{\
				__vInfo.resize(__iCurrentVariable+1);\
				__vInfo[__iCurrentVariable].strVariable = #A;\
			}\
			if (__vInfo[__iCurrentVariable].iFlag == 2)\
				printf("Warning: \"%s\" is reassigned!\n", __pBuff);\
			__vInfo[__iCurrentVariable].iFlag = 2;\
			__fin >> __strEqu;\
			##A = new WCHAR*[C];\
			if (##A == NULL)\
			{\
				printf("Out of memory!\n");\
				throw ("Out of memory!\n");\
			}\
			for (int __i = 0; __i < (int)C; __i++)\
			{\
				std::string __strTmp;\
				__fin >> __strTmp;\
				##A[__i] = new WCHAR[__strTmp.size()+1];\
				MultiByteToWideChar(0, 0, __strTmp.c_str(), (int)__strTmp.size()+1, ##A[__i], (int)__strTmp.size()+1);\
				std::wcout << L#A << L"[" << __i << L"] = " << ##A[__i] << L";\n";\
				__fin.getline(__pBuff, MaxCommentLength);\
			}\
			__bFlag = true;\
		}\
		__iCurrentVariable++;


#define EndReadConfigure() \
		if (__iTotalVariableNumber != __iCurrentVariable )\
		{\
			__iTotalVariableNumber = __iCurrentVariable;\
			__vInfo.resize(__iTotalVariableNumber);\
		}\
		if (!__bFlag)\
		{\
			if (__str[0] != '#')\
				printf("Warning: \"%s\" is never used!\n", __str.c_str());\
			__fin.getline(__pBuff, MaxCommentLength);\
		}\
		__fin >> __str;\
	}\
	__fin.close();\
	for (int __i = 0; __i < (int)__vInfo.size(); __i++)\
	{\
		if (__vInfo[__i].iFlag == 1)\
		{\
			/*printf("Warning: \"%s\" isn't assigned!\n", __vInfo[__i].strVariable.c_str());*/\
		}\
	}\
}


#define DestroyCharString(A) delete[] (A); ##A = NULL;
#define DestroyCharStringArray(A, C) \
	for (int __i = 0; __i < (int)C; __i++)\
	{\
		delete[] (##A[__i]);\
		##A[__i] = NULL;\
	}\
	delete[] (A);\
	##A = NULL;

#define DestroyWcharStringArray DestroyCharStringArray
#define DestroyWcharString DestroyCharString

#endif