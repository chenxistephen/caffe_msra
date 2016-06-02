#define NOMINMAX
#include <Windows.h>
#include "direct.h"
#include <tchar.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <Pathcch.h>
#include <assert.h>
#include <map>
#include <omp.h>
#include <iostream>
#include <array>
#include <map>
#include "HighResolutionTimer.h"
#include "ImageRecognitionApi.h"
wchar_t* CharToWchar(const char* c)
{
	int len = MultiByteToWideChar(CP_ACP, 0, c, strlen(c), NULL, 0);
	wchar_t* m_wchar = new wchar_t[len + 1];
	MultiByteToWideChar(CP_ACP, 0, c, strlen(c), m_wchar, len);
	m_wchar[len] = '\0';
	return m_wchar;
}

vector<wchar_t> CharToWcharVector(const char* c)
{
	int len = MultiByteToWideChar(CP_ACP, 0, c, strlen(c), NULL, 0);
	vector<wchar_t> m_wchar(len + 1);
	MultiByteToWideChar(CP_ACP, 0, c, strlen(c), &m_wchar[0], len);
	m_wchar[len] = '\0';
	return m_wchar;
}

struct CategoryInfo
{
	std::wstring name;
	float threshold;
};

class CategoryInfoParser
{
public:
	static HRESULT Parse(const wchar_t* infoFilePath,
		std::map<int, CategoryInfo>& categoryDict)
	{
		std::wifstream file;
		file.open(infoFilePath);

		if (file.bad())
		{
			return E_FAIL;
		}
		else
		{
			categoryDict.clear();

			std::wstring line;
			int idx = 0;
			while (std::getline(file, line))
			{
				std::wstringstream lineStream(line);

				std::wstring categoryName;
				float threshold;

				lineStream >> categoryName >> threshold;

				CategoryInfo info = { categoryName, threshold };
				categoryDict[idx] = info;

				++idx;
			}
		}

		return S_OK;
	}
};


class ObjectDetectionLibraryWrapper
{
private:
	string proposal_modelpath;
	const wchar_t *modelConfigFile;
	const wchar_t *modelCategoriesFile;
	intptr_t handle;
	std::map<int, CategoryInfo> categoryDict;
public:
	ObjectDetectionLibraryWrapper(const char *modelConfigFilePath, const char *modelCategoriesFilePath);
	~ObjectDetectionLibraryWrapper();
	int LoadModel();
	int ProcessImage(const char* inputFileName, const char* outputFileName);
	int ProcessImage(const char* inputFileName, const char* outputFileName, vector<int>& categCount);
};

ObjectDetectionLibraryWrapper::ObjectDetectionLibraryWrapper(const char *modelConfigFilePath, const char *modelCategoriesFilePath)
{
	//const wchar_t *modelConfigFile = L"..\\..\\..\\model\\MultiTasks\\Zeiler_conv5_c1s5n64c2n192\\Cloth_config_detection.txt";
	//const wchar_t *modelCategoriesFile = L"..\\..\\..\\model\\MultiTasks\\Zeiler_conv5_c1s5n64c2n192\\Cloth_detection_label.txt";
	modelConfigFile = CharToWchar(modelConfigFilePath);
	modelCategoriesFile = CharToWchar(modelCategoriesFilePath);
}


int ObjectDetectionLibraryWrapper::LoadModel()
{
	HRESULT hr = S_OK;
	hr = IUCreateHandle(handle);
	if (hr != S_OK) return 1;
	hr = IUSetNumThreads(handle, 1);
	if (hr != S_OK) return 1;

	hr = IULoadModel(handle, IUTask_Detection, modelConfigFile);
	if (hr != S_OK) return 1;


	CategoryInfoParser::Parse(modelCategoriesFile, categoryDict);
	const int class_num = categoryDict.size();
	std::vector<float> pThresholds;
	pThresholds.resize(class_num);
	for (int i = 0; i < class_num; i++)
	{
		pThresholds[i] = categoryDict[i].threshold;
	}
	IUSetDetectionPerClassThresholds(handle, class_num, &pThresholds[0]);

	return 0;
}

int _main()
{
	const char* modelConfigFilePath = "E:\\chenxi\\ImageRecognitionSDK\\model\\MultiTasks\\Zeiler_conv5_c1s5n64c2n192\\Cloth_config_detection.txt";
	const char* modelCategoriesFilePath = "E:\\chenxi\\ImageRecognitionSDK\\model\\MultiTasks\\Zeiler_conv5_c1s5n64c2n192\\Cloth_detection_label.txt";
	ObjectDetectionLibraryWrapper* ptr = new ObjectDetectionLibraryWrapper(modelConfigFilePath, modelCategoriesFilePath);
	ptr->LoadModel();

	return 0;
}