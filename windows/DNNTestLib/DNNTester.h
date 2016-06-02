
#pragma once

#include "Layer.h"
#include "DataLayer.h"
#include "CPUData.h"

#include <iostream>
#include <istream>
#include <vector>
#include <string>
#include <map>

namespace DNNTestLib
{
	class DNNTester
	{
	public:
		DNNTester();
		~DNNTester();

		void TestImage(const unsigned char *source, int width, int height, 
			int stride, int channels, void *param, __out float *output, int buffer_size);
		void TestImage(const unsigned char *source, int width, int height,
			int stride, int channels, void *param, __out CPUData *output);
		void TestCPUData(const CPUData *input, CPUData *output);
		//add by hg
		void DNNTester::TestCPUData(const CPUData *input, CPUData *output1, CPUData *output2);

		void TestCPUDataDBG(const CPUData *input, CPUData *output);
		void GetOutputDim(int *channels, int *width, int *height) const;
		void GetInputDim(int *channels, int *width, int *height) const;
		void GetLayerNum(int *num) const;
		const char * GetLayerName(const int layer_idx) const;
		void GetLayerResponse(std::string &layer_name, float *buffer, int buffer_size);
		const CPUData * GetLayerResponse(std::string &layer_name);
		void GetLayerDim(std::string &layer_name, int *channels, int *width, int *height, int *num);
		void Dump(std::ostream &stream);

		void Load(const std::wstring &file);
		void Load(const std::string &file);
		void Load(std::istream &stream);
		void Load(const char *buffer, int buffer_size, int &bytes_read);
		void Clear();
		
	public://private:
		// attributes
		std::vector<Layer *> m_layers;
		DataLayer *m_data_layer;
		std::map<std::string, Layer*> m_name2Layer;
		bool m_init;
	};

}