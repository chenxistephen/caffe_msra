
#include "DNNTester.h"
#include "NetUtils.h"

#include <streambuf>
#include <iostream>
#include <atlconv.h>
#include <new>
#include <set>
#include <map>
#include <time.h>
#include <Windows.h>

using namespace DNNTestLib;

// add by hg
#define ShowTime
#ifdef ShowTime
#include "../ImageRecognitionSDK/HighResolutionTimer.h"
#define TimerStart(x) HighResolutionTimer timer_##x; timer_##x.Reset()
#define TimerShow(x) std::cout << #x##" time = " << timer_##x.Duration() << "ms" << endl
#else
#define TimerStart(x) 
#define TimerShow(x)
#endif


static const std::string dump_header = "conv_net_v1.0";

DNNTester::DNNTester() : m_init(false), m_data_layer(NULL)
{
}

DNNTester::~DNNTester()
{
	Clear();
}

void DNNTester::Clear()
{
	m_init = false;
	for (auto iter = m_layers.begin();
		iter != m_layers.end(); ++iter)
	{
		delete *iter;
	}
	m_layers.clear();
	m_name2Layer.clear();
	m_data_layer = NULL;
}

void DNNTester::Load(const std::wstring &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
	Load(fin);
	fin.close();
}

void DNNTester::Load(const std::string &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
	Load(fin);
	fin.close();
}

void DNNTester::Load(const char *buffer, int buffer_size, int &bytes_read)
{
	class mem_stream_buf : public std::streambuf
	{
	public:
		mem_stream_buf(const char *buffer, int buffer_size)
		{
			this->setg(const_cast<char*>(buffer), const_cast<char*>(buffer), const_cast<char*>(buffer)+buffer_size);
		}

		int bytes_read()
		{
			return static_cast<int>(gptr() - eback());
		}
	};

	mem_stream_buf buf(buffer, buffer_size);
	std::istream stream(&buf);
	Load(stream);
	bytes_read = buf.bytes_read();
}

void DNNTester::Load(std::istream &stream)
{
	Clear();
	ThrowException(ConfigFileError, "Illegal Dump file",
		ReadString(stream) == dump_header);

	// build Layer list
	int n = ReadData<int>(stream);
	std::vector<std::vector<std::string> > prev;
	prev.resize(n);

	m_layers.reserve(n);

	std::map<std::string, int> hash;

	for (int i = 0; i < n; i++)
	{
		Layer *p = Layer::Create(stream, prev[i]);
		m_layers.push_back(p);
		hash[p->GetName()] = i;
		m_name2Layer[p->GetName()] = p;
	}

	m_data_layer = dynamic_cast<DataLayer*>(m_layers[0]);
	ThrowException(ConfigFileError, "First Layer should be a kind of data Layer", m_data_layer != NULL);

	// add links
	for (int i = 0; i < n; i++)
	{
		if (prev[i].size() != 0)
		{
			for (auto iter = prev[i].begin();
				iter != prev[i].end(); ++iter)
			{
				assert(hash.count(*iter) != 0);
				int lk = hash[*iter];
				ThrowException(ConfigFileError, "The input of each Layer should before itself", lk < i);
				m_layers[i]->AddInputLayer(m_layers[lk]);
			}
		}
	}

	// initialize Acts matrices
	for (auto iter = m_layers.rbegin(); iter != m_layers.rend(); iter++)
		(*iter)->InitInputMinPadSize();
	for (auto iter = m_layers.begin(); iter != m_layers.end(); iter++)
	{
		(*iter)->InitActs();
	}
	m_init = true;

	// print network infomation
	for (auto iter = m_layers.begin(); iter != m_layers.end(); iter++)
	{
		auto inputs = (*iter)->GetInputs();
		if ((*iter) == m_data_layer)
		{
			printf("input ");
		}
		for (int i = 0; i < static_cast<int>(inputs.size()); i++)
		{
			printf("%s [%d %d %d %d] ", inputs[i]->GetName().c_str(),
				inputs[i]->GetActs()->GetNum(), inputs[i]->GetActs()->GetChannels(),
				inputs[i]->GetActs()->GetWidth(), inputs[i]->GetActs()->GetHeight());
		}
		printf("-> %s [%d %d %d %d] ", (*iter)->GetName().c_str(),
			(*iter)->GetActs()->GetNum(), (*iter)->GetActs()->GetChannels(),
			(*iter)->GetActs()->GetWidth(), (*iter)->GetActs()->GetHeight());
		if (inputs.size() == 1 && inputs[0]->GetActs() == (*iter)->GetActs())
			printf("(in-place)\n");
		else
			printf("\n");
	}
}

void DNNTester::TestImage(const unsigned char *source, int width, int height,
	int stride, int channels, void *param, __out float *output, int buffer_size)
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);

	m_data_layer->FeedImage(source, width, height, channels, stride, param);
	for (auto iter = m_layers.begin() + 1; iter != m_layers.end(); iter++)
	{
		(*iter)->InitActs();	// in case you would like to change input size every time
		(*iter)->Forward();
	}
	m_layers.back()->GetResponse(output, buffer_size);
	return;
}

void DNNTester::TestImage(const unsigned char *source, int width, int height,
	int stride, int channels, void *param, __out CPUData *output)
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);

	m_data_layer->FeedImage(source, width, height, channels, stride, param);
	for (auto iter = m_layers.begin() + 1; iter != m_layers.end(); iter++)
	{
		//clock_t start = clock();
		(*iter)->InitActs();	// in case you would like to change input size every time
		(*iter)->Forward();
		//printf("%s spend %dms\n", (*iter)->GetName().c_str(), clock() - start);
	}
	const CPUData *result = m_layers.back()->GetActs();
	output->Resize(result->GetNum(), result->GetChannels(), result->GetWidth(), result->GetHeight(), 0, 0);
	output->CopyDataFrom(result);
	return;
}




void DNNTester::TestCPUData(const CPUData *input, CPUData *output)
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);

	if (output == NULL)
		return;
	if (input->GetDataStorageType() != CPUData::CWHN)
	{
		CPUData inputCopy;
		inputCopy.Resize(input->GetNum(), input->GetChannels(), input->GetWidth(), input->GetHeight(), 0, 0, CPUData::CWHN);
		inputCopy.CopyDataFrom(input);
		m_data_layer->FeedCPUData(&inputCopy);
	}
	else
	{
		m_data_layer->FeedCPUData(input);
	}
	for (auto iter = m_layers.begin() + 1; iter != /*m_layers.begin() + 4*/m_layers.end(); iter++)
	{
		(*iter)->InitActs();	// in case you would like to change input size every time
		(*iter)->Forward();
	}
	const CPUData *result = /*(*(m_layers.begin() + 3))->GetActs()*/m_layers.back()->GetActs();
	output->Resize(result->GetNum(), result->GetChannels(), result->GetWidth(), result->GetHeight(), 0, 0);
	output->CopyDataFrom(result);
}

void DNNTester::TestCPUData(const CPUData *input, CPUData *output1, CPUData *output2)
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);

	if (output1 == NULL || output2 == NULL)
		return;
	if (input->GetDataStorageType() != CPUData::CWHN)
	{
		CPUData inputCopy;
		inputCopy.Resize(input->GetNum(), input->GetChannels(), input->GetWidth(), input->GetHeight(), 0, 0, CPUData::CWHN);
		inputCopy.CopyDataFrom(input);
		m_data_layer->FeedCPUData(&inputCopy);
	}
	else
	{
		m_data_layer->FeedCPUData(input);
	}
	for (auto iter = m_layers.begin() + 1; iter != m_layers.end(); iter++)
	{
		//clock_t start = clock();
		(*iter)->InitActs();	// in case you would like to change input size every time
		(*iter)->Forward();
		//printf("%s spend %dms\n", (*iter)->GetName().c_str(), clock() - start);
	}
	const CPUData *result = m_layers.back()->GetActs();
	output1->Resize(result->GetNum(), result->GetChannels(), result->GetWidth(), result->GetHeight(), 0, 0);
	output1->CopyDataFrom(result);

	result = m_layers[m_layers.size() - 2]->GetActs();
	output2->Resize(result->GetNum(), result->GetChannels(), result->GetWidth(), result->GetHeight(), 0, 0);
	output2->CopyDataFrom(result);
}

void DNNTester::TestCPUDataDBG(const CPUData *input, CPUData *output)
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);

	if (output == NULL)
		return;
	m_data_layer->FeedCPUData(input);
	for (auto iter = m_layers.begin() + 1; iter != m_layers.end(); iter++)
	{
		//clock_t start = clock();
		(*iter)->InitActs();	// in case you would like to change input size every time
		(*iter)->Forward();
		//printf("%s spend %dms\n", (*iter)->GetName().c_str(), clock() - start);
	}
	const CPUData *result = m_layers.back()->GetActs();
	output->Resize(result->GetNum(), result->GetChannels(), result->GetWidth(), result->GetHeight(), 0, 0);
	output->CopyDataFrom(result);
}


void DNNTester::GetOutputDim(int *channels, int *width, int *height) const
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);
	m_layers.back()->GetResponseDim(channels, width, height, int());
}

void DNNTester::GetInputDim(int *channels, int *width, int *height) const
{
	ThrowException(InitError, "Conv net doese not init", m_init || m_layers.size() == 0);
	m_data_layer->GetResponseDim(channels, width, height, int());
}

void DNNTester::GetLayerNum(int *num) const
{
	*num = static_cast<int>(m_layers.size());
}

const char *DNNTester::GetLayerName(const int layer_idx) const
{
	ThrowException(ParameterError, "The specified Layer does not exist",
		layer_idx >= 0 && layer_idx < static_cast<int>(m_layers.size()));


	return m_layers[layer_idx]->GetNameRef().c_str();
}

void DNNTester::GetLayerResponse(std::string &layer_name, float *buffer, int buffer_size)
{
	ThrowException(ParameterError, "The specified Layer does not exist",
		m_name2Layer.count(layer_name) == 1);

	ThrowException(ParameterError, "Buffer cannot be NULL",
		buffer != NULL);

	m_name2Layer[layer_name]->GetResponse(buffer, buffer_size);
}

const CPUData * DNNTester::GetLayerResponse(std::string &layer_name)
{
	ThrowException(ParameterError, "The specified Layer does not exist",
		m_name2Layer.count(layer_name) == 1);

	return m_name2Layer[layer_name]->GetActs();
}

void DNNTester::GetLayerDim(std::string &layer_name, int *channels, int *width, int *height, int *num)
{
	ThrowException(ParameterError, "The specified Layer does not exist",
		m_name2Layer.count(layer_name) == 1);

	m_name2Layer[layer_name]->GetResponseDim(channels, width, height, num);
}

void DNNTester::Dump(std::ostream &stream)
{
	WriteString(stream, dump_header);

	// save layers
	int n = (int)m_layers.size();
	WriteData(stream, n);
	for (int i = 0; i < n; i++)
	{
		m_layers[i]->Dump(stream);
	}
}

