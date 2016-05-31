
#include "LayerActs.h"
#include "NetUtils.h"

#include "AcrossChannelNormLayer.h"
#include "ConcatLayer.h"
#include "ConvLayer.h"
#include "DataLayer.h"
#include "DropoutLayer.h"
#include "FaceDataLayer.h"
#include "FCLayer.h"
#include "FlattenLayer.h"
#include "DataLayerMinusMean.h"
#include "JointBayesianLayer.h"
#include "Layer.h"
#include "LocalShareLayer.h"
#include "WithinChannelNormLayer.h"
#include "NeuronLayer.h"
#include "PaddingLayer.h"
#include "PoolLayer.h"
#include "PoolLayerWithPadding.h"
#include "SoftmaxLayer.h"
#include "SkipConvLayer.h"
#include "DataLayerMinusConst.h"

#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <new>

using namespace DNNTestLib;


Layer::Layer(std::istream &stream, std::vector<std::string> &input_links)
{
	m_name = ReadString(stream);
	int n = ReadData<int>(stream);
	for (int i = 0; i<n; i++)
	{
		std::string link = ReadString(stream);
		input_links.push_back(link);
	}

	m_layerActs = NULL;
	m_inputs.clear();
	m_min_pad_x = 0;
	m_min_pad_y = 0;
}

Layer::~Layer()
{
	SaveRelease(m_layerActs);
}

int Layer::GetInputNumber() const
{
	return (int)m_inputs.size();
}

const std::vector<Layer*>& Layer::GetInputs() const
{
	return m_inputs;
}

void Layer::AddInputLayer(Layer *f_layer)
{
	m_inputs.push_back(f_layer);
}

CPUData* Layer::GetActs() const
{
	return m_layerActs;
}

std::string Layer::GetType() const
{
	return m_type;
}

std::string Layer::GetName() const
{
	return m_name;
}

std::string &Layer::GetNameRef()
{
	return m_name;
}

void Layer::Forward()
{
}

int Layer::GetInputPadX()
{
	return 0;
}

int Layer::GetInputPadY()
{
	return 0;
}

void Layer::InitInputMinPadSize()
{
	for (size_t i = 0; i < m_inputs.size(); i++)
	{
		if (m_inputs[i]->m_min_pad_x < this->GetInputPadX())
			m_inputs[i]->m_min_pad_x = this->GetInputPadX();
		if (m_inputs[i]->m_min_pad_y < this->GetInputPadY())
			m_inputs[i]->m_min_pad_y = this->GetInputPadY();
	}
}


Layer *Layer::Create(std::istream &stream, std::vector<std::string> &input_links)
{
	std::string type = ReadString(stream);
	std::string attribute = ReadString(stream);
	Layer *p_layer = NULL;

	if (type == "data")
		p_layer = new DataLayer(stream, input_links);
	else if (type == "data_minus_mean")
		p_layer = new DataLayerMinusMean(stream, input_links);
	else if (type == "data_minus_const")
		p_layer = new DataLayerMinusConst(stream, input_links);
	else if (type == "face_data")
		p_layer = new FaceDataLayer(stream, input_links);
	else if (type == "conv")
		p_layer = new ConvLayer(stream, input_links);
	else if (type == "skipconv")
		p_layer = new SkipConvLayer(stream, input_links);
	else if (type == "local_share")
		p_layer = new LocalShareLayer(stream, input_links);
	else if (type == "pool")
	{
		if (attribute == "max")
			p_layer = new PoolLayer<MaxPoolActs>(stream, input_links);
		else if (attribute == "avg")
			p_layer = new PoolLayer<AvgPoolActs>(stream, input_links);
		else
			ThrowException(ConfigFileError, "Unsupported pooling type");
	}
	else if (type == "pool_with_padding")
	{
		if (attribute == "max")
			p_layer = new PoolLayerWithPadding<MaxPoolActs>(stream, input_links);
		else if (attribute == "avg")
			p_layer = new PoolLayerWithPadding<AvgPoolActs>(stream, input_links);
		else
			ThrowException(ConfigFileError, "Unsupported pooling type");
	}
	else if (type == "neuron")
	{
		if (attribute == "relu")
			p_layer = new NeuronLayer<ReluNeuronActs>(stream, input_links);
		else if (attribute == "sigmoid")
			p_layer = new NeuronLayer<SigmoidNeuronActs>(stream, input_links);
		else
			ThrowException(ConfigFileError, "Unsupported neuron type");
	}
	else if (type == "fc")
		p_layer = new FCLayer(stream, input_links);
	else if (type == "flatten")
		p_layer = new FlattenLayer(stream, input_links);
	else if (type == "concat")
		p_layer = new ConcatLayer(stream, input_links);
	else if (type == "dropout")
		p_layer = new DropoutLayer(stream, input_links);
	else if (type == "padding")
		p_layer = new PaddingLayer(stream, input_links);
	else if (type == "across_channel_norm")
		p_layer = new AcrossChannelNormLayer(stream, input_links);
	else if (type == "lrnorm" || type == "within_channel_norm")
		p_layer = new WithinChannelNormLayer(stream, input_links);
	else if (type == "softmax")
		p_layer = new SoftmaxLayer(stream, input_links);
	else if (type == "joint_bayesian")
		p_layer = new JointBayesianLayer(stream, input_links);
	else
		ThrowException(ConfigFileError, "Unsupported Layer type");

	ThrowException(MemoryAllocationError, "Out of memory", p_layer != NULL);
	p_layer->m_type = type;
	p_layer->m_attribute = attribute;
	return p_layer;
}

void Layer::GetResponse(float *buffer, int buffer_size) const
{
	ThrowException(ParameterError, "The specified Layer contains no responses", GetActs() != NULL);

	m_layerActs->GetData(buffer, buffer_size);
}

void Layer::Dump(std::ostream &stream)
{
	WriteString(stream, m_type);
	WriteString(stream, m_attribute);
	WriteString(stream, m_name);

	WriteData<int>(stream, (int)m_inputs.size());
	for (std::vector<Layer *>::iterator iter = m_inputs.begin();
		iter != m_inputs.end(); ++iter)
	{
		WriteString(stream, (*iter)->GetName());
	}
}

void Layer::GetResponseDim(int *channels, int *width, int *height, int *num) const
{
	assert(m_layerActs != NULL);

	if (channels != NULL)
		*channels = m_layerActs->GetChannels();
	if (width != NULL)
		*width = m_layerActs->GetWidth();
	if (height != NULL)
		*height = m_layerActs->GetHeight();
	if (num != NULL)
		*num = m_layerActs->GetNum();
}

void Layer::InitActs()
{
}

