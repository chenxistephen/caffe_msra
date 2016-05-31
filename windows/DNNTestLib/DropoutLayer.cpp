#include "DropoutLayer.h"
#include "NetUtils.h"

using namespace DNNTestLib;

DropoutLayer::DropoutLayer(std::istream &stream, std::vector<std::string> &input_links)
	: Layer(stream, input_links)
{
	m_dropout_ratio = ReadData<float>(stream);
}

DropoutLayer::~DropoutLayer()
{
	m_layerActs = NULL;
}

void DropoutLayer::InitActs()
{
	assert(m_inputs[0]->GetActs()->GetDataStorageType() == CPUData::CWHN);
	assert(GetInputNumber() == 1);
	m_layerActs = m_inputs[0]->GetActs();
}

void DropoutLayer::Forward()
{
	float *data = m_layerActs->GetDataPointer(0, 0, -m_layerActs->GetPaddingX(), -m_layerActs->GetPaddingY());
	const int count = m_layerActs->GetDataSizeIncludePadding();
	float scale = 1.0f - m_dropout_ratio;

	for (int i = 0; i < count; i++)
		data[i] *= scale;

}
