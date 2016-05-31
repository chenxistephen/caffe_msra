#include "JointBayesianLayer.h"
#include "NetUtils.h"

using namespace DNNTestLib;

JointBayesianLayer::JointBayesianLayer(std::istream &stream, std::vector<std::string> &input_links)
: Layer(stream, input_links)
{
	m_sqrt = ReadData<bool>(stream);
	m_mean = CPUData::CreateFromDump(stream);
	m_A = CPUData::CreateFromDump(stream);
	m_G = CPUData::CreateFromDump(stream);
	m_A_buffer = CPUData::Create(1, m_A->GetChannels(), 1, 1, 0, 0);

	assert(m_mean->GetNum() == 1);
	assert(m_mean->GetWidth() == 1);
	assert(m_mean->GetHeight() == 1);
	assert(m_A->GetWidth() == 1);
	assert(m_A->GetHeight() == 1);
	assert(m_G->GetWidth() == 1);
	assert(m_G->GetHeight() == 1);
	assert(m_mean->GetChannels() == m_A->GetChannels());
	assert(m_mean->GetChannels() == m_G->GetChannels());
}


JointBayesianLayer::~JointBayesianLayer()
{
	SaveRelease(m_mean);
	SaveRelease(m_A);
	SaveRelease(m_G);
	SaveRelease(m_A_buffer);
}


void JointBayesianLayer::Dump(std::ostream &stream)
{
	Layer::Dump(stream);
	WriteData(stream, m_sqrt);
	m_mean->Dump(stream);
	m_A->Dump(stream);
	m_G->Dump(stream);
}


void JointBayesianLayer::InitActs()
{
	assert(GetInputNumber() == 1);

	CPUData *input = m_inputs[0]->GetActs();
	assert(input->GetNum() == 1);
	assert(input->GetWidth() == 1);
	assert(input->GetHeight() == 1);
	assert(input->GetChannels() == m_mean->GetChannels());
	assert(m_min_pad_x == 0 && m_min_pad_y == 0);

	if (m_layerActs == NULL)
		m_layerActs = CPUData::Create(1, m_G->GetNum() + 1, 1, 1, 0, 0);
	else
		m_layerActs->Resize(1, m_G->GetNum() + 1, 1, 1, 0, 0);
}


void JointBayesianLayer::Forward()
{
	const float *input = m_inputs[0]->GetActs()->GetDataPointer();
	const float *mean = m_mean->GetDataPointer();
	float *output = m_layerActs->GetDataPointer();
	float *A_buffer = m_A_buffer->GetDataPointer();
	float *G_buffer = output;
	float *p_A = m_A->GetDataPointer();
	float *p_G = m_G->GetDataPointer();

	const int dim_ft = m_mean->GetChannels();
	const int dim_A = m_A->GetNum();
	const int dim_G = m_G->GetNum();
	memset(output, 0, sizeof(float)*(dim_G + 1));
	memset(A_buffer, 0, sizeof(float)*dim_A);

	for (int i = 0; i < dim_ft; i++)
	{
		float v;
		if (m_sqrt)
			v = (input[i] > 0 ? sqrt(input[i]) : -sqrt(-input[i])) - mean[i];
		else
			v = input[i] - mean[i];
		for (int j = 0; j < dim_G; j++, p_G++)
			G_buffer[j] += v * (*p_G);
		for (int j = 0; j < dim_A; j++, p_A++)
			A_buffer[j] += v * (*p_A);
	}
	float f_A = 0;
	for (int i = 0; i < dim_A; i++)
		f_A += A_buffer[i] * A_buffer[i];
	output[dim_G] = f_A;
}




