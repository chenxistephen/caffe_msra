#pragma once

#include "Layer.h"

namespace DNNTestLib
{
	class JointBayesianLayer : public Layer
	{
	public:
		JointBayesianLayer(std::istream &stream, std::vector<std::string> &input_links);
		virtual ~JointBayesianLayer();
		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual void InitActs();

	protected:
		// properties
		CPUData *m_mean, *m_G, *m_A;
		CPUData *m_A_buffer;
		bool m_sqrt;
	};

}
