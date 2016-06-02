
#pragma once

#include "CPUData.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

/*
call order
1. Create
2. GetInputPadX, GetInputPadY
3. InitActs
4. FeedImage
5. Forward
*/


namespace DNNTestLib
{
	class Layer
	{
	public:
		static Layer *Create(std::istream &stream, std::vector<std::string> &input_links);
		int GetInputNumber() const;
		const std::vector<Layer*>& GetInputs() const;
		CPUData* GetActs() const;
		std::string GetType() const;
		std::string GetName() const;
		std::string &Layer::GetNameRef();
		void GetResponse(float *buffer, int buffer_size) const;
		void GetResponseDim(int *channels, int *width, int *height, int *num) const;
		void InitInputMinPadSize();
		void AddInputLayer(Layer *f_layer); // add a former Layer

		virtual void Dump(std::ostream &stream);
		virtual void Forward();
		virtual int GetInputPadX();  // the required padding size of input layers.
		virtual int GetInputPadY();
		virtual void InitActs();
		// Added by Stephen Chen
		virtual CPUData* GetWeight() const { return NULL; }
		virtual CPUData* GetBias() const { return NULL; }
		virtual void CopyWeight(float *buffer, int buffer_size) const{}
		virtual void CopyBias(float *buffer, int buffer_size) const{}
		/////////////////////////////////////
		virtual ~Layer();

	protected:
		Layer(std::istream &stream, std::vector<std::string> &input_links);

	protected:
		std::string m_name;
		std::string m_type;
		std::string m_attribute;
		std::vector<Layer *> m_inputs; // input layers 
		CPUData *m_layerActs;
		int m_min_pad_x;  // the minimal padding size of this Layer
		int m_min_pad_y;
	};


}