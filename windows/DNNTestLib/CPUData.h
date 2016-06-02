
#pragma once 

#include <iostream>
#include <fstream>

namespace DNNTestLib
{
	// cpu data is a set of cubes with order (number, channel, width, height)
	class CPUData
	{
	public:	
		enum DataStorageType
		{
			CWHN = 0,
			NCWH = 1,
		};


		CPUData();
		~CPUData();
		void Clear();
		void Resize(int num, int channels, int width, int height, int padding_x, int padding_y, DataStorageType dst = DataStorageType::CWHN, bool clear_mem = true);
		int GetOffset(int num, int channels, int width, int height) const;
		float *GetDataPointer();
		const float *GetDataPointer() const;
		float *GetDataPointer(int num, int channels, int width, int height);
		const float *GetDataPointer(int num, int channels, int width, int height) const;
		void SetData(const float *buffer, int buffer_size);
		void GetData(float *buffer, int buffer_size) const;
		int GetDataSize() const;
		int GetDataSizeIncludePadding() const;
		int GetStrideH() const;
		int GetNum() const;
		int GetChannels() const;
		int GetWidth() const;
		int GetHeight() const;
		int GetPaddingX() const;
		int GetPaddingY() const;
		void CopyDataFrom(const CPUData *p);
		DataStorageType GetDataStorageType() const;
		void PermuteData(DataStorageType dst);

		// save & Load
		void Dump(std::ostream &stream) const;
		static CPUData* Create(int num, int channels, int width, int height, int padding_x = 0, int padding_y = 0, DataStorageType dst = DataStorageType::CWHN);
		static CPUData* CreateFromDump(std::istream &stream);

	public://private:
		CPUData(const CPUData &);
		const CPUData& operator = (const CPUData &);
		int Offset(int num, int channels, int width, int height) const;
		
		// attributes
		float *m_data;
		int m_num;
		int m_channels;
		int m_width;
		int m_height;
		int m_padding_x;
		int m_padding_y;
		int m_capacity;
		DataStorageType m_dst;
	};

}