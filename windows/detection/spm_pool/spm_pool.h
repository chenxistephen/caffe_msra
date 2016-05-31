#ifndef __SPM_POOL_H__
#define __SPM_POOL_H__

#include "../../DNNTestLib/CPUData.h"
#include <vector>
#include <cstdio>
#include <cmath>
#include <omp.h>


namespace DNNTestLib
{
	class CSpmPooler
	{
	public:
		CSpmPooler() :m_offset0(-1), m_offset(-1), m_step_standard(0), m_standard_img_size(0), m_channels(0), m_is_max_pool(true)
		{
			Setup();
		}

		void SpatialPooling(__out CPUData *output_rsps, __in const CPUData *response_map, int num_boxes, const double *boxes);
		void SpatialPooling_Max(__out CPUData *output_rsps, __in std::vector<const CPUData *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids);
		void SpatialPooling_Ave(__out CPUData *output_rsps, __in std::vector<const CPUData *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids);

		void Load(const std::wstring &file);
		void Load(const std::string &file);
		void Load(std::istream &stream);
		void Load(const char *buffer, int buffer_size, int &bytes_read);
		void Clear();

	private:
		void Setup();

		void NormlizeBoxes(__out int *normed_boxes, __in const std::vector<int> &rsp_widths, const std::vector<int> &rsp_heights, int num_boxes, const double *boxes, const int *use_response_ids);

		double m_offset0;
		double m_offset;
		double m_step_standard; 
		double m_standard_img_size;
		std::vector<int> m_spm_divs;
		int m_channels;
		int m_total_spm_bins;

		static const int m_thread_num = 1;
		std::vector<float *> m_pooled_caches;

		bool m_is_max_pool;
	};

	
}

#endif