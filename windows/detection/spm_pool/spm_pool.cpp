#include "spm_pool.h"
#include "../../DNNTestLib/NetUtils.h"
#include <vector>
#include <cstdio>
#include <cmath>
#include <omp.h>
using std::vector;


namespace DNNTestLib
{
	void CSpmPooler::Load(const char *buffer, int buffer_size, int &bytes_read)
	{
		class mem_stream_buf : public std::streambuf
		{
		public:
			mem_stream_buf(const char *buffer, int buffer_size)
			{
				this->setg(const_cast<char*>(buffer), const_cast<char*>(buffer), const_cast<char*>(buffer) + buffer_size);
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

	void CSpmPooler::Load(std::istream &stream)
	{
		Clear();

		m_offset0 = ReadData<double>(stream);
		m_offset = ReadData<double>(stream);
		m_step_standard = ReadData<double>(stream);
		m_standard_img_size = ReadData<double>(stream);
		int num_divs = ReadData<int>(stream);

		// add for average pool
		if (num_divs == -1)
		{
			m_is_max_pool = false;
			num_divs = ReadData<int>(stream);
		}

		m_spm_divs.resize(num_divs);
		for (int i = 0; i < num_divs; ++i)
		{
			int divs = ReadData<int>(stream);
			m_spm_divs[i] = divs;
		}

		Setup();
	}

	void CSpmPooler::Load(const std::wstring &file)
	{
		std::ifstream fin(file, std::ios::binary | std::ios::in);
		ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
		Load(fin);
		fin.close();
	}

	void CSpmPooler::Load(const std::string &file)
	{
		std::ifstream fin(file, std::ios::binary | std::ios::in);
		ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
		Load(fin);
		fin.close();
	}

	void CSpmPooler::Setup()
	{
		/// parse spm_divs
		m_total_spm_bins = 0;
		for (int i = 0; i < (int)m_spm_divs.size(); i++)
		{
			m_total_spm_bins += m_spm_divs[i] * m_spm_divs[i];
		}

		if ( m_channels % 4 != 0 )
		{
			std::cout << "spm pool: only support channel % 4 == 0" << std::endl;
			throw std::runtime_error("spm pool: only support channel % 4 == 0");
		}

		const int dim_pooled = std::abs(m_total_spm_bins) * m_channels;
		m_pooled_caches.resize(m_thread_num, NULL);
		for (int iThread = 0; iThread < m_thread_num; ++iThread)
		{
			m_pooled_caches[iThread] = (float*)_aligned_malloc(dim_pooled * sizeof(float), 16);
			if (m_pooled_caches[iThread] == NULL)
			{
				std::cout << "spm pool: malloc error" << std::endl;
				throw std::runtime_error("spm pool: malloc error");
			}
		}
	}

	void CSpmPooler::Clear()
	{
		m_offset = -1;
		m_offset0 = -1;
		m_standard_img_size = -1;
		m_step_standard = -1;
		m_spm_divs.clear();
		m_total_spm_bins = 0;
		m_is_max_pool = true;

		for (int iThread = 0; iThread < (int)m_pooled_caches.size(); ++iThread)
		{
			if (m_pooled_caches[iThread])
			{
				_aligned_free(m_pooled_caches[iThread]);
				m_pooled_caches[iThread] = NULL;	
			}
		}
		m_pooled_caches.clear();
	}
		

	void CSpmPooler::SpatialPooling(__out CPUData *output_rsps, __in const CPUData *response_map, int num_boxes, const double *boxes)
	{
		vector<const CPUData *> response_maps(1);
		response_maps[0] = response_map;
		vector<int> use_response_ids(num_boxes, 0);

		if (m_is_max_pool)
			SpatialPooling_Max(output_rsps, response_maps, num_boxes, boxes, &use_response_ids[0]);
		else
			SpatialPooling_Ave(output_rsps, response_maps, num_boxes, boxes, &use_response_ids[0]);
	}


	// usage:
	// void spatial_pooling(__out cpu_data *output_rsps, __in vector<const cpu_data *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids)
	// all in c++ index
	// response_maps in (num, channel, width, height), num is fastest which always == 1

	void CSpmPooler::SpatialPooling_Max(__out CPUData *output_rsps, __in vector<const CPUData *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids)
	{	
		int feats_num = static_cast<int>(response_maps.size());
		vector<int> rsp_channels(feats_num);
		vector<int> rsp_widths(feats_num), rsp_heights(feats_num);
		vector<const float *> feats(feats_num);

		for (int i = 0; i < feats_num; ++i)
		{
			feats[i] = response_maps[i]->GetDataPointer();
		}

		for (int i = 0; i < feats_num; ++i)
		{
			rsp_channels[i] = response_maps[i]->GetChannels();
			rsp_widths[i] = response_maps[i]->GetWidth();
			rsp_heights[i] = response_maps[i]->GetHeight();
		}

		int channel = rsp_channels[0];
		for (int i = 1; i < feats_num; ++i)
		{
			if (rsp_channels[i] != channel)
			{
				std::cout << "spm pool: currently only support feats with the same channel." << std::endl;
				throw std::runtime_error("spm pool: currently only support feats with the same dim.");
			}
		}

		///// normalize box
		int* normed_boxes = new int[num_boxes * 4];
		NormlizeBoxes(normed_boxes, rsp_widths, rsp_heights, num_boxes, boxes, use_response_ids);
		//////////////////////////////////////////////////////////////////////////

		const int dim_pooled = std::abs(m_total_spm_bins) * channel;

		output_rsps->Resize(num_boxes, dim_pooled, 1, 1, 0, 0, CPUData::CWHN, false);
		float* pooled = output_rsps->GetDataPointer();
		//memset((void*)pooled, 0, sizeof(float) * num_boxes * dim_pooled);

		const int dim_pack = channel / 4;

		const int expend_length = 8;
		const int expend_loop = dim_pack / expend_length;
		const int expend_remain = dim_pack % expend_length;

		if ( channel % 4 != 0 )
		{
			std::cout << "spm pool: only support channel % 4 == 0" << std::endl;
			throw std::runtime_error("spm pool: only support channel % 4 == 0");
		}

		const int thread_num = 1;
		vector<float *> pooled_caches(thread_num);
		for (int iThread = 0; iThread < thread_num; ++iThread)
		{
			pooled_caches[iThread] = (float*)_aligned_malloc(dim_pooled * sizeof(float), 16);
			if (pooled_caches[iThread] == NULL)
			{
				std::cout << "spm_pool:: malloc error" << std::endl;
				throw std::runtime_error("spm_pool:: malloc error");
			}
		}

		//omp_set_num_threads(thread_num);
		int num_boxes_per_thread = int(num_boxes / thread_num) + 1;
		//#pragma omp parallel for ordered schedule(dynamic)
		for (int iThread = 0; iThread < thread_num; ++iThread)
		{
			for (int i = num_boxes_per_thread * iThread; i < min(num_boxes_per_thread * (iThread + 1), num_boxes); i ++)
			{
				float *pooled_cache = pooled_caches[iThread];

				int best_feat = use_response_ids[i];

				const int feats_stride = rsp_widths[best_feat] * channel;

				memset((void*)pooled_cache, 0, sizeof(float) * dim_pooled);

				const int* box_norm = normed_boxes + i * 4;

				const int boxwidth = box_norm[3] - box_norm[1] + 1;
				const int boxheight = box_norm[2] - box_norm[0] + 1;

				const int rspwidth = rsp_widths[best_feat];
				const int rspheight = rsp_heights[best_feat];

				float* pooled_this_div_cache = pooled_cache;
				for (int lv = 0; lv < (int)m_spm_divs.size(); lv ++)
				{
					const float bin_divs = (float)m_spm_divs[lv];

					//const float wunit = boxwidth / (float)bin_divs;
					//const float hunit = boxheight / (float)bin_divs;

					for (int yy = 0; yy < bin_divs; yy ++)
					{
						//int y_start = (int)floor((yy - 1) * hunit) + 1;
						//int y_end = (int)ceil(yy * hunit);

						int y_start = max(0, min(rspheight, (int)floor(yy / bin_divs * boxheight) + box_norm[0]));
						int y_end = max(0, min(rspheight, (int)ceil((yy + 1) / bin_divs * boxheight) + box_norm[0]));

						//assert((y_start >= 0) && (y_end <= height));

						for (int xx = 0; xx < bin_divs; xx ++)
						{
							int x_start = max(0, min(rspwidth, (int)floor(xx / bin_divs * boxwidth) + box_norm[1]));
							int x_end = max(0, min(rspwidth, (int)ceil((xx + 1) / bin_divs * boxwidth) + box_norm[1]));

							//assert((x_start >= 0) && (x_end <= width));

							const float* feats_ptr = feats[best_feat] + y_start * feats_stride;
							for (int y = y_start; y < y_end; y ++)
							{
								//const float* feats_this = feats_ptr + x_start * dim;
								const __m128* feats_this_sse = (__m128*)(feats_ptr + x_start * channel);
								__m128* pooled_this_div_sse = (__m128*)pooled_this_div_cache;

								//for (int x = x_start; x < x_end; x++)
								//{
								//	for (int d = 0; d < dim_pack; d ++)
								//		pooled_this_div_sse[d] = _mm_max_ps(pooled_this_div_sse[d], *feats_this_sse ++);
								//}

								for (int x = x_start; x < x_end; x++)
								{
									pooled_this_div_sse = (__m128*)pooled_this_div_cache;
									for (int d = 0; d < expend_loop; d++)
									{
										pooled_this_div_sse[0] = _mm_max_ps(pooled_this_div_sse[0], feats_this_sse[0]);
										pooled_this_div_sse[1] = _mm_max_ps(pooled_this_div_sse[1], feats_this_sse[1]);
										pooled_this_div_sse[2] = _mm_max_ps(pooled_this_div_sse[2], feats_this_sse[2]);
										pooled_this_div_sse[3] = _mm_max_ps(pooled_this_div_sse[3], feats_this_sse[3]);
										pooled_this_div_sse[4] = _mm_max_ps(pooled_this_div_sse[4], feats_this_sse[4]);
										pooled_this_div_sse[5] = _mm_max_ps(pooled_this_div_sse[5], feats_this_sse[5]);
										pooled_this_div_sse[6] = _mm_max_ps(pooled_this_div_sse[6], feats_this_sse[6]);
										pooled_this_div_sse[7] = _mm_max_ps(pooled_this_div_sse[7], feats_this_sse[7]);
										pooled_this_div_sse += expend_length;
										feats_this_sse += expend_length;
									}
									for (int d = 0; d < expend_remain; d++)
										pooled_this_div_sse[d] = _mm_max_ps(pooled_this_div_sse[d], *feats_this_sse++);
								}
									
								feats_ptr += feats_stride;
							}//y
							pooled_this_div_cache += channel;
						}//xx
					}//yy
				}//lv

				{
					// keep ([channel, width, height], num)
					float *pooled_this_box = pooled + i * dim_pooled;
					float *pooled_this_div_cache = pooled_cache;
					memcpy(pooled_this_box, pooled_this_div_cache, dim_pooled * sizeof(float));
				}

				//{
				//	// trans from ([channel, width, height], num) to ([width, height, channel], num)
				//	float *pooled_this_box = pooled + i * dim_pooled;
				//	float *pooled_this_div = pooled_this_box, *pooled_this_div_cache = pooled_cache;
				//	for (int lv = 0; lv < (int)spm_divs.size(); lv ++)
				//	{
				//		const int bin_divs = spm_divs[lv];
				//		for (int ii = 0; ii < bin_divs * bin_divs; ii ++)
				//		{
				//			for (int d = 0; d < channel; d ++)
				//			{
				//				pooled_this_div[(d * bin_divs * bin_divs + ii)] = pooled_this_div_cache[d];
				//			}
				//			pooled_this_div_cache += channel;
				//		}
				//		pooled_this_div += bin_divs * bin_divs * channel; 
				//	}//lv
				//}
			}//i
		}

		for (int iThread = 0; iThread < thread_num; ++iThread)
			_aligned_free(pooled_caches[iThread]);
		pooled_caches.clear();
		delete[] normed_boxes;

	}

	// usage:
	// void spatial_pooling(__out cpu_data *output_rsps, __in vector<const cpu_data *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids)
	// all in c++ index
	// response_maps in (num, channel, width, height), num is fastest which always == 1

	void CSpmPooler::SpatialPooling_Ave(__out CPUData *output_rsps, __in vector<const CPUData *> response_maps, int num_boxes, const double *boxes, const int *use_response_ids)
	{
		int feats_num = static_cast<int>(response_maps.size());
		vector<int> rsp_channels(feats_num);
		vector<int> rsp_widths(feats_num), rsp_heights(feats_num);
		vector<const float *> feats(feats_num);

		for (int i = 0; i < feats_num; ++i)
		{
			feats[i] = response_maps[i]->GetDataPointer();
		}

		for (int i = 0; i < feats_num; ++i)
		{
			rsp_channels[i] = response_maps[i]->GetChannels();
			rsp_widths[i] = response_maps[i]->GetWidth();
			rsp_heights[i] = response_maps[i]->GetHeight();
		}

		int channel = rsp_channels[0];
		for (int i = 1; i < feats_num; ++i)
		{
			if (rsp_channels[i] != channel)
			{
				std::cout << "spm pool: currently only support feats with the same channel." << std::endl;
				throw std::runtime_error("spm pool: currently only support feats with the same dim.");
			}
		}

		///// normalize box
		int* normed_boxes = new int[num_boxes * 4];
		NormlizeBoxes(normed_boxes, rsp_widths, rsp_heights, num_boxes, boxes, use_response_ids);
		//////////////////////////////////////////////////////////////////////////

		const int dim_pooled = std::abs(m_total_spm_bins) * channel;

		output_rsps->Resize(num_boxes, dim_pooled, 1, 1, 0, 0, CPUData::CWHN, false);
		float* pooled = output_rsps->GetDataPointer();
		//memset((void*)pooled, 0, sizeof(float) * num_boxes * dim_pooled);

		const int dim_pack = channel / 4;

		const int expend_length = 8;
		const int expend_loop = dim_pack / expend_length;
		const int expend_remain = dim_pack % expend_length;

		if (channel % 4 != 0)
		{
			std::cout << "spm pool: only support channel % 4 == 0" << std::endl;
			throw std::runtime_error("spm pool: only support channel % 4 == 0");
		}

		const int thread_num = 1;
		vector<float *> pooled_caches(thread_num);
		for (int iThread = 0; iThread < thread_num; ++iThread)
		{
			pooled_caches[iThread] = (float*)_aligned_malloc(dim_pooled * sizeof(float), 16);
			if (pooled_caches[iThread] == NULL)
			{
				std::cout << "spm_pool:: malloc error" << std::endl;
				throw std::runtime_error("spm_pool:: malloc error");
			}
		}

		//omp_set_num_threads(thread_num);
		int num_boxes_per_thread = int(num_boxes / thread_num) + 1;
		//#pragma omp parallel for ordered schedule(dynamic)
		for (int iThread = 0; iThread < thread_num; ++iThread)
		{
			for (int i = num_boxes_per_thread * iThread; i < min(num_boxes_per_thread * (iThread + 1), num_boxes); i++)
			{
				float *pooled_cache = pooled_caches[iThread];

				int best_feat = use_response_ids[i];

				const int feats_stride = rsp_widths[best_feat] * channel;

				memset((void*)pooled_cache, 0, sizeof(float)* dim_pooled);

				const int* box_norm = normed_boxes + i * 4;

				const int boxwidth = box_norm[3] - box_norm[1] + 1;
				const int boxheight = box_norm[2] - box_norm[0] + 1;

				const int rspwidth = rsp_widths[best_feat];
				const int rspheight = rsp_heights[best_feat];

				float* pooled_this_div_cache = pooled_cache;
				for (int lv = 0; lv < (int)m_spm_divs.size(); lv++)
				{
					const float bin_divs = (float)m_spm_divs[lv];

					//const float wunit = boxwidth / (float)bin_divs;
					//const float hunit = boxheight / (float)bin_divs;

					for (int yy = 0; yy < bin_divs; yy++)
					{
						//int y_start = (int)floor((yy - 1) * hunit) + 1;
						//int y_end = (int)ceil(yy * hunit);

						int y_start = max(0, min(rspheight, (int)floor(yy / bin_divs * boxheight) + box_norm[0]));
						int y_end = max(0, min(rspheight, (int)ceil((yy + 1) / bin_divs * boxheight) + box_norm[0]));

						//assert((y_start >= 0) && (y_end <= height));

						for (int xx = 0; xx < bin_divs; xx++)
						{
							int x_start = max(0, min(rspwidth, (int)floor(xx / bin_divs * boxwidth) + box_norm[1]));
							int x_end = max(0, min(rspwidth, (int)ceil((xx + 1) / bin_divs * boxwidth) + box_norm[1]));

							//assert((x_start >= 0) && (x_end <= width));

							const float* feats_ptr = feats[best_feat] + y_start * feats_stride;
							for (int y = y_start; y < y_end; y++)
							{
								//const float* feats_this = feats_ptr + x_start * dim;
								const __m128* feats_this_sse = (__m128*)(feats_ptr + x_start * channel);
								__m128* pooled_this_div_sse = (__m128*)pooled_this_div_cache;

								//for (int x = x_start; x < x_end; x++)
								//{
								//	for (int d = 0; d < dim_pack; d ++)
								//		pooled_this_div_sse[d] = _mm_add_ps(pooled_this_div_sse[d], *feats_this_sse ++);
								//}

								for (int x = x_start; x < x_end; x++)
								{
									pooled_this_div_sse = (__m128*)pooled_this_div_cache;
									for (int d = 0; d < expend_loop; d++)
									{
										pooled_this_div_sse[0] = _mm_add_ps(pooled_this_div_sse[0], feats_this_sse[0]);
										pooled_this_div_sse[1] = _mm_add_ps(pooled_this_div_sse[1], feats_this_sse[1]);
										pooled_this_div_sse[2] = _mm_add_ps(pooled_this_div_sse[2], feats_this_sse[2]);
										pooled_this_div_sse[3] = _mm_add_ps(pooled_this_div_sse[3], feats_this_sse[3]);
										pooled_this_div_sse[4] = _mm_add_ps(pooled_this_div_sse[4], feats_this_sse[4]);
										pooled_this_div_sse[5] = _mm_add_ps(pooled_this_div_sse[5], feats_this_sse[5]);
										pooled_this_div_sse[6] = _mm_add_ps(pooled_this_div_sse[6], feats_this_sse[6]);
										pooled_this_div_sse[7] = _mm_add_ps(pooled_this_div_sse[7], feats_this_sse[7]);
										pooled_this_div_sse += expend_length;
										feats_this_sse += expend_length;
									}
									for (int d = 0; d < expend_remain; d++)
										pooled_this_div_sse[d] = _mm_add_ps(pooled_this_div_sse[d], *feats_this_sse++);
								}

								feats_ptr += feats_stride;
							}//y
							//sum 2 ave
							{
								int blockSize = (y_end - y_start) * (x_end - x_start);
								if (blockSize > 0)
								{
									__m128 blockSizeInv = _mm_set1_ps(1.f / blockSize);
									__m128* pooled_this_div_sse = (__m128*)pooled_this_div_cache;
									for (int d = 0; d < expend_loop; d++)
									{
										pooled_this_div_sse[0] = _mm_mul_ps(pooled_this_div_sse[0], blockSizeInv);
										pooled_this_div_sse[1] = _mm_mul_ps(pooled_this_div_sse[1], blockSizeInv);
										pooled_this_div_sse[2] = _mm_mul_ps(pooled_this_div_sse[2], blockSizeInv);
										pooled_this_div_sse[3] = _mm_mul_ps(pooled_this_div_sse[3], blockSizeInv);
										pooled_this_div_sse[4] = _mm_mul_ps(pooled_this_div_sse[4], blockSizeInv);
										pooled_this_div_sse[5] = _mm_mul_ps(pooled_this_div_sse[5], blockSizeInv);
										pooled_this_div_sse[6] = _mm_mul_ps(pooled_this_div_sse[6], blockSizeInv);
										pooled_this_div_sse[7] = _mm_mul_ps(pooled_this_div_sse[7], blockSizeInv);
										pooled_this_div_sse += expend_length;
									}
									for (int d = 0; d < expend_remain; d++)
										pooled_this_div_sse[d] = _mm_mul_ps(pooled_this_div_sse[d], blockSizeInv);
								}	
							}
							pooled_this_div_cache += channel;
						}//xx
					}//yy
				}//lv

				{
					// keep ([channel, width, height], num)
					float *pooled_this_box = pooled + i * dim_pooled;
					float *pooled_this_div_cache = pooled_cache;
					memcpy(pooled_this_box, pooled_this_div_cache, dim_pooled * sizeof(float));
				}

				//{
				//	// trans from ([channel, width, height], num) to ([width, height, channel], num)
				//	float *pooled_this_box = pooled + i * dim_pooled;
				//	float *pooled_this_div = pooled_this_box, *pooled_this_div_cache = pooled_cache;
				//	for (int lv = 0; lv < (int)spm_divs.size(); lv ++)
				//	{
				//		const int bin_divs = spm_divs[lv];
				//		for (int ii = 0; ii < bin_divs * bin_divs; ii ++)
				//		{
				//			for (int d = 0; d < channel; d ++)
				//			{
				//				pooled_this_div[(d * bin_divs * bin_divs + ii)] = pooled_this_div_cache[d];
				//			}
				//			pooled_this_div_cache += channel;
				//		}
				//		pooled_this_div += bin_divs * bin_divs * channel; 
				//	}//lv
				//}
			}//i
		}

		for (int iThread = 0; iThread < thread_num; ++iThread)
			_aligned_free(pooled_caches[iThread]);
		pooled_caches.clear();
		delete[] normed_boxes;

	}

	//  box and box_norm both in [l, t, r, b] * num_boxes
	void CSpmPooler::NormlizeBoxes(__out int *normed_boxes, __in const vector<int> &rsp_widths, const vector<int> &rsp_heights, int num_boxes, const double *boxes, const int *use_response_ids)
	{
		for (int i = 0; i < num_boxes; i ++)
		{
			int best_rsp = use_response_ids[i];

			const double* box = boxes + i * 4;

			int* box_norm = normed_boxes + i * 4;

			// to Matlab index
			double x0 = box[0] + 1;
			double y0 = box[1] + 1;
			double x1 = box[2] + 1;
			double y1 = box[3] + 1;

			double width = x1 - x0 + 1;
			double height = y1 - y0 + 1;
			double min_box_sz = max(0.0, min(width, height));
			double offset_norm = m_offset * min_box_sz / m_standard_img_size;

			//// For standard rpn mapping, to c++ index
			//int x0_norm = int(round((x0 - 1 - m_offset0 + offset_norm) / m_step_standard));
			//int y0_norm = int(round((y0 - 1 - m_offset0 + offset_norm) / m_step_standard ));

			//int x1_norm = int(round((x1 - 1 - m_offset0 - offset_norm) / m_step_standard ));
			//int y1_norm = int(round((y1 - 1 - m_offset0 - offset_norm) / m_step_standard ));

			int x0_norm = int(floor((x0 - m_offset0 + offset_norm) / m_step_standard + 0.5));
			int y0_norm = int(floor((y0 - m_offset0 + offset_norm) / m_step_standard + 0.5));

			int x1_norm = int(ceil((x1 - m_offset0 - offset_norm) / m_step_standard - 0.5));
			int y1_norm = int(ceil((y1 - m_offset0 - offset_norm) / m_step_standard - 0.5));

			if (x0_norm > x1_norm)
			{
				x0_norm = (x0_norm + x1_norm) / 2;
				x1_norm = x0_norm;
			}

			if (y0_norm > y1_norm)
			{
				y0_norm = (y0_norm + y1_norm) / 2;
				y1_norm = y0_norm;
			}

			box_norm[0] = y0_norm; // top // must not change
			box_norm[2] = y1_norm; // bottom // must not change

			box_norm[1] = x0_norm; // left // must not change
			box_norm[3] = x1_norm; //right // must not change

			//box_norm[0] = min(rsp_heights[best_rsp]-1, max(0, y0_norm)); // top // must not change
			//box_norm[2] = min(rsp_heights[best_rsp]-1, max(0, y1_norm)); // bottom // must not change

			//box_norm[1] = min(rsp_widths[best_rsp]-1, max(0, x0_norm));  // left // must not change
			//box_norm[3] = min(rsp_widths[best_rsp]-1, max(0, x1_norm));  //right // must not change
		}
	}
}