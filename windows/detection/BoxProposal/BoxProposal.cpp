#include "BoxProposal.h"
#include "../../DNNTestLib/NetUtils.h"
#include "../../ImageRecognitionSDK/nms.hpp"

#include <algorithm>
#include <time.h>

namespace DNNTestLib
{
	CBoxProposal::CBoxProposal()
	{
	}

	CBoxProposal::~CBoxProposal()
	{
		Release();
	}

	void CBoxProposal::Release()
	{
		m_boxPredictInverseFun = NULL;
		m_getBoxesFun = NULL;
		m_centroids.clear();
		m_net.Clear();
		m_sppPooler.Clear();
		m_poolData.Clear();
		m_outputData_boxes.Clear();
		m_outputData_confidence.Clear();
		m_outputData.Clear();
		m_srcBoxes.clear();
	}

	void CBoxProposal::LoadModel(const char *buffer, int buffer_size, int &bytes_read)
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
		LoadModel(stream);
		bytes_read = buf.bytes_read();
	}

	void CBoxProposal::LoadModel(std::istream &stream)
	{
		// load spp pooler 
		m_sppPooler.Load(stream);

		// load centroids/anchors and box proposal methond
		int numCentroids = ReadData<int>(stream);
		// determine box proposal method(rpn has written an additional "-1" to the model)
		if (numCentroids == -1)
		{
			numCentroids = ReadData<int>(stream);
			m_boxProposalMethod = BoxProposal_rpn;
			m_getBoxesFun = &DNNTestLib::CBoxProposal::GetBoxes_rpn;
		}
		else if (numCentroids > 0)
		{
			m_boxProposalMethod = BoxProposal_multibox;
			m_getBoxesFun = &DNNTestLib::CBoxProposal::GetBoxes_multibox;
		}			
		else            		   
			ThrowException(ParameterError, "Undefined proposal method");

		m_centroids.resize(numCentroids);		
		for (int i = 0; i < numCentroids; ++i)
		{
			m_centroids[i].left = ReadData<float>(stream);
			m_centroids[i].top = ReadData<float>(stream);
			m_centroids[i].right = ReadData<float>(stream);
			m_centroids[i].bottom = ReadData<float>(stream);
		}

		// load transform method
		const int transMethod = ReadData<int>(stream);
		switch (transMethod)
		{
		case 0:
			m_boxTransformType = BoxTransform_Delta;
			m_boxPredictInverseFun = &DNNTestLib::CBoxProposal::BoxPredictInverse_delta;
			break;
		case 10:
			m_boxTransformType = BoxTransform_bbox;
			m_boxPredictInverseFun = &DNNTestLib::CBoxProposal::BoxPredictInverse_bbox;
			break;
		default:
			ThrowException(ModelPackageFileError, "Undefined boxTransformType");
			break;
		}

		// load net
		m_net.Load(stream);

		// sort m_centroids by size for multibox proposal method.
		if (m_boxProposalMethod == BoxProposal_multibox)
		{
			struct RectRealWithSize
			{
				float left, top, right, bottom;
				float size;
			};
			std::vector<RectRealWithSize> centroidsWithSize(m_centroids.size());
			for (int i = 0; i < static_cast<int>(centroidsWithSize.size()); ++i)
			{
				centroidsWithSize[i].left = m_centroids[i].left;
				centroidsWithSize[i].top = m_centroids[i].top;
				centroidsWithSize[i].right = m_centroids[i].right;
				centroidsWithSize[i].bottom = m_centroids[i].bottom;

				centroidsWithSize[i].size = (centroidsWithSize[i].right - centroidsWithSize[i].left) * (centroidsWithSize[i].bottom - centroidsWithSize[i].top);
				assert(centroidsWithSize[i].size > 0);
			}
			std::sort(centroidsWithSize.begin(), centroidsWithSize.end(), [](RectRealWithSize &i, RectRealWithSize &j){return i.size > j.size; });
			for (int i = 0; i < static_cast<int>(centroidsWithSize.size()); ++i)
			{
				m_centroids[i].left = centroidsWithSize[i].left;
				m_centroids[i].top = centroidsWithSize[i].top;
				m_centroids[i].right = centroidsWithSize[i].right;
				m_centroids[i].bottom = centroidsWithSize[i].bottom;
			}
		}
	}

	void CBoxProposal::LoadModel(const std::wstring &file)
	{
		std::ifstream fin(file, std::ios::binary | std::ios::in);
		ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
		LoadModel(fin);
		fin.close();
	}

	void CBoxProposal::LoadModel(const std::string &file)
	{
		std::ifstream fin(file, std::ios::binary | std::ios::in);
		ThrowException(ModelPackageFileError, "Cannot open model file", fin.is_open());
		LoadModel(fin);
		fin.close();
	}

	void CBoxProposal::GetBoxes(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
		const int iMinBatchSize /*= 200*/)
	{
		(this->*m_getBoxesFun)(rgBoxes, rspMap, iImgWidth, iImgHeight, iMinBatchSize);
	}

	void CBoxProposal::GetBoxes_rpn(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
		const int iMinBatchSize /*= 200*/)
	{

		const int rspChannels = rspMap.GetChannels();
		const int rspWidth = rspMap.GetWidth();
		const int rspHeight = rspMap.GetHeight();

		const int rpn_stride = ((int)(iImgWidth / rspWidth / 4) + 1) * 4;

		// rpn - conv_proposal
		m_net.TestCPUData(&rspMap, &m_outputData_boxes, &m_outputData_confidence);
		assert(m_outputData_boxes.GetPaddingX() == 0 && m_outputData_boxes.GetPaddingY() == 0);
		assert(m_outputData_confidence.GetPaddingX() == 0 && m_outputData_confidence.GetPaddingY() == 0);

		std::vector<ProposalBox> rpn_boxes;
		const int numAnchors = static_cast<int>(m_centroids.size());
		const int numBoxes = numAnchors * rspWidth * rspHeight;
		rpn_boxes.resize(numBoxes);
		m_srcBoxes.resize(numBoxes);

		// get local anchors
		for (int i = 0; i < rspHeight; ++i)
		{
			for (int j = 0; j < rspWidth; ++j)
			{
				for (int k = 0; k < numAnchors; ++k)
				{
					m_srcBoxes[i * rspWidth * numAnchors + j * numAnchors + k].left =
						m_centroids[k].left + j * rpn_stride;
					m_srcBoxes[i * rspWidth * numAnchors + j * numAnchors + k].top =
						m_centroids[k].top + i * rpn_stride;
					m_srcBoxes[i * rspWidth * numAnchors + j * numAnchors + k].right =
						m_centroids[k].right + j * rpn_stride;
					m_srcBoxes[i * rspWidth * numAnchors + j * numAnchors + k].bottom =
						m_centroids[k].bottom + i * rpn_stride;
				}
			}
		}

		// get proposal boxes and confidence
		const float *pOutput_boxes      = m_outputData_boxes.GetDataPointer();
		const float *pOutput_confidence = m_outputData_confidence.GetDataPointer();
		for (int i = 0; i < numBoxes; ++i)
		{
			rpn_boxes[i].box.left   = pOutput_boxes[i * 4 + 0];
			rpn_boxes[i].box.top    = pOutput_boxes[i * 4 + 1];
			rpn_boxes[i].box.right  = pOutput_boxes[i * 4 + 2];
			rpn_boxes[i].box.bottom = pOutput_boxes[i * 4 + 3];
			rpn_boxes[i].confidence = 1 / (1 + exp(pOutput_confidence[static_cast<int>(i/9) * 18 + i % 9] 
				                                   - pOutput_confidence[static_cast<int>(i/9) * 18 + i % 9 + 9]));

			(this->*m_boxPredictInverseFun)(rpn_boxes[i].box, rpn_boxes[i].box, m_srcBoxes[i]);

			// C index.
			rpn_boxes[i].box.left = max(1.f, min(iImgWidth, rpn_boxes[i].box.left)) - 1;
			rpn_boxes[i].box.top = max(1.f, min(iImgHeight, rpn_boxes[i].box.top)) - 1;
			rpn_boxes[i].box.right = max(1.f, min(iImgWidth, rpn_boxes[i].box.right)) - 1;
			rpn_boxes[i].box.bottom = max(1.f, min(iImgHeight, rpn_boxes[i].box.bottom)) - 1;
		}

		// sort and get top m_iMaxEvalutedCandidate for nms		
		std::sort(rpn_boxes.begin(), rpn_boxes.end(), [](ProposalBox &i, ProposalBox &j){return i.confidence > j.confidence; });
		rpn_boxes.resize(min(m_iMaxEvalutedCandidate, static_cast<int>(rpn_boxes.size())));
		// nms
		std::vector<double> pre_nms_boxes(static_cast<int>(rpn_boxes.size() * 4));
		std::vector<double> pre_nms_confidence(static_cast<int>(rpn_boxes.size()));
		for (int i = 0; i < rpn_boxes.size(); ++i) 
		{
			pre_nms_boxes[i * 4]     = rpn_boxes[i].box.left;
			pre_nms_boxes[i * 4 + 1] = rpn_boxes[i].box.top;
			pre_nms_boxes[i * 4 + 2] = rpn_boxes[i].box.right;
			pre_nms_boxes[i * 4 + 3] = rpn_boxes[i].box.bottom;
			pre_nms_confidence[i] = rpn_boxes[i].confidence;
		}
		std::vector<int> nms_winner;
		const double rpn_NMSThreshold = 0.7;
		nms(static_cast<int>(rpn_boxes.size()), &pre_nms_boxes[0], &pre_nms_confidence[0], 1, rpn_NMSThreshold, nms_winner);

		// get final proposal boxes, drop more than m_iMaxOuputNum ones
		const int num_proposal = min(m_iMaxOuputNum, static_cast<int>(nms_winner.size()));
		rgBoxes.resize(num_proposal);
		for (int i = 0; i < num_proposal; ++i)
		{
			rgBoxes[i] = rpn_boxes[nms_winner[i]];
		}
	}

	void CBoxProposal::GetBoxes_multibox(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
		const int iMinBatchSize /*= 200*/)
	{		
		std::vector<double> srcBoxes(iMinBatchSize * 4);
		// only evaluate m_iMaxEvalutedCandidate candidates
		int numBoxes = max(1, min(static_cast<int>(m_centroids.size()), m_iMaxEvalutedCandidate));
		rgBoxes.resize(numBoxes);
		m_srcBoxes.resize(numBoxes);

		int inputChannels = 0;
		int inputWidth = 0;
		int inputHeight = 0;
		m_net.GetInputDim(&inputChannels, &inputWidth, &inputHeight);

		for (int i = 0; i < numBoxes; ++i)
		{
			m_srcBoxes[i].left = m_centroids[i].left * (iImgWidth - 1);
			m_srcBoxes[i].top = m_centroids[i].top * (iImgHeight - 1);
			m_srcBoxes[i].right = m_centroids[i].right * (iImgWidth - 1);
			m_srcBoxes[i].bottom = m_centroids[i].bottom * (iImgHeight - 1);
		}

		for (int i = 0; i < static_cast<int>(ceil(numBoxes / (float)iMinBatchSize)); ++i)
		{
			int sampleIdxStart = i * iMinBatchSize;
			int sampleNum = min(numBoxes - sampleIdxStart, iMinBatchSize);

			for (int j = 0; j < sampleNum; ++j)
			{
				srcBoxes[j * 4] = m_srcBoxes[j + sampleIdxStart].left;
				srcBoxes[j * 4 + 1] = m_srcBoxes[j + sampleIdxStart].top;
				srcBoxes[j * 4 + 2] = m_srcBoxes[j + sampleIdxStart].right;
				srcBoxes[j * 4 + 3] = m_srcBoxes[j + sampleIdxStart].bottom;
			}

			// spp
			m_sppPooler.SpatialPooling(&m_poolData, &rspMap, sampleNum, &srcBoxes[0]);
			m_poolData.Resize(sampleNum, inputChannels, inputWidth, inputHeight, 0, 0, DNNTestLib::CPUData::CWHN, false);

			// fc
			m_net.TestCPUData(&m_poolData, &m_outputData);
			assert(m_outputData.GetPaddingX() == 0 && m_outputData.GetPaddingY() == 0);

			const float *pOutput = m_outputData.GetDataPointer();
			for (int j = 0; j < sampleNum; ++j)
			{
				rgBoxes[j + sampleIdxStart].box.left = pOutput[j * 5 + 0];
				rgBoxes[j + sampleIdxStart].box.top = pOutput[j * 5 + 1];
				rgBoxes[j + sampleIdxStart].box.right = pOutput[j * 5 + 2];
				rgBoxes[j + sampleIdxStart].box.bottom = pOutput[j * 5 + 3];
				rgBoxes[j + sampleIdxStart].confidence = pOutput[j * 5 + 4];
			}
		}

		for (int i = 0; i < static_cast<int>(rgBoxes.size()); ++i)
		{
			(this->*m_boxPredictInverseFun)(rgBoxes[i].box, rgBoxes[i].box, m_srcBoxes[i]);

			rgBoxes[i].box.left = max(0.f, min(iImgWidth - 1, rgBoxes[i].box.left));
			rgBoxes[i].box.top = max(0.f, min(iImgHeight - 1, rgBoxes[i].box.top));
			rgBoxes[i].box.right = max(0.f, min(iImgWidth - 1, rgBoxes[i].box.right));
			rgBoxes[i].box.bottom = max(0.f, min(iImgHeight - 1, rgBoxes[i].box.bottom));
		}

		// sort and drop more than m_iMaxOuputNum ones
		std::sort(rgBoxes.begin(), rgBoxes.end(), [](ProposalBox &i, ProposalBox &j){return i.confidence > j.confidence; });
		rgBoxes.resize(min(m_iMaxOuputNum, static_cast<int>(rgBoxes.size())));
	}

	void CBoxProposal::SetMaxEvalutedCandidate(int n)
	{
		m_iMaxEvalutedCandidate = n;
	}

	void CBoxProposal::SetMaxOuputNum(int n)
	{
		m_iMaxOuputNum = n;
	}

	void CBoxProposal::BoxPredictInverse_delta(RECTReal &output, const RECTReal &predict, const RECTReal &src)
	{
		output.left = src.left + predict.left;
		output.top = src.top + predict.top;
		output.right = src.right + predict.right;
		output.bottom = src.bottom + predict.bottom;
	}

	void CBoxProposal::BoxPredictInverse_bbox(RECTReal &output, const RECTReal &predict, const RECTReal &src)
	{
		const float pred_ctr_x = predict.left;
		const float pred_ctr_y = predict.top;
		const float pred_scl_x = predict.right;
		const float pred_scl_y = predict.bottom;

		float src_w = src.right - src.left + 1;
		float src_h = src.bottom - src.top + 1;
		float src_ctr_x = src.left + 0.5f*(src_w - 1);
		float src_ctr_y = src.top + 0.5f*(src_h - 1);

		float output_ctr_x = (pred_ctr_x * src_w) + src_ctr_x;
		float output_ctr_y = (pred_ctr_y * src_h) + src_ctr_y;
		float output_w = exp(pred_scl_x) * src_w;
		float output_h = exp(pred_scl_y) * src_h;


		output.left = output_ctr_x - 0.5f * (output_w - 1);
		output.top = output_ctr_y - 0.5f * (output_h - 1);
		output.right = output_ctr_x + 0.5f * (output_w - 1);
		output.bottom = output_ctr_y + 0.5f * (output_h - 1);
	}
}