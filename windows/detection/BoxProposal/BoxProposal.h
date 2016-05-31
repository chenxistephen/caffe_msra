#ifndef BOX_PROPOSAL_H_
#define BOX_PROPOSAL_H_

#pragma once
#include <Windows.h>
#include <vector>

#include "../../DNNTestLib/CPUData.h"
#include "../../DNNTestLib/DNNTester.h"
#include "../../DNNTestLib/NetUtils.h"
#include "../spm_pool/spm_pool.h"

namespace DNNTestLib
{
	struct RECTReal
	{
		float    left;
		float    top;
		float    right;
		float    bottom;
	};

	struct ProposalBox
	{
		RECTReal box;
		float confidence;
	};

	class CBoxProposal
	{
	public:
		enum BoxTransformType
		{
			BoxTransform_Delta,
			BoxTransform_bbox,
		};

		enum BoxProposalMethod
		{
			BoxProposal_multibox,
			BoxProposal_rpn,
		};

		typedef void (CBoxProposal::*BoxPredictInverseFun)(RECTReal &, const RECTReal &, const RECTReal &);

	public:
		CBoxProposal();
		~CBoxProposal();

		void Release();

		void LoadModel(const std::wstring &file);
		void LoadModel(const std::string &file);
		void LoadModel(std::istream &stream);
		void LoadModel(const char *buffer, int bufferSize, int &bytesRead);

		void GetBoxes(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
			const int iMinBatchSize = 200);

		BoxProposalMethod GetProposalMethod() { return m_boxProposalMethod; }
		
		// BoxPredictInverse function is also used in frcn stage.
		BoxPredictInverseFun get_m_boxPredictInverseFun() { return m_boxPredictInverseFun; }

		void SetMaxEvalutedCandidate(int n);
		void SetMaxOuputNum(int n);

		int m_iMaxEvalutedCandidate = 2000;
		int m_iMaxOuputNum = 500;

	private:		
		void BoxPredictInverse_delta(RECTReal &output, const RECTReal &predict, const RECTReal &src);
		void BoxPredictInverse_bbox(RECTReal &output, const RECTReal &predict, const RECTReal &src);

		typedef void (CBoxProposal::*GetBoxesFun)(std::vector<ProposalBox>&, const CPUData &, const int, const int, const int);
		void GetBoxes_multibox(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
			const int iMinBatchSize = 200);
		void GetBoxes_rpn(std::vector<ProposalBox>& rgBoxes, const CPUData &rspMap, const int iImgWidth, const int iImgHeight,
			const int iMinBatchSize = 200);

		BoxProposalMethod m_boxProposalMethod;
		BoxTransformType m_boxTransformType;

		GetBoxesFun m_getBoxesFun;		
		BoxPredictInverseFun m_boxPredictInverseFun;

		// stand for centroids for multibox, anchors for rpn;
		std::vector<RECTReal> m_centroids;
		DNNTestLib::DNNTester m_net;
		DNNTestLib::CSpmPooler m_sppPooler;

		//// intermediate data
		// rpn
		DNNTestLib::CPUData m_outputData_boxes, m_outputData_confidence;
		// multibox
		DNNTestLib::CPUData m_poolData, m_outputData;
		// stand for srcBoxes for multibox, loacl anchors for rpn;
		std::vector<RECTReal> m_srcBoxes;
	};
}

#endif