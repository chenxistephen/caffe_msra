#include "stdafx.h"

#include "vt_detector.h"
#include "dhelpers.h"
#include <algorithm>
#include <numeric>

using namespace vt_internal;

namespace vt
{
	void CDetector::TreeData::CheckInvariant(
		const unsigned int numTrees) const 
	{
		VT_ASSERT(packed.size() == numTrees);

		for (size_t i = 0; i < 3; ++i)
			VT_ASSERT(featureThresholds[i].size() == numTrees);

		UNREFERENCED_PARAMETER(numTrees);
	}

	CDetector::CDetector()	
		: isInitialized(false), 		
		width(0),
		height(0),
		numChannels(0),
		numTrees(0),
		scoreThreshold(0),
		numOutScaleShiftBits(0)
	{
	}

	HRESULT CDetector::GetWindowSize(CSize& size) const
	{
		if (!IsInitialized())
			return E_NOINIT;

		size.cx = width;
		size.cy = height;

		return S_OK;
	}

	HRESULT CDetector::GetNumChannels(unsigned int& nChannels) const
	{
		if (!IsInitialized())
			return E_NOINIT;

		nChannels = numChannels;

		return S_OK;			
	}

	HRESULT CDetector::InitializeFromFile(const wchar_t* detectorFileName,
		const CDetectorParams& params)
	{
		FILE* file = nullptr;

		isInitialized = false;

		VT_HR_BEGIN();

		if (_wfopen_s(&file, detectorFileName, L"rt") != 0)
			VT_HR_EXIT(E_NOTFOUND);

		vector<wstring> lines;
		vector<wchar_t> buf;
		VT_HR_EXIT(buf.resize(256));
		while (fgetws(&buf[0], (int)buf.size(), file) != nullptr)
			VT_HR_EXIT(lines.push_back(&buf[0]));

		if (ferror(file))
			VT_HR_EXIT(E_FAIL);

		if (std::any_of(std::begin(lines), std::end(lines), 
			[&](const wstring& l) { return wcslen(l) + 1 == buf.size(); }))
		{
			VT_HR_EXIT(E_BADFORMAT);
		}

		wstring initString;
		for (size_t i = 0; i < lines.size(); ++i)
			VT_HR_EXIT(initString.append(lines[i]));			

		VT_HR_EXIT(InitializeFromString(initString.get_constbuffer(), params));

		VT_HR_EXIT_LABEL();

		if (file != nullptr)
			fclose(file);

		return hr;
	}

	HRESULT CDetector::InitializeFromString(const wchar_t* initString,
		const CDetectorParams& p)
	{
		isInitialized = false;

		params = p;

		trees = TreeData();
		featureOffsets.clear();

		VT_HR_BEGIN();

		vector<wchar_t> buf;
		VT_HR_EXIT(buf.resize(wcslen(initString) + 1));
		if (wcscpy_s(&buf[0], buf.size(), initString) != 0)
			return E_FAIL;

		wchar_t* linePos = &buf[0];
		wchar_t* line = wcstok_s(linePos, L"\n", &linePos);

		if (swscanf_s(line, L"%u %u %u %u %d %u %u", &width, &height, 
			&numChannels, &numTrees, &scoreThreshold, &numOutScaleShiftBits) 
			!= 6)
		{
			return E_BADFORMAT;
		}

		if (width == 0 || height == 0 || numChannels == 0 || numTrees == 0)
			return E_BADFORMAT;

		VT_HR_EXIT(trees.packed.resize(numTrees));

		for (size_t i = 0; i < 3; ++i)
		{
			trees.featureThresholds[i].clear();
		}

		const int maxFeatureId = width * height * numChannels;
		for (size_t i = 0; i < numTrees; ++i)
		{
			line = wcstok_s(linePos, L"\n", &linePos);
			auto& t = trees.packed[i];

			wchar_t* tokenPos = line;
			wchar_t* token = nullptr;

			for (size_t k = 0; k < 3; ++k)
			{
				token = wcstok_s(tokenPos, L" ", &tokenPos);
				if (swscanf_s(token, L"%hu", &t.featureId[k]) != 1 ||
					t.featureId[k] >= maxFeatureId)
				{
					return E_BADFORMAT;
				}

				token = wcstok_s(tokenPos, L" ", &tokenPos);
				int tmp = 0;
				if (swscanf_s(token, L"%d", &tmp) != 1)
					return E_BADFORMAT;

				t.thresh[k] = static_cast<unsigned char>(tmp);
				if (t.thresh[k] != tmp)
					return E_BADFORMAT;

				VT_HR_EXIT(trees.featureThresholds[k].push_back(t.thresh[k]));
			}

			for (size_t k = 0; k < 4; ++k)
			{
				token = wcstok_s(tokenPos, L" ", &tokenPos);
				int tmp = 0;
				if (swscanf_s(token, L"%d", &tmp) != 1)
					return E_BADFORMAT;

				t.score[k] = static_cast<TreeScoreType>(tmp);
				if (t.score[k] != tmp)
					return E_BADFORMAT;
			}
		}

		trees.CheckInvariant(numTrees);

		isInitialized = true;

		VT_HR_END();
	}

	HRESULT CDetector::RegisterInputSize(unsigned int inputStrideWidth, 
		unsigned int inputHeight, unsigned int inputChannels)
	{
		VT_HR_BEGIN();

		if (!IsInitialized())
			return E_NOINIT;

		if (numChannels != inputChannels)
			return E_INVALIDARG;

		if (inputStrideWidth < width || inputHeight < height)
			return E_INVALIDARG;

		// TODO: remove this requirement (right now it is needed for
		// optimizations in tree evaluation)
		if (inputStrideWidth % 16 != 0)
			return E_INVALIDARG;

		CSize inputSize(inputStrideWidth, inputHeight);
		if (FindFeatureOffsets(inputSize) != nullptr)
			return S_OK;

		vt::vector<int> offsetBuffer;
		VT_HR_EXIT(offsetBuffer.resize(inputStrideWidth * inputHeight
			* numChannels));

		unsigned int id = 0;
		for (unsigned int ch = 0; ch < numChannels; ++ch)
		{
			const int imageOffset = ch * inputStrideWidth * inputHeight;

			for (unsigned int y = 0; y < height; ++y)
			{
				for (unsigned int x = 0; x < width; ++x)
				{	
					const int offset = params.inputIsTransposed 
						? y + inputStrideWidth * x
						: x + inputStrideWidth * y;

					VT_ASSERT(id < offsetBuffer.size());
					offsetBuffer[id] = imageOffset + offset;
					++id;				
				}
			}
		}

		FeatureOffsets offsets;
		offsets.inputSize = inputSize;

		VT_HR_EXIT(offsets.packed.resize(numTrees));

		for (size_t k = 0; k < 3; ++k)
		{
			for (size_t i = 0; i < numTrees; ++i)
			{
				offsets.packed[i].v[k] = 
					offsetBuffer[trees.packed[i].featureId[k]];
			}
		}

		VT_HR_EXIT(featureOffsets.push_back(offsets));

		auto maxIt = std::max_element(std::begin(featureOffsets), 
			std::end(featureOffsets),
			[](const FeatureOffsets& a, const FeatureOffsets& b)
		{
			return a.inputSize.cx < b.inputSize.cx;
		});

		VT_HR_EXIT(scoreRowBuffer.Create((*maxIt).inputSize.cx + 16, 1));

		VT_HR_END();
	}

	const CDetector::FeatureOffsets* CDetector::FindFeatureOffsets(
		const CSize& size) const
	{
		auto it = std::find_if(std::begin(featureOffsets), 
			std::end(featureOffsets), [size](const FeatureOffsets& fo)
		{
			return fo.inputSize == size;
		});

		return it != std::end(featureOffsets) ? &(*it) : nullptr;
	}

	bool CDetector::IsInitialized() const
	{
		return isInitialized;
	}

	HRESULT CDetector::SetScoreThreshold(float threshold)
	{
		if (!isInitialized)
			return E_NOINIT;

		scoreThreshold = F2I(threshold * (1 << numOutScaleShiftBits));

		return S_OK;
	}

	CDetector::TreeScoreType CDetector::EvaluateTree(const Tree& tree,
		const unsigned char* srcPtr, const int (&packedOffsets)[3]) const
	{
		if (srcPtr[packedOffsets[0]] >= tree.thresh[0])
		{
			return srcPtr[packedOffsets[2]] >= tree.thresh[2] ? 
				tree.score[3] : tree.score[2];
		}
		else
		{
			return srcPtr[packedOffsets[1]] >= tree.thresh[1] ? 
				tree.score[1] : tree.score[0];
		}
	}

	unsigned char CDetector::ScaleScore(int score) const
	{
		return 128 + static_cast<unsigned char>(Int32Clamper<true, 8>
			::Clamp(score >> numOutScaleShiftBits));
	}

	unsigned int CDetector::EvaluateTreeBatch(short* scoreBuffer, 
		unsigned int scoreBufferLength, const unsigned char* srcPtr, 
		unsigned int srcWidth, unsigned int srcStrideLength, 
		const FeatureOffsets* offsets, unsigned int numBatches) const
	{
		VT_ASSERT(offsets != nullptr);
		VT_ASSERT(scoreBuffer != nullptr);
		VT_ASSERT(scoreBufferLength >= srcWidth + 16);
		VT_ASSERT(srcWidth <= srcStrideLength);
		
		// bail out if we cannot do all 16-byte steps
		if (srcStrideLength % 16 != 0)
		{
			VT_ASSERT(false);
			memset(scoreBuffer, 0, scoreBufferLength * 
				sizeof(scoreBuffer[0]));
			return 0;
		}

		unsigned int numTreesProcessed = 0;

#ifndef _M_ARM
		UNREFERENCED_PARAMETER(numBatches);
		UNREFERENCED_PARAMETER(offsets);
		UNREFERENCED_PARAMETER(srcWidth);
		UNREFERENCED_PARAMETER(srcPtr);
#else
#define NEON_TREE_SETUP_16() \
	const __int64 LUTMaskLo = 3 | (2 << 8) | (3 << 16) | (2 << 24); \
	const __int64 LUTMaskHi = 1 | (1 << 8) | (0 << 16) | (0 << 24); \
	const auto scoreLUTIdx = vcreate_u8(LUTMaskLo | (LUTMaskHi << 32)); \
	const auto bit0mask = vdupq_n_u8(0xfe); \
	const auto bit1mask = vdupq_n_u8(0xfd); \
	const auto bit2mask = vdupq_n_u8(0xfb);			

#define LOAD_NEON_TREE_16(id_, treeIndex) \
	const auto& tr_##id_ = trees.packed[treeIndex]; \
	int32x2_t scores_##id_; scores_##id_.n64_u64[0] = 0; \
	scores_##id_ = vld1_lane_u32(reinterpret_cast<const uint32_t*>( \
	tr_##id_.score), scores_##id_, 0); \
	const auto scoreLUT_##id_ = vtbl1_s8(scores_##id_, scoreLUTIdx); \
	const auto t0_##id_ = vld1q_dup_u8(&tr_##id_.thresh[0]); \
	const auto t1_##id_ = vld1q_dup_u8(&tr_##id_.thresh[1]); \
	const auto t2_##id_ = vld1q_dup_u8(&tr_##id_.thresh[2]); \
	const auto fofs_##id_ = (*offsets).packed[treeIndex].v;

#define EVAL_NEON_TREE_16(id_) \
	const unsigned char* srcPtrX_##id_ = &srcPtr[x]; \
	uint8x16_t v0_##id_, v1_##id_, v2_##id_; \
	uint8x16_t d0_##id_, d1_##id_, d2_##id_, scoreIdx_##id_; \
	v0_##id_ = vld1q_u8(srcPtrX_##id_ + fofs_##id_[0]); \
	d0_##id_ = vbicq_u8(vcltq_u8(v0_##id_, t0_##id_), bit2mask); \
	v1_##id_ = vld1q_u8(srcPtrX_##id_ + fofs_##id_[1]); \
	d1_##id_ = vbicq_u8(vcltq_u8(v1_##id_, t1_##id_), bit1mask); \
	v2_##id_ = vld1q_u8(srcPtrX_##id_ + fofs_##id_[2]); \
	d2_##id_ = vbicq_u8(vcltq_u8(v2_##id_, t2_##id_), bit0mask); \
	scoreIdx_##id_ = vorrq_u8(vorrq_u8(d0_##id_, d1_##id_), d2_##id_); \
	const auto slo_##id_ = vtbl1_s8(scoreLUT_##id_, \
		vget_low_u8(scoreIdx_##id_)); \
	const auto shi_##id_ = vtbl1_s8(scoreLUT_##id_, \
		vget_high_u8(scoreIdx_##id_)); 

#define STORE_NEON_TREE_16() \
	vst1q_s16(&scoreBuffer[x], vaddw_s8(vaddl_s8(slo_0, slo_1), slo_2)); \
	vst1q_s16(&scoreBuffer[x+8], vaddw_s8(vaddl_s8(shi_0, shi_1), shi_2)); 

#define ACC_NEON_TREE_16() \
	const auto buflo = vld1q_s16(&scoreBuffer[x]); \
	const auto bufhi = vld1q_s16(&scoreBuffer[x+8]); \
	const auto sumlo = vaddq_s16(buflo, \
	vaddw_s8(vaddl_s8(slo_0, slo_1), slo_2)); \
	vst1q_s16(&scoreBuffer[x], sumlo); \
	const auto sumhi = vaddq_s16(bufhi, \
	vaddw_s8(vaddl_s8(shi_0, shi_1), shi_2)); \
	vst1q_s16(&scoreBuffer[x+8], sumhi); 

		NEON_TREE_SETUP_16();

		for (unsigned int batchNr = 0; batchNr < numBatches; ++batchNr)
		{
			LOAD_NEON_TREE_16(0, numTreesProcessed);
			LOAD_NEON_TREE_16(1, numTreesProcessed + 1);
			LOAD_NEON_TREE_16(2, numTreesProcessed + 2);
			numTreesProcessed += 3;

			if (batchNr == 0)
			{
				for (unsigned int x = 0; x < srcWidth; x += 16)
				{
					EVAL_NEON_TREE_16(0);
					EVAL_NEON_TREE_16(1);
					EVAL_NEON_TREE_16(2);

					STORE_NEON_TREE_16();
				}		
			}
			else
			{
				for (unsigned int x = 0; x < srcWidth; x += 16)
				{
					EVAL_NEON_TREE_16(0);
					EVAL_NEON_TREE_16(1);
					EVAL_NEON_TREE_16(2);

					ACC_NEON_TREE_16();
				}				
			}
		}
#endif

		if (numTreesProcessed == 0)
		{
			memset(scoreBuffer, 0, scoreBufferLength * 
				sizeof(scoreBuffer[0]));
		}
		
		return numTreesProcessed;
	}

	HRESULT CDetector::Evaluate(vt::CLumaByteImg& output, 
		const vt::CPlanarByteImg& channels) const
	{
		VT_HR_BEGIN();

		if (!IsInitialized())
			return E_NOINIT;

		if (!channels.IsValid())
			return E_INVALIDSRC;

		if (!output.IsValid())
			return E_INVALIDDST;

		if (!AreSameSize(channels, output))
			return E_INVALIDARG;

		CSize size(channels.StrideBytes(), channels.Height());
		const auto* offsets = FindFeatureOffsets(size);
		if (offsets == nullptr)
			return E_INVALIDARG;

		// Map upper left corner of detector window to center pixel, 
		// also accounting for Matlab 1-based image coordinates
		const int dstXOfs = -1 + width / 2;
		const int dstYOfs = -1 + height / 2;

		const int defaultScore = scoreThreshold - 1;

		VT_HR_EXIT(output.Fill(ScaleScore(defaultScore)));

		const int numRows = channels.Height() - height + 1;
		const int numCols = channels.Width() - width + 1;

		for (int y = 0; y < numRows; ++y)
		{
			const unsigned char* srcPtr = channels.Band(0).Ptr(y);

			short* scoreBuf = 
				reinterpret_cast<int16_t*>(scoreRowBuffer.Ptr());
			const unsigned int scoreBufLength = scoreRowBuffer.Width();

			const auto numTreesProcessed = 
				EvaluateTreeBatch(scoreBuf, scoreBufLength, srcPtr,
				numCols, channels.StrideBytes(), offsets, 4);

			const auto* treePtr = &trees.packed[0];
			const auto* ofsPtr = &(*offsets).packed[0];
			unsigned char* dstPtr = output.Ptr(y + dstYOfs) + dstXOfs;

			for (int x = 0; x < numCols; ++x)
			{
				const unsigned char* srcPtrX = &srcPtr[x];

				unsigned int i = numTreesProcessed;
				int score = scoreBuf[x];
				do
				{
					if (score < scoreThreshold)
					{
						score = defaultScore;
						break;
					}

					VT_ASSERT(i < numTrees);
					score += EvaluateTree(treePtr[i], srcPtrX, ofsPtr[i].v);
					++i;					
				}
				while (i < numTrees);

				dstPtr[x] = ScaleScore(score);				
			}
		}

		VT_HR_END();
	}

	bool CDetector::IsInputSizeRegistered(unsigned int numStrideBytes, 
		unsigned int inputHeight, unsigned int inputChannels)
	{
		if (!IsInitialized())
			return false;
		
		CSize size(numStrideBytes, inputHeight);
		
		return FindFeatureOffsets(size) != nullptr 
			&& numChannels == inputChannels;
	}
}