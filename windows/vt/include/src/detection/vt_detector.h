#pragma once

#include "vtcore.h"

#include "vt_planarimage.h"

namespace vt
{
	struct CDetectorParams
	{
		CDetectorParams()
			: inputIsTransposed(false)
		{
		}

		bool inputIsTransposed;
	};

	class CDetector
	{
	public:
		CDetector();

		HRESULT InitializeFromFile(const wchar_t* detectorFileName, 
			const CDetectorParams& params);
		HRESULT InitializeFromString(const wchar_t* initString,
			const CDetectorParams& params);

		bool IsInitialized() const;

		HRESULT GetWindowSize(CSize& size) const;
		HRESULT GetNumChannels(unsigned int& numChannels) const;

		template <typename ImageType>
		HRESULT RegisterInputSize(const ImageType& image);
		HRESULT RegisterInputSize(unsigned int numStrideBytes, 
			unsigned int inputHeight, unsigned int inputChannels);
		bool IsInputSizeRegistered(unsigned int numStrideBytes, 
			unsigned int inputHeight, unsigned int inputChannels);
		
		HRESULT SetScoreThreshold(float threshold);

		HRESULT Evaluate(CLumaByteImg& output, 
			const CPlanarByteImg& channels) const;

	private:
		CDetectorParams params;

		bool isInitialized;		

		unsigned int width;
		unsigned int height;
		unsigned int numChannels;
		unsigned int numTrees;

		int scoreThreshold;
		unsigned int numOutScaleShiftBits;

	private:
		typedef char TreeScoreType;

		struct Tree
		{
			unsigned short featureId[3];
			unsigned char thresh[3];
			TreeScoreType score[4];
		};
				
		struct TreeData
		{
			vector<Tree> packed;		
			vector<unsigned char> featureThresholds[3];

			void CheckInvariant(const unsigned int numTrees) const;
		};

		TreeData trees;

		struct FeatureOffsets
		{			
			CSize inputSize;
			struct PackedOffset { int v[3]; };
			vector<PackedOffset> packed;
		};

		vector<FeatureOffsets> featureOffsets;
		const FeatureOffsets* FindFeatureOffsets(
			const CSize& inputSize) const;

		unsigned char ScaleScore(int score) const;

		TreeScoreType EvaluateTree(const Tree& tree,
			const unsigned char* srcPtr, const int (&packedOffsets)[3]) const;

		mutable CLumaShortImg scoreRowBuffer;

		unsigned int EvaluateTreeBatch(short* scoreBuffer, 
			unsigned int scoreBufferLength, const unsigned char* srcPtr, 
			unsigned int srcWidth, unsigned int srcStrideLength, 
			const FeatureOffsets* offsets, unsigned int numBatches) const;
	};

	template <typename ImageType>
	HRESULT CDetector::RegisterInputSize(const ImageType& image)
	{
		return RegisterInputSize(image.StrideBytes(), image.Height(),
			image.Bands());
	}

	struct CDetection
	{
		CRect rect;
		int score;
	};
}
