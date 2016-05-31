#pragma once

#include "vtcore.h"
#include "vt_planarimage.h"

namespace vt
{
	struct CDetectionPyramidParams
	{
		CDetectionPyramidParams()
			: minLevelWidth(1),
			minLevelHeight(1),
			maxNumLevels(128),
			maxNumOctaves(32),
			numLevelsPerOctave(4),
			levelSizeMultipleOf(1)
		{
		}

		unsigned int minLevelWidth;
		unsigned int minLevelHeight;
		unsigned int maxNumLevels;
		unsigned int maxNumOctaves;		
		unsigned int numLevelsPerOctave;        
		unsigned int levelSizeMultipleOf;
	};

	struct CDetectionPyramidLevelInfo
	{
		CDetectionPyramidLevelInfo()
			: scale(0), actualScaleX(0), actualScaleY(0),
			actualWidth(0), actualHeight(0)
		{				
		}

		float scale;
		float actualScaleX;
		float actualScaleY;
		unsigned int actualWidth;
		unsigned int actualHeight;
	};

	class CDetectionPyramid
	{
	public:
		CDetectionPyramid();

	private:
		CDetectionPyramid(const CDetectionPyramid&);
		CDetectionPyramid& operator=(const CDetectionPyramid&);

	public:
		HRESULT Initialize(unsigned int baseWidth, unsigned int baseHeight, 
			unsigned int numBands, const CDetectionPyramidParams& params);
		bool IsInitialized() const;
		
		unsigned int NumLevels() const;
		unsigned int NumLevelsPerOctave() const;

		HRESULT Share(CDetectionPyramid& other) const;
		HRESULT ShareLevel(CPlanarByteImg& img, unsigned int level) const;
		HRESULT GetLevelInfo(CDetectionPyramidLevelInfo& levelInfo, 
			unsigned int level) const;

	private:
		HRESULT MakeLevelInfo(unsigned int baseWidth, unsigned int baseHeight, 
			const CDetectionPyramidParams& params);

	private:
		bool isInitialized;
		unsigned int numLevelsPerOctave;

		vector<CDetectionPyramidLevelInfo> linfo;
		vector<CPlanarByteImg> levels;
	};
}