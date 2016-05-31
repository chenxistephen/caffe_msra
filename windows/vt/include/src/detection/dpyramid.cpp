#include "stdafx.h"

#include "vt_dpyramid.h"

namespace vt
{
	CDetectionPyramid::CDetectionPyramid()
		: isInitialized(false),
		numLevelsPerOctave(0)
	{
	}

	HRESULT CDetectionPyramid::Initialize(unsigned int baseWidth, 
		unsigned int baseHeight, unsigned int numBands, 
		const CDetectionPyramidParams& params)
	{
		VT_HR_BEGIN();
				
		isInitialized = false;

		VT_HR_EXIT(MakeLevelInfo(baseWidth, baseHeight, params));

		VT_HR_EXIT(levels.resize(linfo.size()));

		for (size_t i = 0; i < levels.size(); ++i)
		{
			VT_HR_EXIT(levels[i].Create(linfo[i].actualWidth,
				linfo[i].actualHeight, numBands));
		}

		numLevelsPerOctave = params.numLevelsPerOctave;

		isInitialized = true;

		VT_HR_END();	
	}

	bool CDetectionPyramid::IsInitialized() const
	{
		return isInitialized;
	}

	unsigned int CDetectionPyramid::NumLevels() const
	{
		return static_cast<unsigned int>(linfo.size());
	}

	unsigned int CDetectionPyramid::NumLevelsPerOctave() const
	{
		return numLevelsPerOctave;
	}

	HRESULT CDetectionPyramid::GetLevelInfo(CDetectionPyramidLevelInfo& li, 
		unsigned int level) const
	{
		if (level >= levels.size())
			return E_INVALIDARG;

		li = linfo[level];

		return S_OK;
	}
		
	HRESULT CDetectionPyramid::Share(CDetectionPyramid& other) const
	{
		VT_HR_BEGIN();
				
		other.isInitialized = false;

		other.linfo = linfo;
		other.levels.resize(levels.size());
		other.numLevelsPerOctave = numLevelsPerOctave;

		for (size_t i = 0; i < levels.size(); ++i)
		{
			VT_HR_EXIT(levels[i].Share(other.levels[i]));
		}

		other.isInitialized = true;

		VT_HR_END();	
	}

	HRESULT CDetectionPyramid::ShareLevel(CPlanarByteImg& dst, 
		unsigned int level) const
	{
		if (level >= levels.size())
			return E_INVALIDARG;

		return levels[level].Share(dst);
	}

	HRESULT CDetectionPyramid::MakeLevelInfo(unsigned int baseWidth, 
		unsigned int baseHeight, 
		const CDetectionPyramidParams& params) 
	{
		VT_HR_BEGIN();

		linfo.clear();		

		for (unsigned int level = 0; 
			level < params.maxNumLevels && 
			level / params.numLevelsPerOctave < params.maxNumOctaves;
			++level)
		{
			CDetectionPyramidLevelInfo li;

			li.scale = powf(2.0f, -static_cast<float>(level) / 
				params.numLevelsPerOctave);

			li.actualWidth = 
				F2I(baseWidth * li.scale / params.levelSizeMultipleOf) * 
				params.levelSizeMultipleOf;
			if (li.actualWidth < params.minLevelWidth)
				break;

			li.actualHeight = 
				F2I(baseHeight * li.scale / params.levelSizeMultipleOf) * 
				params.levelSizeMultipleOf;
			if (li.actualHeight < params.minLevelHeight)
				break;

			li.actualScaleX = static_cast<float>(li.actualWidth) / 
				baseWidth;
			li.actualScaleY = static_cast<float>(li.actualHeight) / 
				baseHeight;
			
			VT_HR_EXIT(linfo.push_back(li));
		}

		VT_HR_END();
	}
}