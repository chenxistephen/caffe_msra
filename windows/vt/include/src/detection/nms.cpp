#include "stdafx.h"

#include "vt_nms.h"

#include <algorithm>

using namespace vt;

namespace vt_internal
{
	// Axis-aligned bounding box for FTD detector result
	struct BoundingBox : public CDetection
	{		
		int id;
		bool suppressed;

		// Constructor with box extents
		BoundingBox(const CDetection& detection)
			: CDetection(detection), id(-1), suppressed(false)  
		{ 
		}

		// Returns the area
		int Area() const
		{
			return rect.Width() * rect.Height();
		}

		// Compute intersection area of this and another Box
		int IntersectionArea(const BoundingBox& box) const
		{
			CRect i;
			i.IntersectRect(&rect, &box.rect);

			return i.Width() * i.Height();
		}

		// Compute union area of this and another Box
		int UnionArea(const BoundingBox& box) const
		{
			return Area() + box.Area() - IntersectionArea(box);
		}

		// Operator used to sort boxes by descending score
		bool operator<(const BoundingBox& box) const
		{
			return score > box.score;
		}

		// Computes overlap with another BoundingBox using the "union" formula
		float UnionOverlap(const BoundingBox& box) const
		{
			return float(IntersectionArea(box)) / UnionArea(box);
		}

		// Computes overlap with another BoundingBox using the "min" formula
		float MinOverlap(const BoundingBox& box) const
		{
			int minArea = (Area() < box.Area()) ? Area() : box.Area();
			return float(IntersectionArea(box)) / minArea;
		}
	};
		
	HRESULT ExtractNonSuppressed(const vector<BoundingBox>& input,
		vector<BoundingBox>& output)
	{
		VT_HR_BEGIN();
		
		output.clear();
		for (size_t k = 0; k < input.size(); k++)
		{
			if (!input[k].suppressed)
			{
				VT_HR_EXIT(output.push_back(input[k]));
			}
		}

		VT_HR_END();
	}

	// Non-maximum supression with "exhaustive-min" algorithm
	HRESULT ExhaustiveMin(vector<BoundingBox> &input,
		vector<BoundingBox> &output, float overlapThresh)
	{
		VT_HR_BEGIN();

		// iterate over all pairs of input boxes
		for (size_t i = 0; i < input.size() - 1; i++)
		{
			for (size_t j = i + 1; j < input.size(); j++)
			{			
				float overlap = input[i].MinOverlap(input[j]);
				if (overlap >= overlapThresh)
				{
					// for every box pair with overlap of (overlap x 100)%
					// "suppress" the box with the lower score
					if (input[i].score > input[j].score)
					{
						input[j].suppressed = true;
					}
					else if (input[j].score > input[i].score)
					{
						input[i].suppressed = true;
					}
				}
			}
		}

		// add all boxes which are not suppressed to the output
		VT_HR_EXIT(ExtractNonSuppressed(input, output));

		// sort output by descending score
		std::sort(output.begin(), output.end());

		VT_HR_END();
	}

	// Non-maximum supression with "exhaustive-union" algorithm
	HRESULT ExhaustiveUnion(vector<BoundingBox> &input,
		vector<BoundingBox> &output, float overlapThresh)
	{
		VT_HR_BEGIN();

		// iterate over all pairs of input boxes
		for (size_t i = 0; i < input.size() - 1; i++)
		{
			for (size_t j = i + 1; j < input.size(); j++)
			{			
				float overlap = input[i].UnionOverlap(input[j]);
				if (overlap >= overlapThresh)
				{
					// for every box pair with overlap of (overlap x 100)%
					// "suppress" the box with the lower score
					if (input[i].score > input[j].score)
					{
						input[j].suppressed = true;
					}
					else if (input[j].score > input[i].score)
					{
						input[i].suppressed = true;
					}
				}
			}
		}

		// add all boxes which are not suppressed to the output
		VT_HR_EXIT(ExtractNonSuppressed(input, output));

		// sort output by descending score
		std::sort(output.begin(), output.end());

		VT_HR_END();
	}

	// Non-maximum supression with "greedy-min" algorithm
	HRESULT GreedyMin(vector<BoundingBox> &input,
		vector<BoundingBox> &output, float overlapThresh)
	{
		VT_HR_BEGIN();

		// sort input by descending score
		std::sort(input.begin(), input.end());
		
		// iterate over all pairs of input boxes
		for (size_t i = 0; i < input.size() - 1; i++)
		{
			// skip any previously suppressed boxes
			if (input[i].suppressed)
				continue;

			for (size_t j = i + 1; j < input.size(); j++)
			{			
				float overlap = input[i].MinOverlap(input[j]);
				if (overlap >= overlapThresh)
				{
					// for every box pair with overlap of (overlap x 100)%
					// "suppress" the box with the lower score
					if (input[i].score > input[j].score)
					{
						input[j].suppressed = true;
					}
					else if (input[j].score > input[i].score)
					{
						input[i].suppressed = true;
					}
				}
			}
		}

		// add all boxes which are not suppressed to the output
		VT_HR_EXIT(ExtractNonSuppressed(input, output));

		VT_HR_END();
	}

	// Non-maximum supression with "greedy-union" algorithm
	HRESULT GreedyUnion(vector<BoundingBox> &input,
		vector<BoundingBox> &output, float overlapThresh)
	{
		VT_HR_BEGIN();

		// sort input by descending score
		std::sort(input.begin(), input.end());

		// iterate over all pairs of input boxes
		for (size_t i = 0; i < input.size() - 1; i++)
		{
			// skip any previously suppressed boxes
			if (input[i].suppressed)
				continue;

			for (size_t j = i + 1; j < input.size(); j++)
			{			
				float overlap = input[i].UnionOverlap(input[j]);
				if (overlap >= overlapThresh)
				{
					// for every box pair with overlap of (overlap x 100)%
					// "suppress" the box with the lower score
					if (input[i].score > input[j].score)
					{
						input[j].suppressed = true;
					}
					else if (input[j].score > input[i].score)
					{
						input[i].suppressed = true;
					}
				}
			}
		}

		// add all boxes which are not suppressed to the output
		VT_HR_EXIT(ExtractNonSuppressed(input, output));

		VT_HR_END();
	}

	HRESULT NonMaxSuppression(vector<BoundingBox>& inputBBs,
		vector<BoundingBox>& outputBBS, float overlapThresh, 
		bool greedyType, bool minDenom)
	{
		VT_HR_BEGIN();

		// Init ids
		for (size_t i = 0; i < inputBBs.size(); i++)
		{
			inputBBs[i].id = (int)i;
		}

		if (greedyType)        
		{
			// Greedy
			if (minDenom)
			{
				VT_HR_EXIT(GreedyMin(inputBBs, outputBBS, overlapThresh));
			}
			else
			{
				VT_HR_EXIT(GreedyUnion(inputBBs, outputBBS, overlapThresh));
			}
		}
		else
		{
			// Exhaustive
			if (minDenom)
			{
				VT_HR_EXIT(ExhaustiveMin(inputBBs, outputBBS, 
					overlapThresh));
			}
			else
			{
				VT_HR_EXIT(ExhaustiveUnion(inputBBs, outputBBS,
					overlapThresh));
			}
		}

		VT_HR_END();
	}
}

using namespace vt_internal;

HRESULT vt::VtNonMaximumSuppression(CDetection* inOutDetections, 
	unsigned int numDetectionsIn, unsigned int* numDetectionsOut, 
	int* outDetectionIndex, float overlapThreshold, 
	NMSMethod method, NMSDenominator denominator)
{
	VT_HR_BEGIN();

	if (numDetectionsIn > 0 && inOutDetections == nullptr)
		return E_INVALIDARG;

	if (overlapThreshold < 0.f || overlapThreshold > 1.f)
		return E_INVALIDARG;
	if (method != NMSMExhaustive && method != NMSMGreedy)
		return E_INVALIDARG;
	if (denominator != NMSDMin && denominator != NMSDUnion)
		return E_INVALIDARG;

	vector<BoundingBox> in;
	for (unsigned int i = 0; i < numDetectionsIn; ++i)
	{
		BoundingBox bb(inOutDetections[i]);
		VT_HR_EXIT(in.push_back(bb));
	}

	vector<BoundingBox> out;
	if (in.size() > 0) // NonMaxSuppression() requires non-empty input
	{
		VT_HR_EXIT(NonMaxSuppression(in, out, overlapThreshold, 
			method == NMSMGreedy, denominator == NMSDMin));
	}

	VT_ASSERT(out.size() <= in.size());

	for (size_t i = 0; i < out.size(); ++i)
	{
		inOutDetections[i].rect.left = out[i].rect.left;
		inOutDetections[i].rect.top = out[i].rect.top;
		inOutDetections[i].rect.right = out[i].rect.left + out[i].rect.Width();
		inOutDetections[i].rect.bottom = out[i].rect.top + out[i].rect.Height();
		inOutDetections[i].score = out[i].score;
	}

	if (numDetectionsOut != nullptr)
		*numDetectionsOut = (unsigned int)out.size();

	if (outDetectionIndex != nullptr)
	{
		for (size_t i = 0; i < out.size(); ++i)
			outDetectionIndex[i] = out[i].id;
	}

	VT_HR_END();
}
