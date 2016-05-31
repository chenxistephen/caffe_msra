#ifndef NMS_HPP_
#define NMS_HPP_

#include <vector>
#include <map>
using namespace std;

//struct score {
//	double s;
//	int idx;
//	bool operator() (score i, score j) { return (i.idx < j.idx);}
//} score;

template <typename T1, typename T2>
void nms(int num_boxes, T1 *boxes, T2 *scores, int scores_stride, double overlap, vector<int> &vPick)
{
	vPick.resize(num_boxes);

	//int nSample = (int)mxGetM(input_boxes);
	//int nDim_boxes = (int)mxGetN(input_boxes);

	//T *pBoxes = (T*)mxGetData(input_boxes);

	vector<double> vArea(num_boxes);
	for (int i = 0; i < num_boxes; ++i)
	{
		vArea[i] = double(boxes[2 + i*4] - boxes[0 + i*4] + 1) 
			* (boxes[3 + i*4] - boxes[1 + i*4] + 1);
		if (vArea[i] < 0)
		{
			std::cout << "nms:: Boxes area must >= 0" << std::endl;
			throw std::runtime_error("nms:: Boxes area must >= 0");
		}
	}

	std::multimap<T2, int> vScores;
	for (int i = 0; i < num_boxes; ++i)
		vScores.insert(std::pair<T2,int>(scores[i*scores_stride], i));

	int nPick = 0;

	do 
	{
		int last = vScores.rbegin()->second;
		vPick[nPick] = last;
		nPick += 1;

		for (std::multimap<T2, int>::iterator it = vScores.begin(); it != vScores.end();)
		{
			int it_idx = it->second;
			T1 xx1 = max(boxes[0 + last*4], boxes[0 + it_idx*4]);
			T1 yy1 = max(boxes[1 + last*4], boxes[1 + it_idx*4]);
			T1 xx2 = min(boxes[2 + last*4], boxes[2 + it_idx*4]);
			T1 yy2 = min(boxes[3 + last*4], boxes[3 + it_idx*4]);

			double w = max(0.0, xx2-xx1+1), h = max(0.0, yy2-yy1+1);

			double ov = w*h / (vArea[last] + vArea[it_idx] - w*h);

			if (ov > overlap)
			{
				it = vScores.erase(it);
			}
			else
			{
				it++;
			}
		}

	} while (vScores.size() != 0);

	vPick.resize(nPick);
}

#endif