#include "ImageObjectDetector.h"
#include "nms.hpp"
#include <algorithm>
#include <iostream>

#define ShowTime

#ifdef ShowTime
#include "HighResolutionTimer.h"
#define TimerStart(x) HighResolutionTimer timer_##x; timer_##x.Reset()
#define TimerShow(x) std::cout << #x##" time = " << timer_##x.Duration() << "ms" << endl
#else
#define TimerStart(x) 
#define TimerShow(x)
#endif

CImageObjectDetector::CImageObjectDetector()
	: m_fNMSThreshold(0.3f)
{

}

CImageObjectDetector::~CImageObjectDetector()
{
	Release();
}

void CImageObjectDetector::Release()
{
	m_getObjectInfoCountFun = NULL;
	
	m_netFcs.Clear();
	m_sppPooler.Clear();
	m_boxProposal.Release();

	m_classThresholds.clear();
	m_objInfos.clear();

	m_poolData.Clear();
	m_outputData.Clear();
	m_outputScores.Clear();
	m_outputBoxes.Clear();
}

HRESULT CImageObjectDetector::SetThresholds(const float *pThresholds,
	const int count)
{
	if (count < 0)
	{
		return E_INVALIDARG;
	}

	for (int i = 0; i < min(count, static_cast<int>(m_classThresholds.size())); ++i)
		m_classThresholds[i] = pThresholds[i];

	return S_OK;
}

HRESULT CImageObjectDetector::LoadModel(const char *buffer, int buffer_size, int &bytes_read)
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

	return S_OK;
}

HRESULT CImageObjectDetector::LoadModel(std::istream &stream)
{
	Release();

	m_sppPooler.Load(stream);
	m_netFcs.Load(stream);
	m_boxProposal.LoadModel(stream);

	// TODO: read in thresholds for each class
	int num_classes = 0;
	m_netFcs.GetOutputDim(&num_classes, nullptr, nullptr);

	// Determine getObjectInfoCountFun by boxProposalMethod.
	using DNNTestLib::CBoxProposal;
	CBoxProposal::BoxProposalMethod boxProposalMethod = m_boxProposal.GetProposalMethod();
	if (boxProposalMethod == CBoxProposal::BoxProposal_multibox)
	{
		m_getObjectInfoCountFun = &CImageObjectDetector::GetObjectInfoCount_rcnn;
		m_classThresholds.resize(num_classes, 0);
	}		
	else if (boxProposalMethod == CBoxProposal::BoxProposal_rpn)
	{
		m_getObjectInfoCountFun = &CImageObjectDetector::GetObjectInfoCount_frcn;
		// 0.6 by default for frcn
		m_classThresholds.resize(num_classes - 1, (const float)0.6);
	}		
	else
		DNNTestLib::ThrowException(DNNTestLib::ParameterError, "Undefined proposal method");

	return S_OK;
}

HRESULT CImageObjectDetector::LoadModel(const std::wstring &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

HRESULT CImageObjectDetector::LoadModel(const std::string &file)
{
	std::ifstream fin(file, std::ios::binary | std::ios::in);
	DNNTestLib::ThrowException(DNNTestLib::ModelPackageFileError, "Cannot open model file", fin.is_open());
	LoadModel(fin);
	fin.close();

	return S_OK;
}

bool SortByScore(std::vector<double> &v1, std::vector<double> &v2)
{
	return v1[4]<v2[4];
}

void nms_interclass(std::vector<std::vector<int>> &nms_winners, std::vector<std::vector<DNNTestLib::RECTReal> > &predicted_boxes, DNNTestLib::CPUData &m_outputScores, int num_classes, vector<int> &class_ID)
{
	vector<int>::iterator iter;
	std::vector<std::vector<double>> Cloth_boxes;
	int Cloth_boxes_num = 0;
	for (iter = class_ID.begin(); iter != class_ID.end(); iter++)
	{
		int i = *iter;
		Cloth_boxes_num = Cloth_boxes_num + nms_winners[i].size();
	}
	Cloth_boxes.resize(Cloth_boxes_num);
	float * scoreP = m_outputScores.GetDataPointer();
	Cloth_boxes_num = 0;
	for (iter = class_ID.begin(); iter != class_ID.end(); iter++)
	{
		int i = *iter;
		for (int j = 0; j < nms_winners[i].size(); j++)
		{
			Cloth_boxes[Cloth_boxes_num].resize(7);//[left top right bottom score class_ID nms_winners_index]
			Cloth_boxes[Cloth_boxes_num][0] = (double)predicted_boxes[i][nms_winners[i][j]].left;
			Cloth_boxes[Cloth_boxes_num][1] = (double)predicted_boxes[i][nms_winners[i][j]].top;
			Cloth_boxes[Cloth_boxes_num][2] = (double)predicted_boxes[i][nms_winners[i][j]].right;
			Cloth_boxes[Cloth_boxes_num][3] = (double)predicted_boxes[i][nms_winners[i][j]].bottom;
			Cloth_boxes[Cloth_boxes_num][4] = scoreP[i + 1 + num_classes*nms_winners[i][j]];
			Cloth_boxes[Cloth_boxes_num][5] = i;
			Cloth_boxes[Cloth_boxes_num][6] = nms_winners[i][j];
			Cloth_boxes_num += 1;
		}
	}
	std::sort(Cloth_boxes.begin(), Cloth_boxes.end(), SortByScore);
	
	vector<int> vPick;
	vPick.resize(Cloth_boxes.size());
	vector<double> vArea(Cloth_boxes.size());
	for (int i = 0; i < Cloth_boxes.size(); ++i)
	{
		vArea[i] = double(Cloth_boxes[i][2] - Cloth_boxes[i][0] + 1)
			* (Cloth_boxes[i][3] - Cloth_boxes[i][1] + 1);
	}

	std::multimap<float, int> vScores;
	for (int i = 0; i < Cloth_boxes.size(); ++i)
		vScores.insert(std::pair<float, int>(Cloth_boxes[i][4], i));

	int nPick = 0;

	do
	{
		int last = vScores.rbegin()->second;
		vPick[nPick] = last;
		nPick += 1;

		for (std::multimap<float, int>::iterator it = vScores.begin(); it != vScores.end();)
		{
			int it_idx = it->second;
			float xx1 = max(Cloth_boxes[last][0], Cloth_boxes[it_idx][0]);
			float yy1 = max(Cloth_boxes[last][1], Cloth_boxes[it_idx][1]);
			float xx2 = min(Cloth_boxes[last][2], Cloth_boxes[it_idx][2]);
			float yy2 = min(Cloth_boxes[last][3], Cloth_boxes[it_idx][3]);

			double w = max(0.0, xx2 - xx1 + 1), h = max(0.0, yy2 - yy1 + 1);

			double ov = w*h / (vArea[last] + vArea[it_idx] - w*h);

			if (ov > 0.3)
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
	std::vector<std::vector<double>> Cloth_boxes_tmp;
	Cloth_boxes_tmp.resize(vPick.size());
	for (int i = 0; i < vPick.size(); i++)
	{
		Cloth_boxes_tmp[i] = Cloth_boxes[vPick[i]];
	}
	Cloth_boxes = Cloth_boxes_tmp;

	std::vector<int> tmpIndex;
	for (iter = class_ID.begin(); iter != class_ID.end(); iter++)
	{
		int i = *iter;
		tmpIndex.swap(vector<int>());
		for (int k = 0; k < Cloth_boxes.size(); ++k)
		{
			if (Cloth_boxes[k][5] == i)
			{
				tmpIndex.push_back(Cloth_boxes[k][6]);
			}
		}
		nms_winners[i] = tmpIndex;
	}
}


HRESULT CImageObjectDetector::GetObjectInfoCount(
	const int imgWidth,
	const int imgHeight,
	const DNNTestLib::CPUData &convFeatureMap,
	int& count,
	const float coordinateScale)
{
	return (this->*m_getObjectInfoCountFun)(imgWidth, imgHeight, convFeatureMap, count, coordinateScale);
}


HRESULT CImageObjectDetector::GetObjectInfoCount_frcn(
	const int imgWidth,
	const int imgHeight,
	const DNNTestLib::CPUData &convFeatureMap,
	int& count,
	const float coordinateScale)
{
	HRESULT hr = S_OK;

	// box proposal
	std::vector<DNNTestLib::ProposalBox> boxes;
	TimerStart(boxproposal);
	m_boxProposal.GetBoxes(boxes, convFeatureMap, imgWidth, imgHeight);
	TimerShow(boxproposal);
	if (FAILED(hr)) return hr;

	const int num_rects = (int)boxes.size();
	std::vector<double> rects(4 * num_rects);
	for (int i = 0; i < num_rects; ++i)
	{
		rects[4 * i + 0] = (double)boxes[i].box.left;
		rects[4 * i + 1] = (double)boxes[i].box.top;
		rects[4 * i + 2] = (double)boxes[i].box.right;
		rects[4 * i + 3] = (double)boxes[i].box.bottom;
	}

	int num_classes = 0;
	m_netFcs.GetOutputDim(&num_classes, nullptr, nullptr);

	int fc_input_channels = 0;
	int fc_input_width = 0;
	int fc_input_height = 0;
	m_netFcs.GetInputDim(&fc_input_channels, &fc_input_width, &fc_input_height);

	m_outputScores.Resize(num_rects, num_classes, 1, 1, 0, 0);
	m_outputBoxes.Resize(num_rects, 4 * num_classes, 1, 1, 0, 0);

	TimerStart(spp_fc);
	const int minbatch = 200;
	for (int i = 0; i < static_cast<int>(ceil(num_rects / (float)minbatch)); ++i)
	{
		int sample_idx_start = i * minbatch;
		int sample_num = min(num_rects - sample_idx_start, minbatch);
		
		// spp
		m_sppPooler.SpatialPooling(&m_poolData, &convFeatureMap, sample_num, &rects[0] + 4 * sample_idx_start);
		m_poolData.Resize(sample_num, fc_input_channels, fc_input_width, fc_input_height, 0, 0, DNNTestLib::CPUData::CWHN, false);

		// fc
		DNNTestLib::CPUData outputData_scores, outputData_boxes;
		m_netFcs.TestCPUData(&m_poolData, &outputData_scores, &outputData_boxes);
		assert(outputData_scores.GetPaddingX() == 0 && outputData_scores.GetPaddingY() == 0);
		assert(outputData_boxes.GetPaddingX() == 0 && outputData_boxes.GetPaddingY() == 0);
		memcpy(m_outputScores.GetDataPointer(sample_idx_start, 0, 0, 0), outputData_scores.GetDataPointer(), outputData_scores.GetDataSize() * sizeof(float));
		memcpy(m_outputBoxes.GetDataPointer(sample_idx_start, 0, 0, 0), outputData_boxes.GetDataPointer(), outputData_boxes.GetDataSize() * sizeof(float));
	}
	TimerShow(spp_fc);

	//// get predicted boxes
	
	for (int j = 0; j < num_rects; ++j)
	{
		// scale back, convert to Matlab index.
		boxes[j].box.left = (boxes[j].box.left) * coordinateScale + 1;
		boxes[j].box.top = (boxes[j].box.top) * coordinateScale + 1;
		boxes[j].box.right = (boxes[j].box.right) * coordinateScale + 1;
		boxes[j].box.bottom = (boxes[j].box.bottom) * coordinateScale + 1;
	}

	const float *pOutput_boxes = m_outputBoxes.GetDataPointer();
	std::vector<std::vector<DNNTestLib::RECTReal> > predicted_boxes(num_classes - 1, vector<DNNTestLib::RECTReal>(num_rects));
	DNNTestLib::CBoxProposal::BoxPredictInverseFun fptr = m_boxProposal.get_m_boxPredictInverseFun();
	for (int i = 0; i < num_classes - 1; ++i)
	{		
		for (int j = 0; j < num_rects; ++j)
		{
			predicted_boxes[i][j].left = pOutput_boxes[j * num_classes * 4 + i * 4 + 4];
			predicted_boxes[i][j].top = pOutput_boxes[j * num_classes * 4 + i * 4 + 5];
			predicted_boxes[i][j].right = pOutput_boxes[j * num_classes * 4 + i * 4 + 6];
			predicted_boxes[i][j].bottom = pOutput_boxes[j * num_classes * 4 + i * 4 + 7];


			(m_boxProposal.*fptr)(predicted_boxes[i][j], predicted_boxes[i][j], boxes[j].box);

			// C index.
			predicted_boxes[i][j].left = max(0.f, min(imgWidth * coordinateScale - 1, predicted_boxes[i][j].left - 1));
			predicted_boxes[i][j].top = max(0.f, min(imgHeight * coordinateScale - 1, predicted_boxes[i][j].top - 1));
			predicted_boxes[i][j].right = max(0.f, min(imgWidth * coordinateScale - 1, predicted_boxes[i][j].right - 1));
			predicted_boxes[i][j].bottom = max(0.f, min(imgHeight * coordinateScale - 1, predicted_boxes[i][j].bottom - 1));

			//rects[4 * j + 0] = (double)predicted_boxes[j].left;
			//rects[4 * j + 1] = (double)predicted_boxes[j].top;
			//rects[4 * j + 2] = (double)predicted_boxes[j].right;
			//rects[4 * j + 3] = (double)predicted_boxes[j].bottom;
		}
	}
	

	// nms, pay attention that the output dimension of frcn is actual_class_num + 1.
	TimerStart(nms);
	std::vector<std::vector<int>> nms_winners;
	nms_winners.resize(num_classes - 1);

	for (int i = 0; i < num_classes -  1; ++i)
	{
		for (int j = 0; j < num_rects; ++j)
		{
			rects[4 * j + 0] = (double)predicted_boxes[i][j].left;
			rects[4 * j + 1] = (double)predicted_boxes[i][j].top;
			rects[4 * j + 2] = (double)predicted_boxes[i][j].right;
			rects[4 * j + 3] = (double)predicted_boxes[i][j].bottom;
		}
		nms(num_rects, &rects[0], m_outputScores.GetDataPointer() + i + 1, num_classes, m_fNMSThreshold, nms_winners[i]);
	}

	// NMS between each class write by Jiankangdeng
	//int selected_classes[6] = { 0, 1, 2, 3, 4, 5 };//nms between top wear{0dress,1coatjacket,2sweater,3shirt_long,4shirt,5vest}
	//vector<int> class_id(selected_classes, selected_classes + 6);
	//int selected_classes[5] = { 0, 6, 7, 8, 9};//nms between bottom wear{0dress,6pants,7shorts,8skirt,9skirt_long}
	//vector<int> class_id(selected_classes, selected_classes + 5);
	int selected_classes[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };//nms between top and bottom wear
	vector<int> class_ID(selected_classes, selected_classes + 10);
	nms_interclass(nms_winners, predicted_boxes, m_outputScores, num_classes, class_ID);
	TimerShow(nms);

	// get results, no need to scale back for rpn
	m_objInfos.clear();

	for (int i = 0; i < num_classes - 1; ++i)
	{
		for (int j = 0; j < nms_winners[i].size(); ++j)
		{
			const int idx = nms_winners[i][j];

			const float conf = *(m_outputScores.GetDataPointer() + i + 1 + idx * num_classes);

			if (conf >= m_classThresholds[i])
			{
				const DNNTestLib::RECTReal& box = predicted_boxes[i][idx]/*boxes[idx].box*/;

				IUObjectInfo info;

				info.categoryId = i;
				info.confidence = conf;
				info.position = { box.left, box.top, box.right, box.bottom };

				m_objInfos.push_back(info);
			}
			else
			{
				break;  // nms_winners is sorted w.r.t confidence
			}
		}   // j
	}   // i

	if (!m_objInfos.empty())
	{
		std::sort(m_objInfos.begin(), m_objInfos.end(), [](IUObjectInfo &i, IUObjectInfo &j){return i.confidence > j.confidence; });
	}

	count = (int)m_objInfos.size();

	return S_OK;
}

HRESULT CImageObjectDetector::GetObjectInfoCount_rcnn(
	const int imgWidth,
	const int imgHeight,
	const DNNTestLib::CPUData &convFeatureMap,
    int& count,
	const float coordinateScale)
{ 
	HRESULT hr = S_OK;

	// box proposal
	std::vector<DNNTestLib::ProposalBox> boxes;
	TimerStart(boxproposal);
	m_boxProposal.GetBoxes(boxes, convFeatureMap, imgWidth, imgHeight);
	TimerShow(boxproposal);
	if (FAILED(hr)) return hr;

	const int num_rects = (int)boxes.size();
	std::vector<double> rects(4 * num_rects);
	for (int i = 0; i < num_rects; ++i)
	{
		rects[4 * i + 0] = (double)boxes[i].box.left;
		rects[4 * i + 1] = (double)boxes[i].box.top;
		rects[4 * i + 2] = (double)boxes[i].box.right;
		rects[4 * i + 3] = (double)boxes[i].box.bottom;
	}

    int num_classes = 0;
    m_netFcs.GetOutputDim(&num_classes, nullptr, nullptr);

    int fc_input_channels = 0;
    int fc_input_width = 0;
    int fc_input_height = 0;
    m_netFcs.GetInputDim(&fc_input_channels, &fc_input_width, &fc_input_height);

	m_outputScores.Resize(num_rects, num_classes, 1, 1, 0, 0);

	TimerStart(spp_fc);
	const int minbatch = 200;
	for (int i = 0; i < static_cast<int>(ceil(num_rects / (float)minbatch)); ++i)
	{
		int sample_idx_start = i * minbatch;
		int sample_num = min(num_rects - sample_idx_start, minbatch);

		// spp
		m_sppPooler.SpatialPooling(&m_poolData, &convFeatureMap, sample_num, &rects[0] + 4 * sample_idx_start);
		m_poolData.Resize(sample_num, fc_input_channels, fc_input_width, fc_input_height, 0, 0, DNNTestLib::CPUData::CWHN, false);

		// fc
		m_netFcs.TestCPUData(&m_poolData, &m_outputData);
		assert(m_outputData.GetPaddingX() == 0 && m_outputData.GetPaddingY() == 0);
		memcpy(m_outputScores.GetDataPointer(sample_idx_start, 0, 0, 0), m_outputData.GetDataPointer(), m_outputData.GetDataSize() * sizeof(float));
	}
	TimerShow(spp_fc);

    // nms
	TimerStart(nms);
    std::vector<std::vector<int>> nms_winners;
    nms_winners.resize(num_classes);
    
    for (int i = 0; i < num_classes; ++i)
    {
        nms(num_rects, &rects[0], m_outputScores.GetDataPointer() + i, num_classes, m_fNMSThreshold, nms_winners[i]);
    }
	TimerShow(nms);

    // get results and scale back
    m_objInfos.clear();

    for (int i = 0; i < num_classes; ++i)
    {
        for (int j = 0; j < nms_winners[i].size(); ++j)
        {
            const int idx = nms_winners[i][j];

            const float conf = *(m_outputScores.GetDataPointer() + i + idx * num_classes);

            if (conf >= m_classThresholds[i])
            {    
				const DNNTestLib::RECTReal& box = boxes[idx].box;

				IUObjectInfo info;

                info.categoryId = i;
                info.confidence = conf;
				info.position = { box.left * coordinateScale, box.top * coordinateScale, 
					              box.right * coordinateScale, box.bottom * coordinateScale };

                m_objInfos.push_back(info);
            }
            else
            {
                break;  // nms_winners is sorted w.r.t confidence
            }
        }   // j
    }   // i

    if (!m_objInfos.empty())
    {
		std::sort(m_objInfos.begin(), m_objInfos.end(), [](IUObjectInfo &i, IUObjectInfo &j){return i.confidence > j.confidence; });
    }

    count = (int)m_objInfos.size();

    return S_OK;
}

HRESULT CImageObjectDetector::FillObjectInfos(IUObjectInfo* pObjInfos,
    const int count)
{
    if (m_objInfos.empty())
    {
        return E_FAIL;
    }

    if (pObjInfos == nullptr
        || count > (int)m_objInfos.size())
    {
        return E_INVALIDARG;
    }

	memcpy(pObjInfos, &m_objInfos[0], count * sizeof(IUObjectInfo));

    return S_OK;
}


HRESULT CImageObjectDetector::SetMaxNumPreFilter(const int num)
{
	m_boxProposal.SetMaxEvalutedCandidate(num);

	return S_OK;
}
HRESULT CImageObjectDetector::SetMaxNumPostFilter(const int num)
{
	m_boxProposal.SetMaxOuputNum(num);

	return S_OK;
}


HRESULT CImageObjectDetector::SetNMSThreshold(const float thres)
{
	m_fNMSThreshold = thres;

	return S_OK;
}