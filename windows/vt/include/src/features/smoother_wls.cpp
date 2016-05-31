
#include "StdAfx.h"
#include "vt_smoother_wls.h"
#include "transforms.h"

using namespace vt;

/////////////////////////////////////////////////////////////////
///Compute gaussian weights based on kernel size.
/////////////////////////////////////////////////////////////////
void ComputeSpatialGaussianWeight(vector<float>& weights, int kernelSize, float sigma)
{
    // Compute Gaussian weight for path smoothing

    float sum = 0.0;
    const int num = kernelSize;
    memset(&weights[0], 0, weights.size() * sizeof(float));

    const float alpha = -1.0f / (2.0f * sigma * sigma);
    for (int i = -num; i <= num; i ++)
    {
        float dist = (float)i * (float)i;
        dist = exp(dist * alpha);
        weights[i + num] = dist;
        sum += dist;
    }

    float sum_inv = 1.0f / sum;
    const int totalNum = 2 * num + 1;
    for ( int i = 0; i < totalNum; i ++)
    {
        weights[i] *= sum_inv;
    }
}

/////////////////////////////////////////////////////////////////
/// Declaration of basic path planner class.
/////////////////////////////////////////////////////////////////
class BasicPath
{
public:

    HRESULT Initialize(const vt::vector<vt::CMtx3x3f>& motionSet, int kernelsize, int width, int height);
    HRESULT BasicSmoothing(int iteration, const vt::vector<float>& lamdaList);
    const vt::vector<vt::CMtx3x3f>& GetUpdateSet() const { return m_updateSet; }

protected:
    int		m_iKernelSize;		// radius of smoothing kernel
    int		m_iFrameWidth;		// frame width
    int		m_iFrameHeight;		// frame height	
    int		m_frameNum;			// buffer size of frames used for smoothing

    vt::vector<vt::CMtx3x3f> m_motionSet; // measures the transform between two neighboring frames
    vt::vector<vt::CMtx3x3f> m_updateSet; // warp matrix from input frame to output frame	
    vt::vector<float> m_weightSet;			// weight for path smoothing
    int		m_weightSetSize;

private:
    HRESULT OneUpdate(OUT vt::vector<vt::CMtx3x3f>& ratioUpdateSet, IN const vt::vector<vt::CMtx3x3f>& updateSet, 
        IN const vt::vector<vt::CMtx3x3f>& FList, IN int kernelsize);
};

/////////////////////////////////////////////////////////////////
///Implementation of Class BasicPath
/////////////////////////////////////////////////////////////////

HRESULT BasicPath::Initialize(
    const vector<CMtx3x3f>& motionSet, 
    int kernelsize,
    int width, int height
    )
{
    VT_HR_BEGIN()

    m_motionSet = motionSet;
    m_iKernelSize = kernelsize;
    m_iFrameWidth = width;
    m_iFrameHeight = height;

    m_frameNum = (int) motionSet.size();
    VT_HR_EXIT( m_updateSet.resize(m_frameNum) );
    for(int i = 0; i < m_frameNum; i++)
    {
        m_updateSet[i] = CMtx3x3f().MakeI();
    }

    m_weightSetSize = 2 * m_iKernelSize + 1;	
    VT_HR_EXIT( m_weightSet.resize(m_weightSetSize) );

    const float sigma = sqrt((float) m_iKernelSize);
    ComputeSpatialGaussianWeight(m_weightSet, m_iKernelSize, sigma);

    VT_HR_END()
}

HRESULT BasicPath::OneUpdate(  // one iteration of smoothing path using Gaussian weight
    vector<CMtx3x3f>& ratioUpdateSet, // {\delta}B(t), which will iteratively updates B(t) by B(t) = B(t)*{\delta}B(t)
    const vector<CMtx3x3f>& updateSet, // B(t), mapping original camera pose C(t) to filtered pose P(t) by B(t)*C(t)=P(t)
    const vector<CMtx3x3f>& FList,  // F(t), mapping original camera pose C(t) to C(t-1) by F(t)*C(t)=C(t-1)	
    int kernelsize
    )
{	
    VT_HR_BEGIN()

    vector<CMtx3x3f> HList;	// H(t), mapping filtered camera pose P(t) to P(t-1) by H(t)*P(t)=P(t-1)	
    VT_HR_EXIT( HList.resize(m_frameNum) ); 
    HList[0] = CMtx3x3f().MakeI();
    
    vector<CMtx3x3f> HInvList;	// H(t)^{-1}, which is inverse transform of H(t)
    VT_HR_EXIT( HInvList.resize(m_frameNum) ); 
    HInvList[0] = CMtx3x3f().MakeI();

    CMtx3x3f ivt;
    for(int i = 1; i < m_frameNum; i++)
    {
        HList[i] = updateSet[i].Inv() * FList[i] * updateSet[i - 1];
        HInvList[i] = HList[i].Inv();
    }

    const int windowsize = kernelsize * 2 + 1;
    vector<CMtx3x3f> srcMotionSet;  // from curFrame to the nbFrames
    VT_HR_EXIT( srcMotionSet.resize(windowsize) ); 
    int index0 = kernelsize;
    srcMotionSet[index0] = CMtx3x3f().MakeI();
    
    for (int t = 0; t < m_frameNum; t ++)	//for all the frames
    {	
        CMtx3x3f sumMotion(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        float sumWeight = 0;

        // firstly, from current frame to left neighboring frames
        int iMostLeftIdx = VtMax(t - kernelsize, 0) - t;
        for(int i = -1; i > iMostLeftIdx; i --)
        {
            int ngbIdx = index0 + i;
            srcMotionSet[ngbIdx] = srcMotionSet[ngbIdx + 1] * HList[t + i + 1];
            sumMotion += m_weightSet[ngbIdx] * srcMotionSet[ngbIdx];
            sumWeight += m_weightSet[ngbIdx]; 
        }

        // secondly, from current frame to right neighboring frames
        int iMostRightIdx = VtMin(t + kernelsize, m_frameNum - 1) - t;
        for(int i = 1; i < iMostRightIdx; i ++)
        {
            int ngbIdx = index0 + i;
            srcMotionSet[ngbIdx] = srcMotionSet[ngbIdx - 1] * HInvList[t + i];
            sumMotion += m_weightSet[ngbIdx] * srcMotionSet[ngbIdx];
            sumWeight += m_weightSet[ngbIdx]; 
        }

        sumMotion += m_weightSet[index0] * srcMotionSet[index0];
        sumWeight += m_weightSet[index0]; 

        ratioUpdateSet[t] = sumMotion / sumWeight;
    }

    VT_HR_END()
}

HRESULT BasicPath::BasicSmoothing(int iteration, const vector<float>& lamdaList)
{
    VT_HR_BEGIN()

    CMtx3x3f matE = CMtx3x3f().MakeI(); // identity matrix

    vector<CMtx3x3f> ratioUpdateSet;
    VT_HR_EXIT( ratioUpdateSet.resize(m_frameNum) );

    for(int i = 0; i < iteration; i++)
    {
        VT_HR_EXIT( OneUpdate(ratioUpdateSet, m_updateSet, m_motionSet, m_iKernelSize) );

        if(lamdaList.size() == 0)
        {
            for(int t = 0; t < m_frameNum; t++)
            {
                // Update B(t) by B(t) = B(t)*{\delta}B(t)
                m_updateSet[t] = m_updateSet[t] * ratioUpdateSet[t];
            }
        }
        else
        {		
            for(int t = 0; t < m_frameNum; t++)
            {   
                // Update B(t) by B(t) = B(t)*{\delta}B(t)*(1-lamda) + I*lamda
                // Here, the optimization requires result partially close to original frame
                CMtx3x3f matNewUpdate = m_updateSet[t] * ratioUpdateSet[t];
                m_updateSet[t] = matNewUpdate + lamdaList[t] * (matE - matNewUpdate);
            }
        }	
    }

    VT_HR_END()
}

/////////////////////////////////////////////////////////////////
/// Declaration of active path planner class.
/////////////////////////////////////////////////////////////////
class AdaptivePath
{
public:
    HRESULT Initialize(const vt::vector<vt::CMtx3x3f>& motionSet, int kernelsize, int width, int height);
    HRESULT InitialSmoothingLoop(int iteration);
    HRESULT IterativeSmoothingLoop(const vt::vector<float>& lamdaList, int iteration);
    const vt::vector<vt::CMtx3x3f>& GetUpdateSet() const { return m_updateSet; }
    void SetSmoothingContext(int bilateral_iteration, float sigma_result2guidedpath, float sigma_guided2original);
private:
    HRESULT ComputeAdaptiveWeight();
    HRESULT GuidedPathSmoothing(IN OUT vt::vector<float>& guidedPath);

    HRESULT OneUpdate(OUT vt::vector<vt::CMtx3x3f>& ratioUpdateSet, IN const vt::vector<vt::CMtx3x3f>& updateSet, 
        IN const vt::vector<vt::CMtx3x3f>& FList, IN const vt::vector<vt::vector<float>>& weights, IN int kernelsize,
        IN OUT vt::vector<vt::CMtx3x3f>& HList, IN OUT vt::vector<vt::CMtx3x3f>& HInvList, IN OUT vt::vector<vt::CMtx3x3f>& srcMotionSet);
        
private:
    int		m_iKernelSize;		// radius of smoothing kernel
    int		m_iFrameWidth;		// frame width
    int		m_iFrameHeight;		// frame height	
    int		m_frameNum;			// buffer size of frames used for smoothing

    vt::vector<vt::CMtx3x3f> m_motionSet; // measures the transform between two neighboring frames
    vt::vector<vt::CMtx3x3f> m_updateSet; // warp matrix from input frame to output frame	

    vt::vector<float> m_spatialWeight; // Gaussian weights for temporal neighborhoods
    vt::vector<vt::vector<float>> m_bilateralWeights; // Bilateral weights from guided path
    vt::vector<vt::CVec2f> m_CList;  // translation components used for camera pose representation

    float m_sigma_Result2GuidedPath;  // value = 20.  smaller, result path more similar to guided path
    float m_sigma_Guided2OriginalPath; // value = 10. smaller, result more similar to original path
    int	   m_bilateral_iteration;		// value = 20. smaller, result more similar to original path

    // temporary buffers to avoid repetitively memory allocation in core loops.
    // these buffers will be created on the initialize function
    vt::vector<vt::CMtx3x3f> m_ratioUpdateSet;  // same length as input montin set
    vt::vector<vt::CMtx3x3f> m_HList;			// same length as input montin set
    vt::vector<vt::CMtx3x3f> m_HInvList;		// same length as input montin set
    vt::vector<vt::CMtx3x3f> m_srcMotionSet;	// 2 * kernel_size + 1
};

/////////////////////////////////////////////////////////////////
///Implementation of Class AdaptivePath
/////////////////////////////////////////////////////////////////
HRESULT AdaptivePath::Initialize(
    const vector<CMtx3x3f>& motionSet,
    int kernelsize, 
    int width, 
    int height
    )
{
    VT_HR_BEGIN()

    m_motionSet = motionSet;
    m_iKernelSize = kernelsize;
    m_iFrameWidth = width;
    m_iFrameHeight = height;

    m_frameNum = (int) motionSet.size();
    VT_HR_EXIT( m_updateSet.resize(m_frameNum) );
    for(int i = 0; i < m_frameNum; i++)
    {
        m_updateSet[i] = CMtx3x3f().MakeI();
    }	

    const int windowSize = 2 * m_iKernelSize + 1;

#ifdef SIGMA_MODIFIED
    const float sigma = VtMax(m_iKernelSize / 1.5f, 0.1f);
    const float sigma_result2guidedpath = 10; 
    const float sigma_guided2original = 20; 
    const int bilateral_iteration = 5;
#else
    const float sigma = sqrt((float) m_iKernelSize); 
    const float sigma_result2guidedpath = 20; 
    const float sigma_guided2original = 10; 
    const int bilateral_iteration = 20;
#endif

    // Compute spatial Gaussian weights, only related to relationship between frames indices
    VT_HR_EXIT( m_spatialWeight.resize(windowSize) );
    ComputeSpatialGaussianWeight(m_spatialWeight, m_iKernelSize, sigma);

    // Compute bilateral range weights from the difference in guided camera path
    VT_HR_EXIT( m_bilateralWeights.resize(m_frameNum) );
    for(int i = 0; i < m_frameNum; i++)
    {
        VT_HR_EXIT( m_bilateralWeights[i].resize(windowSize) );
        memset(&m_bilateralWeights[i][0], 0, m_bilateralWeights[i].size() * sizeof(float));
    }

    VT_HR_EXIT( m_CList.resize(m_frameNum) );
    m_CList[0].x = m_motionSet[0](0,2);
    m_CList[0].y = m_motionSet[0](1,2);
    for(int i = 1; i < m_frameNum; i++)
    {
        m_CList[i].x = m_CList[i-1].x + m_motionSet[i](0,2);
        m_CList[i].y = m_CList[i-1].y + m_motionSet[i](1,2);
    }

    SetSmoothingContext(bilateral_iteration, sigma_result2guidedpath, sigma_guided2original);

    VT_HR_EXIT( ComputeAdaptiveWeight() );

    // Initialize temporary buffers.
    VT_HR_EXIT( m_ratioUpdateSet.resize(m_frameNum) ); 
    VT_HR_EXIT( m_HList.resize(m_frameNum) );		
    VT_HR_EXIT( m_HInvList.resize(m_frameNum) );
    VT_HR_EXIT( m_srcMotionSet.resize(m_iKernelSize * 2 + 1) );

    VT_HR_END()
}

void AdaptivePath::SetSmoothingContext(
    int bilateral_iteration, 
    float sigma_result2guidedpath,
    float sigma_guided2original
    )
{
    m_sigma_Result2GuidedPath = sigma_result2guidedpath;  
    m_sigma_Guided2OriginalPath = sigma_guided2original; 
    m_bilateral_iteration = bilateral_iteration;
}

HRESULT AdaptivePath::InitialSmoothingLoop(int iteration)
{
    VT_HR_BEGIN()

    for(int i = 0; i < iteration; i++)
    {
        VT_HR_EXIT( OneUpdate(m_ratioUpdateSet, m_updateSet, m_motionSet, m_bilateralWeights, m_iKernelSize, m_HList, m_HInvList, m_srcMotionSet) );

        for(size_t t = 0; t < m_updateSet.size(); t++)
        {  // Update B(t) by B(t) = B(t)*{\delta}B(t)
            m_updateSet[t] = m_updateSet[t] * m_ratioUpdateSet[t];
        }
    }

    VT_HR_END()
}

HRESULT AdaptivePath::IterativeSmoothingLoop(
    const vector<float>& lamdaList, // {\lamda}(t) for tradeoff between original path and smoothed path
    int iteration
    )
{
    VT_HR_BEGIN()

    for(size_t i = 0; i < m_updateSet.size(); i++)
    {
        m_updateSet[i] = CMtx3x3f().MakeI();
    }

    CMtx3x3f matE = CMtx3x3f().MakeI(); // identity matrix

    for(int i = 0; i < iteration; i++)
    {
        VT_HR_EXIT( OneUpdate(m_ratioUpdateSet, m_updateSet, m_motionSet, m_bilateralWeights, m_iKernelSize, m_HList, m_HInvList, m_srcMotionSet) );

        for(size_t t = 0; t < m_updateSet.size(); t++)
        {   
            // Update B(t) by B(t) = B(t)*{\delta}B(t)*(1-lamda) + I*lamda
            // Here, the optimization requires result partially close to original frame
            CMtx3x3f matNewUpdate = m_updateSet[t] * m_ratioUpdateSet[t];
            m_updateSet[t] = matNewUpdate + lamdaList[t] * (matE - matNewUpdate);
        }	
    }

    VT_HR_END()
}

HRESULT AdaptivePath::OneUpdate( // one iteration of smoothing path using Adaptive weight
    vector<CMtx3x3f>& ratioUpdateSet, // {\delta}B(t), which will iteratively updates B(t) by B(t) = B(t)*{\delta}B(t)
    const vector<CMtx3x3f>& updateSet, // B(t), mapping original camera pose C(t) to filtered pose P(t) by B(t)*C(t)=P(t)
    const vector<CMtx3x3f>& FList,  // F(t), mapping original camera pose C(t) to C(t-1) by F(t)*C(t)=C(t-1)
    const vector<vector<float>>& weights, // bilateral weight W(t, i) for each frame t and its neigbors i
    int kernelsize,
    vector<CMtx3x3f>& HList, // H(t), mapping filtered camera pose P(t) to P(t-1) by H(t)*P(t)=P(t-1)
    vector<CMtx3x3f>& HInvList,  // H(t)^{-1}, which is inverse transform of H(t)
    vector<CMtx3x3f>& srcMotionSet // from curFrame to the nbFrames
    )
{
    VT_HR_BEGIN()

    VT_ASSERT(HList.size() == m_frameNum);
    VT_ASSERT(HInvList.size() == m_frameNum);
    VT_ASSERT(srcMotionSet.size() == (kernelsize * 2 + 1));

    HList[0] = CMtx3x3f().MakeI();
    HInvList[0] = CMtx3x3f().MakeI();
    CMtx3x3f ivt;

    for(size_t i = 1; i < updateSet.size(); i++)
    {
        HList[i] = updateSet[i].Inv() * FList[i] * updateSet[i - 1];
        HInvList[i] = HList[i].Inv();
    }

    const int windowsize = kernelsize * 2 + 1;
    int index0 = kernelsize;
    srcMotionSet[index0] = CMtx3x3f().MakeI();

    for (int t = 0; t < m_frameNum; t ++)	//for all the frames
    {	
        CMtx3x3f sumMotion(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
        float sumWeight = 0;

        // firstly, from current frame to left neighboring frames
        int iMostLeftIdx = VtMax(t - kernelsize, 0) - t;
        for(int i = -1; i > iMostLeftIdx; i --)
        {
            int ngbIdx = index0 + i;
            srcMotionSet[ngbIdx] = srcMotionSet[ngbIdx + 1] * HList[t + i + 1];
            float weight = weights[t][ngbIdx];
            sumMotion += weight * srcMotionSet[ngbIdx];
            sumWeight += weight; 
        }

        // secondly, from current frame to right neighboring frames
        int iMostRightIdx = VtMin(t + kernelsize, m_frameNum - 1) - t;
        for(int i = 1; i < iMostRightIdx; i ++)
        {
            int ngbIdx = index0 + i;
            srcMotionSet[ngbIdx] = srcMotionSet[ngbIdx - 1] * HInvList[t + i];
            float weight = weights[t][ngbIdx];
            sumMotion += weight * srcMotionSet[ngbIdx];
            sumWeight += weight; 
        }

        sumMotion += weights[t][index0] * srcMotionSet[index0];
        sumWeight += weights[t][index0]; 

        ratioUpdateSet[t] = sumMotion * (1.0f / sumWeight);
    }

    VT_HR_END()
}

HRESULT AdaptivePath::ComputeAdaptiveWeight()
{
    VT_HR_BEGIN()

    const int frameNum = (int) m_CList.size();
    vector<float> srcX;  
    VT_HR_EXIT( srcX.resize(frameNum) );
    vector<float> srcY;  
    VT_HR_EXIT( srcY.resize(frameNum) );
    for(int i = 0; i < frameNum; i++)
    {
        srcX[i] = m_CList[i].x;
        srcY[i] = m_CList[i].y;
    }

    VT_HR_EXIT( GuidedPathSmoothing(srcX) );
    VT_HR_EXIT( GuidedPathSmoothing(srcY) );

    vector<CVec2f> filteredPath;
    VT_HR_EXIT( filteredPath.resize(frameNum) );
    for(int i = 0; i < frameNum; i++)
    {
        filteredPath[i].x = srcX[i];
        filteredPath[i].y = srcY[i];
    }

    const int radius = m_iKernelSize;
    const int windowSize = 2 * radius + 1;
    const float camsigma = m_sigma_Result2GuidedPath;

    //compute camera distance weight
    const float alpha = -1.0f / (2.0f * camsigma * camsigma);
    for(int t = 0; t < frameNum; t++)
    {
        float weightSum = 0.0f;
        for(int i = -radius; i <= radius; i ++)
        {
            // compute range weight
            int curNbIdx = t + i;		
            float rangeWeight = 0.0f;
            if(curNbIdx >= 0 && curNbIdx < frameNum)
            {
                // compute camera distance according to translations
                float distX = filteredPath[curNbIdx].x - filteredPath[t].x;
                distX = distX > 0 ? distX : -distX;
                float distY = filteredPath[curNbIdx].y - filteredPath[t].y;
                distY = distY > 0 ? distY : -distY;
                float camdis = (distX + distY) * 100.0f;
                rangeWeight = exp(camdis * camdis * alpha);  
            }

            // combine bilateral weight from range and spatial Gaussian weight
            int idx = i + radius;
            float jointWeight = rangeWeight * m_spatialWeight[idx];
            m_bilateralWeights[t][idx] = jointWeight;
            weightSum += jointWeight;
        }

        // bilateral weight normalization
        float sumInv = 1.0f / weightSum;
        for(int i = 0; i < windowSize; i++)
        {
            m_bilateralWeights[t][i] *= sumInv;
        }
    }

    VT_HR_END()
}

HRESULT AdaptivePath::GuidedPathSmoothing(vector<float>& guidedPath)
{
    VT_HR_BEGIN()

    const int frameNum = (int) guidedPath.size();
    const int radius = m_iKernelSize;
    const int windowSize = 2 * radius + 1;

    vector<vector<float>> weights;
    VT_HR_EXIT( weights.resize(frameNum) );
    for(int t = 0; t < frameNum; t ++)
    {
        VT_HR_EXIT( weights[t].resize(windowSize) );
    }

    float camsigma = m_sigma_Guided2OriginalPath;
    const float alpha = -(100.0f * 100.0f) / (2 * camsigma * camsigma);
    for(int t = 0; t < frameNum; t++)
    {
        float weightSum = 0.0f;
        for(int i = -radius; i <= radius; i ++)
        {
            // compute range weight
            int curNbIdx = t + i;		
            float rangeWeight = 0.0f;
            if(curNbIdx >= 0 && curNbIdx < frameNum)
            {
                float camdis = guidedPath[curNbIdx] - guidedPath[t];
                rangeWeight = exp(camdis * camdis * alpha);  
            }

            // combine bilateral weight from range and spatial Gaussian weight
            int idx = i + radius;
            float jointWeight = rangeWeight * m_spatialWeight[idx];
            weights[t][idx] = jointWeight;
            weightSum += jointWeight;
        }

        // bilateral weight normalization
        float sumInv = 1.0f / weightSum;
        for(int i = 0; i < windowSize; i++)
        {
            weights[t][i] *= sumInv;
        }
    }

    // Iteratively smooth guided Path
    vector<float> guidedPath_backup;
    VT_HR_EXIT( guidedPath_backup.resize(frameNum) );

    for(int Iter = 0; Iter < m_bilateral_iteration; Iter ++)
    {
        memcpy(&guidedPath_backup[0], &guidedPath[0], sizeof(float) * frameNum);
        memset(&guidedPath[0], 0, sizeof(float) * frameNum);

        for(int i = 0; i < frameNum; i++)
        {
            for(int j = -radius; j <= radius; j ++)
            {
                int NgbIdx = i + j;
                if(NgbIdx > 0 && NgbIdx < frameNum)
                {
                    guidedPath[i] += weights[i][j + radius] * guidedPath_backup[NgbIdx];   
                }
            }
        }
    }

    VT_HR_END()
}

/////////////////////////////////////////////////////////////////
/// Declaration of weighter updator class.
/////////////////////////////////////////////////////////////////
class LamdaUpdator
{
public:
    struct node
    {
        int index;
        float crop;
        float wobble;
        float lamda;
        node(int index_,float crop_,float wobble_,float lamda_)
        {
            this->index = index_;
            this->crop = crop_;
            this->wobble = wobble_;
            this->lamda = lamda_;
        }
    };

private:
    int m_iFrameWidth;		// frame width
    int m_iFrameHeight;		// frame height	
    int m_iKernelSize;		// radius of smoothing kernel

    vt::vector<node> m_unsatisfied_frames;   // information of unsatisfied frames (large crop, wobble)
    vt::vector<float> m_lamdaList;          // trade-off between origianl and smoothed, {\lamda}(t)   

    float m_fCropThres;		// crop ratio constraint
    float m_fWobbleThres;      // wobble constraint
    float m_fLamdaStep;		// increment step of lamda update at iterations
private:
    float CalcWobbleRatio(const vt::CMtx3x3f &update);
    float CalcCropSizeRatio(const vt::CMtx3x3f &update, int width, int height);
    HRESULT SmoothingLamdas(vt::vector<float>& lamdaList, int kernelsize, int iteration);

public:
    LamdaUpdator(){};
    ~LamdaUpdator(){};

    HRESULT Initialize(const vt::vector<float> &lamdaList, int width, int height, int kernelsize, float fCrop = 0.9, float fWobble = 1.04, float fLamdaStep = 0.1);
    HRESULT InitUpdateLamdaList(const vt::vector<vt::CMtx3x3f>& updateSet, vt::vector<float>& lamdaList);
    HRESULT UpdateLamdaList(const vt::vector<vt::CMtx3x3f>& updateSet, vt::vector<float>& lamdaList);
    HRESULT UpdateLamdaByCropHardConstraint(const vt::vector<vt::CMtx3x3f>& updateSet, vt::vector<float>& lamdaList);
};

/////////////////////////////////////////////////////////////////
///Implementation of LamdaUpdator
/////////////////////////////////////////////////////////////////
HRESULT LamdaUpdator::Initialize(
    const vector<float>& lamdaList,
    int width, int height,
    int kernelsize,
    float fCrop, 
    float fWobble,
    float fLamdaStep
    )
{
    VT_HR_BEGIN()

    m_iFrameWidth = width;
    m_iFrameHeight = height;
    m_iKernelSize = kernelsize;
    VT_HR_EXIT( m_unsatisfied_frames.reserve(lamdaList.size()) );
    
    m_lamdaList = lamdaList;

    m_fCropThres = fCrop;
    m_fWobbleThres = fWobble;
    m_fLamdaStep = fLamdaStep;

    VT_HR_END()
}

HRESULT LamdaUpdator::InitUpdateLamdaList(const vector<CMtx3x3f>& updateSet, vector<float>& lamdaList)
{
    VT_HR_BEGIN()

    VT_ASSERT(m_lamdaList.size() == updateSet.size());
    VT_ASSERT(lamdaList.size() == updateSet.size());

    // clear unsatisfied frame lists
    m_unsatisfied_frames.clear();

    for(size_t i = 0; i < updateSet.size(); i++)
    {
        // compute crop ratio and wobble for current warping matrix UpdateSet (B(t))
        float crop_ratio = CalcCropSizeRatio(updateSet[i], m_iFrameWidth, m_iFrameHeight);
        float wobble_ratio = CalcWobbleRatio(updateSet[i]);

        // when crop ratio or wobble beyond the given threshold
        if(crop_ratio < m_fCropThres || wobble_ratio > m_fWobbleThres)
        {
            // mapping crop ratio to lamda
            float lamdaCrop = VtMax(0.0f, m_fCropThres - crop_ratio);
            lamdaCrop /= VtMax(0.001f, 1.0f - crop_ratio); 

            // mapping wobble to lamda
            float lamdaWobble = VtMax(0.0f, wobble_ratio - m_fWobbleThres);
            lamdaWobble /= VtMax(0.001f,  wobble_ratio - 1.0f); 
            
            // combine both to final lamda value for current frame
            float lamda = 0.5f * (lamdaCrop + lamdaWobble * 2.5f);
            m_lamdaList[i] = VtClamp(lamda, 0.0f, 1.0f);

            // push unsatisfied frames (large crop size, wobble distortion) to the list:
            VT_HR_EXIT( m_unsatisfied_frames.push_back(node(static_cast<int>(i), crop_ratio, wobble_ratio, m_lamdaList[i])) );
        }
    }

    if(m_unsatisfied_frames.size() == 0)
    {
        hr = S_FALSE;
    }
    else
    {
        VT_HR_EXIT( SmoothingLamdas(m_lamdaList, m_iKernelSize, 3) );
        for(size_t i = 0; i < m_unsatisfied_frames.size(); i++)
        {
            int index = m_unsatisfied_frames[i].index;
            m_unsatisfied_frames[i].lamda = m_lamdaList[index];
        }

        for(size_t i = 0; i < m_lamdaList.size(); i ++)
        {
            lamdaList[i] = m_lamdaList[i];
        }
    }

    VT_HR_END()
}

HRESULT LamdaUpdator::UpdateLamdaList(const vector<CMtx3x3f>& updateSet, vector<float>& lamdaList)
{
    VT_HR_BEGIN()

    // clear unsatisfied frame lists
    m_unsatisfied_frames.clear();

    for(size_t i = 0; i < updateSet.size(); i++)
    {
        // compute crop ratio and wobble for current warping matrix UpdateSet (B(t))
        float crop_ratio = CalcCropSizeRatio(updateSet[i], m_iFrameWidth, m_iFrameHeight);
        float wobble_ratio = CalcWobbleRatio(updateSet[i]);

        // when crop ratio or wobble beyond the given threshold
        if(crop_ratio < m_fCropThres || wobble_ratio > m_fWobbleThres)
        {
            //get the lamda value from the last iteration
            float lamdaOld = m_lamdaList[i];

            // mapping crop ratio to lamda
            float lamdaCrop = VtMax(0.0f, m_fCropThres - crop_ratio);
            lamdaCrop /= VtMax(0.001f, 1.0f - crop_ratio); 

            // mapping wobble to lamda
            float lamdaWobble = VtMax(0.0f, wobble_ratio - m_fWobbleThres);
            lamdaWobble /= VtMax(0.001f,  wobble_ratio - 1.0f); 

            // combine both to new lamda value for current frame
            float lamdaNew = VtMax(lamdaCrop, lamdaWobble * 2.5f);
            lamdaNew = VtClamp(lamdaNew, 0.0f, 1.0f);

            // compute the increment of lamda value after comparing its old and new values
            float step = m_fLamdaStep;
            if(lamdaNew > lamdaOld){
                step = VtMax(0.5f * (lamdaNew - lamdaOld), m_fLamdaStep);
            }

            // increase lamda to enforce the result close to original
            m_lamdaList[i] = VtMin(lamdaOld + step, 1.0f);

            // push unsatisfied frames (large crop size, wobble distortion) to the list:
            VT_HR_EXIT( m_unsatisfied_frames.push_back(node(static_cast<int>(i), crop_ratio, wobble_ratio, m_lamdaList[i])) );
        }
    }

    if(m_unsatisfied_frames.size() == 0)
    {
        hr = S_FALSE;
    }
    else
    {
        for(size_t i = 0; i < m_unsatisfied_frames.size(); i++)
        {
            int index = m_unsatisfied_frames[i].index;
            m_unsatisfied_frames[i].lamda = m_lamdaList[index];
        }

        for(size_t i = 0; i < m_lamdaList.size(); i ++)
        {
            lamdaList[i] = m_lamdaList[i];
        }
    }

    VT_HR_END()
}

HRESULT LamdaUpdator::UpdateLamdaByCropHardConstraint(const vector<CMtx3x3f>& updateSet, vector<float>& lamdaList)
{
    VT_HR_BEGIN()

    // enforce crop hard constraint for post-smoothing
    for(size_t i = 0; i < lamdaList.size(); i ++)
    {
        float crop_ratio = CalcCropSizeRatio(updateSet[i], m_iFrameWidth, m_iFrameHeight); 
        if(crop_ratio < m_fCropThres)
        {
            // set lamda to make crop size within given crop size threshold
            lamdaList[i] = VtMin(1.0f, lamdaList[i] + 4.0f * (m_fCropThres - crop_ratio));
        }
        else
        {
            // let lamda be when crop size is reasonable
            lamdaList[i] = 0.0f;
        }
    }

    VT_HR_EXIT( SmoothingLamdas(lamdaList, m_iKernelSize, 1) );

    VT_HR_END()
}

float LamdaUpdator::CalcCropSizeRatio(
    const CMtx3x3f &update, 
    int width, 
    int height
    )
{
    CMtx3x3f mtxH;
    CMtx3x3f mtxN;	mtxN.MakeI();
    float maxDim = (float)VtMax(width, height);
    mtxN(0,0) = 1.0f / maxDim;
    mtxN(1,1) = 1.0f / maxDim;

    // Get warping matrix to compute frame crop size
    mtxH = mtxN.Inv() * update * mtxN;

    // Get four coordinates of frame corners
    CVec3f src00 = CVec3f(0.0f, 0.0f, 1.0f);
    CVec3f src01 = CVec3f((float)width, 0.0f, 1.0f);
    CVec3f src10 = CVec3f(0.0f, (float)height, 1.0f);
    CVec3f src11 = CVec3f((float)width, (float)height, 1.0f);

    // Get warped coordinates of frame corners
    CVec3f dst00 = mtxH * src00; dst00.x /= dst00.z; dst00.y /= dst00.z;
    CVec3f dst01 = mtxH * src01; dst01.x /= dst01.z; dst01.y /= dst01.z;
    CVec3f dst10 = mtxH * src10; dst10.x /= dst10.z; dst10.y /= dst10.z;
    CVec3f dst11 = mtxH * src11; dst11.x /= dst11.z; dst11.y /= dst11.z;

    float ptL = VtMax(dst00.x, dst10.x);
    float ptR = VtMin(dst01.x, dst11.x);
    float ptT = VtMax(dst00.y, dst01.y);
    float ptB = VtMin(dst10.y, dst11.y);

    // Get valid rectangle after cropping, crop center is original frame center
    float rectW = VtMin(width * 0.5f - ptL, ptR - width * 0.5f);
    float rectH = VtMin(height * 0.5f - ptT, ptB - height * 0.5f);

    // Compute crop size ratio
    float crop_ratio = VtMin(rectW * 2.0f / width, rectH * 2.0f / height);
    crop_ratio = VtClamp(crop_ratio, 0.0f, 1.0f);
    return crop_ratio;	
}

float LamdaUpdator::CalcWobbleRatio(const CMtx3x3f &update)
{
    CSolveSVD<float> svd(update);
    float s1 = svd.D()(0,0);
    float s2 = svd.D()(1,1);

    float wobble_ratio = s1 > s2 ? (s1/s2) : (s2/s1);
    return wobble_ratio;
}

HRESULT LamdaUpdator::SmoothingLamdas(vector<float>& lamdaList, int kernelsize, int iteration)
{
    VT_HR_BEGIN()

    float sigma = sqrt(1.0f*kernelsize);
    float sum = 0;
    const int windowSize = 2 * kernelsize + 1;
    vector<float> weights;
    VT_HR_EXIT( weights.resize(windowSize) );
    float alpha = -1.0f / (2.0f * sigma * sigma);
    for ( int i = -kernelsize; i <= kernelsize; i ++)
    {
        float dist = (float)i * (float)i;
        dist = exp(dist * alpha);
        weights[i+kernelsize] = dist;
        sum += dist;
    }
    float sum_inv = 1.0f / sum;
    for ( int i = 0; i < windowSize; i ++)
    {
        weights[i] *= sum_inv;
    }

    const int frameNum = (int) lamdaList.size();

    vector<float> lamdaList_backup;
    VT_HR_EXIT( lamdaList_backup.resize(frameNum) );

    //smooth signal by iterative Gaussian filtering
    for(int k = 0; k < iteration; k++)
    {
        memcpy(&lamdaList_backup[0], &lamdaList[0], sizeof(float) * frameNum);
        memset(&lamdaList[0], 0, sizeof(float) * frameNum);

        for(int i = 0; i < frameNum; i++)
        {
            float cumweight = 0;
            for(int j = -kernelsize; j <= kernelsize; j++)
            {
                int curIdx = i + j;
                if(curIdx > 0 && curIdx < frameNum)
                {
                    lamdaList[i] += weights[j+kernelsize] * lamdaList_backup[curIdx];   
                    cumweight += weights[j+kernelsize];
                }
            }
            lamdaList[i] /= cumweight;
        }
    }

    VT_HR_END()
}

/////////////////////////////////////////////////////////////////
///Implementation of CWLSSmoother
/////////////////////////////////////////////////////////////////
HRESULT CWLSSmoother::Begin()
{
    VT_HR_BEGIN()

    m_optimizationCount = 0;

    m_updateSet.clear();
    m_warpsetCache.clear();
    VT_HR_EXIT( m_warpsetCache.resize(m_params.cacheSize) );

    // allocate the rolling buffers
    VT_HR_EXIT( m_bufTrackerData.resize(2*m_params.iBufferSize) );
    VT_HR_EXIT( m_bufResults.resize(m_params.iBufferSize) );

    VT_HR_END()
}

BUFFER_RANGE CWLSSmoother::GetResultsRange()
{
    BUFFER_RANGE r;
    r.frame_count = m_bufResults.get_available_count();
    r.first_frame = m_bufResults.get_first_id();
    return r;
}

HRESULT CWLSSmoother::PushFeatures(const vt::PointMatch* pMatches, int iMatchCount)
{
    VT_HR_BEGIN()

    // compute the transformation
    CMtx3x3d mtxCur_d;
    if ( iMatchCount == 0 || m_bufTrackerData.get_total_count()==0 )
    {
        mtxCur_d.MakeI();
    }
    else
    {
        if (m_params.motionModel == 4)
        {
            VT_HR_EXIT( VtHomographyFromPointMatches2D(mtxCur_d, pMatches, iMatchCount) );
        }
        else
        {
            VT_HR_EXIT( VtAffineFromPointMatches2D(mtxCur_d, pMatches, iMatchCount) );
        }
    }

    CMtx3x3f mtxCur;
    VtConvertMtx(mtxCur_d, mtxCur);

    m_bufTrackerData.advance(mtxCur);

    if ( m_bufTrackerData.get_total_count() >= m_params.iBufferSize )
    {
        m_bufResults.advance();
        VT_HR_EXIT( AdaptiveSmoothTransform(
            m_bufResults.last(), m_bufResults.get_last_id()) );
    }

    VT_HR_END()
}

HRESULT CWLSSmoother::GetResult(vt::CMtx3x3f& xfrm, int frameNumber)
{
    if( !GetResultsRange().InRange(frameNumber) )
    {
        return E_INVALIDARG;
    }
    xfrm = m_bufResults[frameNumber];
    return S_OK;
}

HRESULT CWLSSmoother::End()
{
    VT_HR_BEGIN()
    
    while( m_bufResults.get_total_count() < m_bufTrackerData.get_total_count() )
    {
        m_bufResults.advance();
        VT_HR_EXIT( AdaptiveSmoothTransform(
            m_bufResults.last(), m_bufResults.get_last_id()) );
    }

    VT_HR_END()
}

HRESULT CWLSSmoother::AdaptiveSmoothTransform(CMtx3x3f& result, int iDst)
{
    VT_HR_BEGIN()
    
    vt::CMtx3x3f mtxC2LT = MatrixCenterToTopLeft(m_params.frameWidth, m_params.frameHeight);
    vt::CMtx3x3f mtxLT2C = MatrixTopLeftToCenter(m_params.frameWidth, m_params.frameHeight);
    float max_sizedim = (float)VtMax(m_params.frameWidth, m_params.frameHeight);
    vt::CMtx3x3f mtxScale = vt::CMtx3x3f().MakeDiag(vt::CVec3f(1.0f/max_sizedim, 1.0f/max_sizedim, 1.0f));
    vt::CMtx3x3f mtxScaleInv = mtxScale.Inv();

    if(m_bufTrackerData.get_total_count() <= m_params.iBufferSize)
    {
        // for short video (frame length < buffer size), do whole sequence smoothing
        if(m_updateSet.size() == 0)
        {
            // prepare motion list
            vt::vector<vt::CMtx3x3f> FList;
            VT_HR_EXIT( FList.resize(m_bufTrackerData.get_total_count()) );
            for(int i = 0; i < m_bufTrackerData.get_total_count(); i ++)
            {
                const vt::CMtx3x3f& motionF = m_bufTrackerData[i];
                FList[i] = mtxScale * motionF * mtxScaleInv; // use normalized transform for our smoothing
            }

            // smooth the path
            VT_HR_EXIT( Smooth(FList, 0, 0) );
        }

        result = m_updateSet[iDst];
    }
    else if(iDst <= m_bufTrackerData.get_last_id() - m_params.iBufferSize + m_params.outFramesPerOptimization)
    {
        // video length greater than buffer size.

        // do smoothing at every several frames
        // for tail frames, do batch smoothing for each one.
        if( (iDst % m_params.outFramesPerOptimization == 0) )
        {
            // decide the frame range to apply path smoothing.
            // we'll use at most 2 * buffer_size frames to do path smoothing.
            int start = VtMax(iDst - m_params.iBufferSize, 0);
            int end   = VtMin(iDst + m_params.iBufferSize - 1, m_bufTrackerData.get_last_id());

            // prepare motion list
            vt::vector<vt::CMtx3x3f> FList;
            VT_HR_EXIT( FList.resize(end - start + 1) );
            for(int i = start; i <= end; i ++)
            {
                const vt::CMtx3x3f& motionF = m_bufTrackerData[i];
                FList[i - start] = mtxScale * motionF * mtxScaleInv; // use normalized transform for our smoothing
            }

            // smooth the path
            VT_HR_EXIT( Smooth(FList, start, iDst - start) );
        }

        // retrieve warping matrix
        result = GetUpdate(iDst);
    }
    else
    {
        VT_ASSERT(m_warpsetCache.size() > 0);

        // for tail frames, we reuse the smoothed path from last optimization.
        size_t cacheIndex = VtMin(m_optimizationCount, m_warpsetCache.size()) - 1;
        size_t index = iDst - m_warpsetCache[cacheIndex].start_index();
        result = m_updateSet[index];		
    }

    result = mtxScaleInv * result * mtxScale;	// restore original scale by denormalization
    
    float crop_ratio = m_params.cropRatio;
    CMtx3x3f trans = vt::CMtx3x3f(1.0f/crop_ratio, 0.0f, (m_params.frameWidth/2.0f)*(1.0f - 1.0f/crop_ratio), 0.0f, 
        1.0f/crop_ratio, (m_params.frameHeight/2.0f)*(1.0f - 1.0f/crop_ratio), 0.0f, 0.0f, 1.0f);
    result = trans * result; // apply cropping scaling

    result = result.Inv();
    
    VT_HR_END()
}

const CMtx3x3f& CWLSSmoother::GetUpdate(size_t frameIndex) const
{
    VT_ASSERT(m_optimizationCount > 0 && m_warpsetCache.size() > 0);
    size_t cacheIndex = VtMin(m_optimizationCount - 1, m_warpsetCache.size() - 1);

    VT_ASSERT(m_warpsetCache[cacheIndex].is_in_range(frameIndex));

    size_t index = frameIndex - m_warpsetCache[cacheIndex].start_index();
    return m_updateSet[index];
}

HRESULT CWLSSmoother::Smooth(const vector<CMtx3x3f>& motionList, size_t frameIndexInVideo, size_t outputFrameIndexInMotionList)
{
    VT_HR_BEGIN()

    if(motionList.size() == 0 || outputFrameIndexInMotionList >= motionList.size())
        hr = E_INVALIDARG;

    // 1. calculate the initial warping matrix collection.
    AdaptivePath clsAdapt;
    VT_HR_EXIT( clsAdapt.Initialize(motionList, m_params.kernelSize, m_params.frameWidth, m_params.frameHeight) );

    // initial smoothing step guided by guided path without considering crop and wobble constraint
    VT_HR_EXIT( clsAdapt.InitialSmoothingLoop(m_nIterations) );
    m_updateSet = clsAdapt.GetUpdateSet();

    LamdaUpdator clsLamdaUpd;
    vector<float> lamdaList;
    VT_HR_EXIT( lamdaList.resize(motionList.size()) );
    memset(&lamdaList[0], 0, sizeof(float) * motionList.size());

    VT_HR_EXIT( clsLamdaUpd.Initialize(lamdaList, m_params.frameWidth, m_params.frameHeight, m_params.kernelSize) );

    // iterate the smoothing and increase lamda until all frames satisfy crop and wobble constraints
    HRESULT flag = clsLamdaUpd.InitUpdateLamdaList(m_updateSet, lamdaList);
    if(flag == S_OK)
    {
        VT_HR_EXIT( clsAdapt.IterativeSmoothingLoop(lamdaList, m_nIterations) );
        m_updateSet = clsAdapt.GetUpdateSet();
    }

    flag = clsLamdaUpd.UpdateLamdaList(m_updateSet, lamdaList);
    while(flag == S_OK)
    {
        VT_HR_EXIT( clsAdapt.IterativeSmoothingLoop(lamdaList, m_nIterations) );
        m_updateSet = clsAdapt.GetUpdateSet();
        flag = clsLamdaUpd.UpdateLamdaList(m_updateSet, lamdaList);
    }

    VT_HR_EXIT( flag );

    // 2. average them with history for post-smoothing.
    {
        size_t cacheToUpdate = vt::VtMin(m_optimizationCount, m_warpsetCache.size() - 1);
        size_t startIndexOfCurrentSet = frameIndexInVideo;
        size_t endIndexOfCurrentSet = startIndexOfCurrentSet + m_updateSet.size();

        for(size_t i = startIndexOfCurrentSet; i < endIndexOfCurrentSet; i ++)
        {
            int count = 1;
            vt::CMtx3x3f& mat = m_updateSet[i - frameIndexInVideo];
            for(size_t c = 0; c < cacheToUpdate; c ++)
            {
                const WarpsetWrapper& wrapper = m_warpsetCache[c];
                if(wrapper.is_in_range(i))
                {
                    mat += wrapper[i];
                    count ++;
                }
            }

            if(count > 1) mat /= (float)count;
        }
    }

    // backup current updates
    vector<CMtx3x3f> updateSet_Backup(m_updateSet);

    // 3. apply post-smoothing to further remove residual jitters
    vector<float> lamdaList_postFilt;
    VT_HR_EXIT( PostSmoothing(updateSet_Backup, motionList, lamdaList_postFilt) );

    // update lamda value in post-smoothing according to hard crop constraint
    VT_HR_EXIT( lamdaList_postFilt.resize(motionList.size()) );
    memset(&lamdaList_postFilt[0], 0, sizeof(float) * lamdaList_postFilt.size());
    VT_HR_EXIT( clsLamdaUpd.UpdateLamdaByCropHardConstraint(updateSet_Backup, lamdaList_postFilt) );

    // apply post-smoothing again with hard crop constraint
    VT_HR_EXIT( PostSmoothing(m_updateSet, motionList, lamdaList_postFilt) );

    // 4. update cache
    {
        size_t cacheToUpdate = m_optimizationCount;
        if(m_optimizationCount >= m_warpsetCache.size())
        {
            // the oldest cache should be expired
            for(size_t i = 0; i < m_warpsetCache.size() - 1; i ++)
            {
                m_warpsetCache[i].swap(m_warpsetCache[i + 1]);
            }

            // append updates from this optimization to end of the cache set.
            cacheToUpdate = m_warpsetCache.size() - 1;
        }

        m_warpsetCache[cacheToUpdate].initialize(m_updateSet, frameIndexInVideo);

        m_optimizationCount ++;
    }

    VT_HR_END()
}

HRESULT CWLSSmoother::PostSmoothing(
    vector<CMtx3x3f>& updateSet,
    const vector<CMtx3x3f>& motionList,
    const vector<float>& lamdaList
    )
{
    VT_HR_BEGIN()

    vector<CMtx3x3f> HList; // restore H_i
    VT_HR_EXIT( HList.resize(motionList.size()) ); 
    HList[0] = CMtx3x3f().MakeI();

    CMtx3x3f ivt;
    for(size_t i = 1; i < motionList.size(); i++)
    {
        HList[i] = updateSet[i].Inv() * motionList[i] * updateSet[i - 1];
    }

    BasicPath clsBasic;
    const int kernelsize = VtMax(1, m_params.kernelSize * 2 / 3);
    const int iteration = 3;
    VT_HR_EXIT( clsBasic.Initialize(HList, kernelsize, m_params.frameWidth, m_params.frameHeight) );
    VT_HR_EXIT( clsBasic.BasicSmoothing(iteration, lamdaList) );
    const vector<CMtx3x3f>& updateSetRatio = clsBasic.GetUpdateSet();

    for(size_t i = 0; i < motionList.size(); i++)
    {
        updateSet[i] = updateSet[i] * updateSetRatio[i];
    }

    VT_HR_END()
}
