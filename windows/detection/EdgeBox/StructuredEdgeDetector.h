#pragma once

#include "vtcore.h"
#include "MexHelper.hpp"

namespace edgebox
{
    class CStructuredEdgeDetector
    {
    public:
        CStructuredEdgeDetector();

        HRESULT LoadDetectorModel(const wchar_t* szModelPath);

        HRESULT DetectEdges(const vt::CRGBImg& imgColor, 
            MexImageFloat& imgE, 
            MexImageFloat& imgO);            

    private:        
        BOOL m_seInitialized;   // if structured edge detector has been initialized

    private:        
        void ValidateModel();

        enum PadMethod
        {
            Replicate = 0,
            Symmetric = 1,
            Circular = 2
        };

        HRESULT PadMexImage(const MexImageByte& imgSrc, 
            const int padAmount[],
            PadMethod method,
            MexImageByte& imgDst);

        HRESULT ComputeEdgeFeatures(const MexImageByte& imgSrc,
            MexImageFloat& chnsReg,
            MexImageFloat& chnsSim);

        HRESULT DetectEdgesImpl(const MexImageFloat& imgSrc,
            const MexImageFloat& chnsReg,
            const MexImageFloat& chnsSim,
            MexImageFloat& imgE,
            MexImage<UINT32>& imgInds, 
            UINT8** pImgSegs = nullptr);    

        HRESULT GetImgO(const MexImageFloat& imgE,
            double r,
            MexImageFloat& imgO);

    private:
        struct ModelOptions
        {
            double imWidth;
            double gtWidth;
            double nPos;
            double nNeg;
            double nImgs;
            double nTrees;
            double fracFtrs;
            double minCount;
            double minChild;
            double maxDepth;
            std::string discretize;
            double nSamples;
            double nClasses;
            std::string split;
            double nOrients;
            double grdSmooth;
            double chnSmooth;
            double simSmooth;
            double normRad;
            double shrink;
            double nCells;
            double rgbd;
            double stride;
            double multiscale;
            double sharpen;
            double nTreesVal;
            double nThreads;
            double nms;
            double seed;
            double useParfor;
            std::string modelDir;
            std::string modelFnm;
            std::string bsdsDir;
            double nChns;
            double nChnFtrs;
            double nSimFtrs;
            double nTotFtrs;
        };

        struct DetectorModel
        {
            ModelOptions opts;

            std::vector<float> thrs;
            std::vector<UINT32> fids;
            std::vector<UINT32> child;
            std::vector<UINT32> count;
            std::vector<UINT32> depth;
            std::vector<byte> segs;
            std::vector<byte> nSegs;
            std::vector<UINT16> eBins;
            std::vector<UINT32> eBnds;
        };        

        DetectorModel m_detectorModel;
        int m_fidsDim[4];
    };
}