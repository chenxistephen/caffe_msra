#pragma once

#include "vtcore.h"
#include <vector>
#include "StructuredEdgeDetector.h"

namespace edgebox
{
    class CEdgeBox
    {
    public:
        // passing in nullptr -> using simple edge detector
        CEdgeBox(const wchar_t* szEdgeModel = nullptr);

        HRESULT GetBoxes(const vt::CRGBImg& imgSrc,
            std::vector<RECT>& rgBoxes);

    public:
        struct Context
        {
            double alpha;                   // step size of sliding windows search
            double beta;                    // nms threshold for object proposals
            double minScore;                // min score of boxes to detect
            double maxBoxes;                // max number of boxes to detect
            double edgeMinMag;
            double edgeMergeThr;
            double clusterMinMag;
            double maxAspectRatio;
            double minBoxArea;
            double gamma;
            double kappa;

            int maxDim;
        };

    public:
        void GetContext(Context* pContext) const { memcpy(pContext, &m_context, sizeof(Context)); }

        void SetContext(const Context* pContext) { memcpy(&m_context, pContext, sizeof(Context)); }

    private:
        Context m_context;

        CStructuredEdgeDetector m_edgeDetector;

    private:
        void FillDefaultContext();
    };
}