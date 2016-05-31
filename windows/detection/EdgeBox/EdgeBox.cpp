#include "EdgeBox.h"
#include "EdgeBoxGenerator.h"

namespace edgebox
{
    CEdgeBox::CEdgeBox(const wchar_t* szEdgeModel)
    {
        FillDefaultContext();

        if (szEdgeModel != nullptr)
        {
            m_edgeDetector.LoadDetectorModel(szEdgeModel);
        }
    }

    void CEdgeBox::FillDefaultContext()
    {
        // alpha, beta are crucial for making trade-off between speed and accuracy
        // EdgeBox50: 0.65, 0.55
        // EdgeBox70: 0.65, 0.75    (best trade-off)
        // EdgeBox90: 0.85, 0.95    (slowest)
        m_context.alpha = 0.65;
        m_context.beta = 0.75;
        
        m_context.minScore = 0.01;
        m_context.maxBoxes = 2000;
        m_context.edgeMinMag = 0.1;
        m_context.edgeMergeThr = 0.5;
        m_context.clusterMinMag = 0.5;
        m_context.maxAspectRatio = 3;
        m_context.minBoxArea = 1000;
        m_context.gamma = 2;
        m_context.kappa = 1.5;

        m_context.maxDim = 256; //384;  // 512  // 256  // 200
    }

    HRESULT CEdgeBox::GetBoxes(const vt::CRGBImg& imgSrc, 
        std::vector<RECT>& edgeBoxes)
    {
        const int width = imgSrc.Width();
        const int height = imgSrc.Height();

        vt::CRGBImg imDs;
        float dsFactor = 1.f;
        if (__max(width, height) > m_context.maxDim)
        {
            dsFactor = float(__max(width, height)) / m_context.maxDim;

            const int dsWidth = int(width / dsFactor);
            const int dsHeight = int(height / dsFactor);

            imDs.Create(dsWidth, dsHeight);

            vt::VtResizeImage(imDs, imDs.Rect(), imgSrc);
        }
        else
        {
            imgSrc.Share(imDs);
        }        

        MexImageFloat imgE, imgO;

        HRESULT hr = m_edgeDetector.DetectEdges(imDs, imgE, imgO);
        if (FAILED(hr)) return hr;        

        EdgeBoxGenerator edgeBoxGen; 
        Boxes boxes;

        edgeBoxGen._alpha = float(m_context.alpha);
        edgeBoxGen._beta = float(m_context.beta);
        edgeBoxGen._minScore = float(m_context.minScore);
        edgeBoxGen._maxBoxes = int(m_context.maxBoxes);
        edgeBoxGen._edgeMinMag = float(m_context.edgeMinMag);
        edgeBoxGen._edgeMergeThr = float(m_context.edgeMergeThr);
        edgeBoxGen._clusterMinMag = float(m_context.clusterMinMag);
        edgeBoxGen._maxAspectRatio = float(m_context.maxAspectRatio);
        edgeBoxGen._minBoxArea = float(m_context.minBoxArea);
        edgeBoxGen._gamma = float(m_context.gamma);
        edgeBoxGen._kappa = float(m_context.kappa);
        
        arrayf E;
        E._x = &imgE.buffer[0];
        E._w = imgE.width;
        E._h = imgE.height;

        arrayf O;
        O._x = &imgO.buffer[0];
        O._w = imgO.width;
        O._h = imgO.height;

        arrayf V;   // dummy, for visualization

        edgeBoxGen.generate(boxes, E, O, V);

        edgeBoxes.resize(boxes.size());

        for (size_t i = 0; i < boxes.size(); ++i)
        {
            edgeBoxes[i] = 
            { 
                int(dsFactor * boxes[i].c),
                int(dsFactor * boxes[i].r),
                int(dsFactor * (boxes[i].c + boxes[i].w - 1)), 
                int(dsFactor * (boxes[i].r + boxes[i].h - 1)) 
            };
        }

        // debug
        //float temp = 0;
        //for (int i = 0; i < 100; ++i)
        //    temp += boxes[i].s;
        // end of debug

        return S_OK;
    }
}
