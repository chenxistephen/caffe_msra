#include <Windows.h>
#include <vector>
#include <string>
#include <math.h>
#include "StructuredEdgeDetector.h"

#include "LineSegmentDetector.h"

//#include "vtfileio.h"   // DEBUG ONLY

namespace edgebox
{
    CStructuredEdgeDetector::CStructuredEdgeDetector() : m_seInitialized(FALSE)
    {
    }    

    HRESULT CStructuredEdgeDetector::DetectEdges(const vt::CRGBImg& imgColor, 
        MexImageFloat& imgE,
        MexImageFloat& imgO)
    {
        if (!imgColor.IsValid())
        {
            return E_INVALIDARG;
        }

        if (m_seInitialized)
        {
            ValidateModel();

            const auto& opts = m_detectorModel.opts;

            if ((BOOL)opts.multiscale)
            {
                // TODO: translate matlab code for multiscale case, the model file currently turns this case off
                return E_UNEXPECTED;
            }
            else
            {
                // pad image (new image size to be divisible by 4)
                const int cx = imgColor.Width();
                const int cy = imgColor.Height();

                const int pad = (int)opts.imWidth / 2;

                int padAmount[4] = { pad, pad, pad, pad };   // top, bottom, left, right

                padAmount[1] += (4 - ((cy + 2 * pad) % 4)) % 4;
                padAmount[3] += (4 - ((cx + 2 * pad) % 4)) % 4;

                MexImageByte imageMex;

                CMexHelper::ConvertVt2Mex(imgColor, imageMex);

                MexImageByte paddedImageMex;

                PadMexImage(imageMex, padAmount, PadMethod::Symmetric, paddedImageMex);

                MexImageFloat chnsReg;
                MexImageFloat chnsSim;

                ComputeEdgeFeatures(paddedImageMex, chnsReg, chnsSim);

                MexImageFloat E;    // range[0, 255]
                MexImage<UINT32> inds;

                const double s = opts.sharpen;

                if (!IS_ZERO(s))
                {
                    MexImageFloat paddedFloat;

                    CMexHelper::rgbConvert(paddedImageMex, L"rgb", paddedFloat);

                    MexImageFloat paddedFloatConv;

                    CMexHelper::convTri(paddedFloat, 1, paddedFloatConv);

                    DetectEdgesImpl(paddedFloatConv,
                        chnsReg,
                        chnsSim,
                        E,
                        inds);
                }
                else
                {
                    // TODO: handle s is 0 case, but edgeDetectMex seems to only handle floating point input image
                }

                float t = float(pow(opts.stride, 2) / (pow(opts.gtWidth, 2) * opts.nTreesVal));
                const int r = int(opts.gtWidth / 2);

                if (IS_ZERO(s))
                {
                    t = t * 2;
                }
                else if (IS_EQUAL(s, 1))
                {
                    t = t * 1.8f;
                }
                else
                {
                    t = t * 1.66f;
                }

                // crop and scale
                MexImageFloat croppedE;
                croppedE.Initialize(E.width - 2 * r, E.height - 2 * r, E.channels);

                // TODO: handle multiple channels, here E.channels is 1!
                for (int i = 0; i < croppedE.buffer.size(); ++i)
                {
                    const int ix = i / croppedE.height;
                    const int iy = i - ix * croppedE.height;

                    croppedE.buffer[i] = E.buffer[(ix + r) * E.height + (iy + r)] * t;
                }   // i

                CMexHelper::convTri(croppedE, 1, imgE);

                if (IS_EQUAL(opts.nms, -1))
                {

                }
                else //todo: imgO should probably be optional in interface
                {
                    GetImgO(imgE, 4, imgO);
                }

                if (opts.nms > 0)
                {
                    // TODO: E=edgesNmsMex(E,O,1,5,1.01,opts.nThreads);
                }

                // TODO: not sure why matlab has such code
                //if (false)
                //{
                //    E = gradientMag(convTri(single(I), 4));
                //    E = E / max(E(:));
                //}

                CMexHelper::edgesNms(imgE, imgO, 2, 0, 1, (int)opts.nThreads);
            }
        }
        else
        {
            vt::CByteImg imGray;
            imGray.Create(imgColor.Width(), imgColor.Height());

            vt::VtConvertImageToLuma(imGray, imgColor);

            WhiteboardCleanup::CLineSegmentDetector lineDetector;

            std::vector<POINT> edgePts;

            lineDetector.DetectEdges(imGray, true, edgePts);

            imgE.Initialize(imgColor.Width(), imgColor.Height(), 1, true);

            const int height = imGray.Height();

            for (size_t i = 0; i < edgePts.size(); ++i)
            {
                const int x = edgePts[i].x;
                const int y = edgePts[i].y;

                imgE.buffer[x * height + y] = 1.f;
            }

            GetImgO(imgE, 4, imgO);
        }

        return S_OK;
    }

    HRESULT CStructuredEdgeDetector::GetImgO(const MexImageFloat& imgE,
        double r,
        MexImageFloat& imgO)
    {
        MexImageFloat imgEFilt;

        CMexHelper::convTri(imgE, r, imgEFilt);

        MexImageFloat Ox, Oy;

        CMexHelper::gradient2(imgEFilt, Ox, Oy);

        MexImageFloat Oxx, temp;

        CMexHelper::gradient2(Ox, Oxx, temp);

        MexImageFloat Oxy, Oyy;

        CMexHelper::gradient2(Oy, Oxy, Oyy);

        imgO.Initialize(imgE.width, imgE.height, imgE.channels);

        for (size_t i = 0; i < imgO.buffer.size(); ++i)
        {
            const float num = IS_ZERO(Oxy.buffer[i]) ? 0 : (Oxy.buffer[i] > 0 ? -Oyy.buffer[i] : Oyy.buffer[i]);
            const float denorm = Oxx.buffer[i] + 1e-5f;
            const float val = atan2(num, denorm);

            imgO.buffer[i] = val - floor(val / PI) * PI;    // different behavior of C++::fmod and matlab::mod                    
        }   // it

        return S_OK;
    }

    HRESULT CStructuredEdgeDetector::DetectEdgesImpl(const MexImageFloat& imgSrc,
        const MexImageFloat& chnsReg,
        const MexImageFloat& chnsSim,
        MexImageFloat& imgE,
        MexImage<UINT32>& imgInds, 
        UINT8** pImgSegs)
    {       
        float *I = (float*)&imgSrc.buffer[0];
        float *chns = (float*)&chnsReg.buffer[0];
        float *chnsSs = (float*)&chnsSim.buffer[0];

        // extract relevant fields from model and options
        float *thrs = &m_detectorModel.thrs[0]; 
        UINT32 *fids = &m_detectorModel.fids[0];
        UINT32 *child = &m_detectorModel.child[0];
        UINT8 *segs = &m_detectorModel.segs[0];
        UINT8 *nSegs = &m_detectorModel.nSegs[0];
        UINT16 *eBins = &m_detectorModel.eBins[0];
        UINT32 *eBnds = &m_detectorModel.eBnds[0];;
        
        const auto& opts = m_detectorModel.opts;
        
        const int shrink = (int)opts.shrink;
        const int imWidth = (int)opts.imWidth;
        const int gtWidth = (int)opts.gtWidth;
        const int nChns = (int)opts.nChns;
        const int nCells = (int)opts.nCells;
        const UINT32 nChnFtrs = (UINT32)opts.nChnFtrs;
        const int stride = (int)opts.stride;
        const int nTreesEval = (int)opts.nTreesVal;
        int sharpen = (int)opts.sharpen;
        int nThreads = (int)opts.nThreads;
        const int nBnds = int(m_detectorModel.eBnds.size() - 1) /
            int(m_detectorModel.thrs.size());
        //const char *msgSharpen = "Model supports sharpening of at most %i pixels!\n";
        if (sharpen > nBnds - 1) 
        { 
            sharpen = nBnds - 1; 
            //mexPrintf(msgSharpen, sharpen); 
        }

        // get dimensions and constants
        const int h = imgSrc.height;
        const int w = imgSrc.width;
        const int Z = imgSrc.channels;
        const int nTreeNodes = m_fidsDim[0];
        const int nTrees = m_fidsDim[1];
        const int h1 = (int)ceil(double(h - imWidth) / stride);
        const int w1 = (int)ceil(double(w - imWidth) / stride);
        const int h2 = h1*stride + gtWidth;
        const int w2 = w1*stride + gtWidth;
        const int imgDims[3] = { h, w, Z };
        const int chnDims[3] = { h / shrink, w / shrink, nChns };
        const int indDims[3] = { h1, w1, nTreesEval };
        const int outDims[3] = { h2, w2, 1 };
        const int segDims[5] = { gtWidth, gtWidth, h1, w1, nTreesEval };

        // construct lookup tables
        UINT32 *iids, *eids, *cids, *cids1, *cids2;
        iids = CMexHelper::buildLookup((int*)imgDims, gtWidth);
        eids = CMexHelper::buildLookup((int*)outDims, gtWidth);
        cids = CMexHelper::buildLookup((int*)chnDims, imWidth / shrink);
        CMexHelper::buildLookupSs(cids1, cids2, (int*)chnDims, imWidth / shrink, nCells);

        // create outputs
        imgE.Initialize(outDims[1], outDims[0], outDims[2]);
        
        float *E = &imgE.buffer[0];

        imgInds.Initialize(indDims[1], indDims[0], indDims[2]);
        
        UINT32 *ind = &imgInds.buffer[0];

        if (pImgSegs != nullptr)
        {
            // TODO: create the array as in matlab code (5D, should probably re-design the MexImage structure)
        }
                
        UINT8 *segsOut = nullptr; 
        if (pImgSegs != nullptr)
        {
            //TODO: segsOut = (uint8*)mxGetData(pl[2]);
        }

        // apply forest to all patches and store leaf inds
#ifdef USEOMP
        nThreads = min(nThreads, omp_get_max_threads());
#pragma omp parallel for num_threads(nThreads)
#endif
        for (int c = 0; c < w1; c++) for (int t = 0; t < nTreesEval; t++) {
            for (int r0 = 0; r0 < 2; r0++) for (int r = r0; r < h1; r += 2) {
                int o = (r*stride / shrink) + (c*stride / shrink)*h / shrink;
                // select tree to evaluate
                int t1 = ((r + c) % 2 * nTreesEval + t) % nTrees; UINT32 k = t1*nTreeNodes;
                while (child[k]) {
                    // compute feature (either channel or self-similarity feature)
                    UINT32 f = fids[k]; float ftr;
                    if (f < nChnFtrs) ftr = chns[cids[f] + o]; else
                        ftr = chnsSs[cids1[f - nChnFtrs] + o] - chnsSs[cids2[f - nChnFtrs] + o];
                    // compare ftr to threshold and move left or right accordingly
                    if (ftr < thrs[k]) k = child[k] - 1; else k = child[k];
                    k += t1*nTreeNodes;
                }
                // store leaf index and update edge maps
                ind[r + c*h1 + t*h1*w1] = k;
            }
        }

        // compute edge maps (avoiding collisions from parallel executions)
        if (!sharpen) for (int c0 = 0; c0 < gtWidth / stride; c0++) {
#ifdef USEOMP
#pragma omp parallel for num_threads(nThreads)
#endif
            for (int c = c0; c < w1; c += gtWidth / stride) {
                for (int r = 0; r < h1; r++) for (int t = 0; t < nTreesEval; t++) {
                    UINT32 k = ind[r + c*h1 + t*h1*w1];
                    float *E1 = E + (r*stride) + (c*stride)*h2;
                    int b0 = eBnds[k*nBnds], b1 = eBnds[k*nBnds + 1]; if (b0 == b1) continue;
                    for (int b = b0; b < b1; b++) E1[eids[eBins[b]]]++;
                    if (pImgSegs != nullptr)
                    {
                        // TODO: memcpy(segsOut + (r + c*h1 + t*h1*w1)*gtWidth*gtWidth, segs + k*gtWidth*gtWidth, gtWidth*gtWidth*sizeof(uint8));
                    }
                }
            }
        }

        // computed sharpened edge maps, snapping to local color values
        if (sharpen) {
            // compute neighbors array
            const int g = gtWidth; UINT16 N[4096 * 4];
            for (int c = 0; c < g; c++) for (int r = 0; r < g; r++) {
                int i = c*g + r; UINT16 *N1 = N + i * 4;
                N1[0] = c > 0 ? i - g : i; N1[1] = c < g - 1 ? i + g : i;
                N1[2] = r > 0 ? i - 1 : i; N1[3] = r < g - 1 ? i + 1 : i;
            }
#ifdef USEOMP
#pragma omp parallel for num_threads(nThreads)
#endif
            for (int c = 0; c < w1; c++) for (int r = 0; r < h1; r++) {
                for (int t = 0; t < nTreesEval; t++) {
                    // get current segment and copy into S
                    UINT32 k = ind[r + c*h1 + t*h1*w1];
                    int m = nSegs[k]; if (m == 1) continue;
                    UINT8 S0[4096], *S = (pImgSegs == nullptr) ? S0 : segsOut + (r + c*h1 + t*h1*w1)*g*g;
                    memcpy(S, segs + k*g*g, g*g*sizeof(UINT8));
                    // compute color model for each segment using every other pixel
                    int ci, ri, s, z; float ns[100], mus[1000];
                    const float *I1 = I + (c*stride + (imWidth - g) / 2)*h + r*stride + (imWidth - g) / 2;
                    for (s = 0; s < m; s++) { ns[s] = 0; for (z = 0; z < Z; z++) mus[s*Z + z] = 0; }
                    for (ci = 0; ci < g; ci += 2) for (ri = 0; ri < g; ri += 2) {
                        s = S[ci*g + ri]; ns[s]++;
                        for (z = 0; z < Z; z++) mus[s*Z + z] += I1[z*h*w + ci*h + ri];
                    }
                    for (s = 0; s < m; s++) for (z = 0; z < Z; z++) mus[s*Z + z] /= ns[s];
                    // update segment S according to local color values
                    int b0 = eBnds[k*nBnds], b1 = eBnds[k*nBnds + sharpen];
                    for (int b = b0; b < b1; b++) {
                        float vs[10], d, e, eBest = 1e10f; int i, sBest = -1, ss[4];
                        for (i = 0; i < 4; i++) ss[i] = S[N[eBins[b] * 4 + i]];
                        for (z = 0; z < Z; z++) vs[z] = I1[iids[eBins[b]] + z*h*w];
                        for (i = 0; i < 4; i++) {
                            s = ss[i]; if (s == sBest) continue;
                            e = 0; for (z = 0; z < Z; z++) { d = mus[s*Z + z] - vs[z]; e += d*d; }
                            if (e < eBest) { eBest = e; sBest = s; }
                        }
                        S[eBins[b]] = sBest;
                    }
                    // convert mask to edge maps (examining expanded set of pixels)
                    float *E1 = E + c*stride*h2 + r*stride; b1 = eBnds[k*nBnds + sharpen + 1];
                    for (int b = b0; b < b1; b++) {
                        int i = eBins[b]; UINT8 s = S[i]; UINT16 *N1 = N + i * 4;
                        if (s != S[N1[0]] || s != S[N1[1]] || s != S[N1[2]] || s != S[N1[3]])
                            E1[eids[i]]++;
                    }
                }
            }
        }

        // free memory
        delete[] iids; delete[] eids;
        delete[] cids; delete[] cids1; delete[] cids2;

        return S_OK;
    }

    HRESULT CStructuredEdgeDetector::ComputeEdgeFeatures(const MexImageByte& imgSrc,
        MexImageFloat& chnsReg,
        MexImageFloat& chnsSim)
    {
        const auto& opts = m_detectorModel.opts;

        const double shrink = opts.shrink;

        const int nTypes = 1; // TODO: translate nTypes = 2 case in the matlab code

        int k = 0;

        MexImageFloat chns; // the cell in matlab

        for (int t = 1; t <= nTypes; ++t)
        {
            const wchar_t* cs = imgSrc.channels == 1 ? L"gray" : L"luv";

            //std::vector<std::vector<float>> chns;

            MexImageFloat imgLuv;

            CMexHelper::rgbConvert(imgSrc, cs, imgLuv);

            MexImageFloat imgLuvResample;

            CMexHelper::imResample(imgLuv, 1 / shrink, imgLuvResample);

            ++k;
            CMexHelper::AppendImage(chns, imgLuvResample);

            for (int i = 1; i <= 2; ++i)
            {
                const double s = pow(2, i - 1);

                MexImageFloat I1;

                if (IS_EQUAL(s, shrink))
                {
                    I1 = imgLuvResample;
                }
                else
                {
                    CMexHelper::imResample(imgLuv, 1 / s, I1);
                }

                MexImageFloat convI1;

                CMexHelper::convTri(I1, opts.grdSmooth, convI1);

                MexImageFloat M, O;

                CMexHelper::gradientMag(convI1, M, &O, 0, opts.normRad, 0.01);

                MexImageFloat H;

                CMexHelper::gradientHist(M, O, H, (int)max(1, shrink / s), (int)opts.nOrients, 0);                

                MexImageFloat M_resample;
                CMexHelper::imResample(M, s / shrink, M_resample);

                ++k;
                CMexHelper::AppendImage(chns, M_resample);                

                MexImageFloat H_resample;
                CMexHelper::imResample(H, max(1, s / shrink), H_resample);

                ++k;
                CMexHelper::AppendImage(chns, H_resample);
            }
        }   // t

        double chnSm = opts.chnSmooth / shrink;
        if (chnSm > 1)
        {
            chnSm = double(int(chnSm + 0.5));
        }

        CMexHelper::convTri(chns, chnSm, chnsReg);

        double simSm = opts.simSmooth / shrink;
        if (simSm > 1)
        {
            simSm = double(int(simSm + 0.5));
        }

        CMexHelper::convTri(chns, simSm, chnsSim);

        return S_OK;
    }

    HRESULT CStructuredEdgeDetector::PadMexImage(const MexImageByte& imgSrc,        
        const int padAmount[],
        PadMethod method,
        MexImageByte& imgDst)
    {
        imgDst.Initialize(imgSrc.width + padAmount[2] + padAmount[3],
            imgSrc.height + padAmount[0] + padAmount[1],
            imgSrc.channels);

        int flag = 0;

        switch (method)
        {
        case PadMethod::Replicate:
            flag = 1;
            break;
        case PadMethod::Symmetric:
            flag = 2;
            break;
        case PadMethod::Circular:
            flag = 3;
            break;
        }

        CMexHelper::imPad((byte*)&imgSrc.buffer[0], 
            &imgDst.buffer[0], 
            imgSrc.height, 
            imgSrc.width, 
            imgSrc.channels, 
            padAmount[0], 
            padAmount[1], 
            padAmount[2], 
            padAmount[3], 
            flag, 
            (byte)0);

        return S_OK;
    }

    void CStructuredEdgeDetector::ValidateModel()
    {
        m_detectorModel.opts.nTreesVal = min(m_detectorModel.opts.nTreesVal, m_detectorModel.opts.nTrees);

        // TODO: in matlab, we checked the existence of "sharpen" and "segs"

        m_detectorModel.opts.stride = max(m_detectorModel.opts.stride, m_detectorModel.opts.shrink);
    }

    HRESULT CStructuredEdgeDetector::LoadDetectorModel(const wchar_t* szModelPath)
    {
        FILE* pFile = nullptr;

        _wfopen_s(&pFile, szModelPath, L"rb");

        // refer matlab variable window for details

        // read opts
        fread(&m_detectorModel.opts.imWidth, sizeof(double), 10, pFile);

        int nLength = 0;

        fread(&nLength, sizeof(int), 1, pFile);
        m_detectorModel.opts.discretize.resize(nLength);
        fread(&m_detectorModel.opts.discretize[0], sizeof(char), nLength, pFile);

        fread(&m_detectorModel.opts.nSamples, sizeof(double), 2, pFile);

        fread(&nLength, sizeof(int), 1, pFile);
        m_detectorModel.opts.split.resize(nLength);
        fread(&m_detectorModel.opts.split[0], sizeof(char), nLength, pFile);

        fread(&m_detectorModel.opts.nOrients, sizeof(double), 16, pFile);

        fread(&nLength, sizeof(int), 1, pFile);
        m_detectorModel.opts.modelDir.resize(nLength);
        fread(&m_detectorModel.opts.modelDir[0], sizeof(char), nLength, pFile);

        fread(&nLength, sizeof(int), 1, pFile);
        m_detectorModel.opts.modelFnm.resize(nLength);
        fread(&m_detectorModel.opts.modelFnm[0], sizeof(char), nLength, pFile);

        fread(&nLength, sizeof(int), 1, pFile);
        m_detectorModel.opts.bsdsDir.resize(nLength);
        fread(&m_detectorModel.opts.bsdsDir[0], sizeof(char), nLength, pFile);

        fread(&m_detectorModel.opts.nChns, sizeof(double), 4, pFile);

        // read arrays
        int dim[4];
        
        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.thrs.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.thrs[0], sizeof(float), m_detectorModel.thrs.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        memcpy(m_fidsDim, dim, 4 * sizeof(int));
        m_detectorModel.fids.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.fids[0], sizeof(UINT32), m_detectorModel.fids.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.child.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.child[0], sizeof(UINT32), m_detectorModel.child.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.count.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.count[0], sizeof(UINT32), m_detectorModel.count.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.depth.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.depth[0], sizeof(UINT32), m_detectorModel.depth.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.segs.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.segs[0], sizeof(byte), m_detectorModel.segs.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.nSegs.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.nSegs[0], sizeof(byte), m_detectorModel.nSegs.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.eBins.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.eBins[0], sizeof(UINT16), m_detectorModel.eBins.size(), pFile);

        fread(dim, sizeof(int), 4, pFile);
        m_detectorModel.eBnds.resize(dim[0] * dim[1] * dim[2] * dim[3]);
        fread(&m_detectorModel.eBnds[0], sizeof(UINT32), m_detectorModel.eBnds.size(), pFile);

        fclose(pFile);

        m_seInitialized = TRUE;

        return S_OK;
    }
}