#pragma once

#include "vtcore.h"
#include "sse.hpp"
#include "MexWrappers.h"
#include <omp.h>
#include <vector>
#include <assert.h>

#define IS_ZERO(x) (abs(x) < 1e-6)
#define IS_EQUAL(x, y) (abs(x - y) < 1e-6)
#define PI 3.14159265f

namespace edgebox
{
    template<typename T>
    struct MexImage
    {
        std::vector<T> buffer;  // column majored
        
        int width;
        int height;        
        int channels;

        MexImage()
        {
            width = height = channels = 0;
        }

        MexImage& operator=(const MexImage& rhs)
        {
            if (&rhs != this)
            {
                buffer = rhs.buffer;
                height = rhs.height;
                width = rhs.width;
                channels = rhs.channels;
            }

            return *this;
        }

        void Initialize(int _width, int _height, int _channels, bool clearInitials = false)
        {
            width = _width;
            height = _height;
            channels = _channels;

            if (clearInitials)
            {
                buffer.resize(width * height * channels, 0);
            }            
            else
            {
                buffer.resize(width * height * channels);
            }
        }
    };

    typedef MexImage<byte> MexImageByte;
    typedef MexImage<float> MexImageFloat;

    class CMexHelper
    {
    public:
        template<typename T>
        static HRESULT AppendImage(MexImage<T>& imgDst,
            const MexImage<T>& imgSrc)
        {
            if (imgDst.buffer.empty())
            {
                imgDst = imgSrc;

                return S_OK;
            }

            if (imgDst.width != imgSrc.width
                || imgDst.height != imgSrc.height)
            {
                return E_INVALIDARG;
            }

            imgDst.channels += imgSrc.channels;

            imgDst.buffer.reserve(imgDst.buffer.size() + imgSrc.buffer.size());
            imgDst.buffer.insert(imgDst.buffer.end(), imgSrc.buffer.begin(), imgSrc.buffer.end());

            return S_OK;
        }

        static HRESULT ConvertVt2Mex(const vt::CRGBImg& imVt,
            MexImageByte& imMex)
        {
            if (!imVt.IsValid())
            {
                return E_INVALIDARG;
            }

            const int cx = imVt.Width();
            const int cy = imVt.Height();

            imMex.Initialize(cx, cy, 3);            

            for (int y = 0; y < cy; ++y)
            {
                const byte* ptrRow = imVt.BytePtr(y);

                for (int x = 0; x < cx; ++x)
                {
                    const byte* ptrPix = &(ptrRow[x * 3]);

                    imMex.buffer[0 * cx * cy + cy * x + y] = ptrPix[2];    // r    
                    imMex.buffer[1 * cx * cy + cy * x + y] = ptrPix[1];    // g
                    imMex.buffer[2 * cx * cy + cy * x + y] = ptrPix[0];    // b
                }   // x
            }   // y

            return S_OK;
        }

        static HRESULT ConvertMex2Vt(const MexImageByte& imMex,
            vt::CRGBImg& imVt)
        {
            if (imMex.channels != 3)
            {
                return E_INVALIDARG;
            }

            const int width = imMex.width;
            const int height = imMex.height;

            HRESULT hr = imVt.Create(width, height);
            if (FAILED(hr)) return hr;

            for (int y = 0; y < height; ++y)
            {
                byte* ptrRow = imVt.BytePtr(y);

                for (int x = 0; x < width; ++x)
                {
                    byte* ptrPix = &(ptrRow[x * 3]);

                    ptrPix[2] = imMex.buffer[0 * width * height + height * x + y];
                    ptrPix[1] = imMex.buffer[1 * width * height + height * x + y];
                    ptrPix[0] = imMex.buffer[2 * width * height + height * x + y];
                }   // x
            }   // y

            return S_OK;
        }

        template<class T>
        static void imPad(T *A, T *B, int h, int w, int d, int pt, int pb,
            int pl, int pr, int flag, T val)
        {
            int h1 = h + pt, hb = h1 + pb, w1 = w + pl, wb = w1 + pr, x, y, z, mPad;
            int ct = 0, cb = 0, cl = 0, cr = 0;
            if (pt<0) { ct = -pt; pt = 0; } if (pb < 0) { h1 += pb; cb = -pb; pb = 0; }
            if (pl < 0) { cl = -pl; pl = 0; } if (pr < 0) { w1 += pr; cr = -pr; pr = 0; }
            int *xs = nullptr, *ys = nullptr; x = pr > pl ? pr : pl; y = pt > pb ? pt : pb; mPad = x > y ? x : y;
            bool useLookup = ((flag == 2 || flag == 3) && (mPad>h || mPad>w))
                || (flag == 3 && (ct || cb || cl || cr));
            // helper macro for padding
#define PAD(XL,XM,XR,YT,YM,YB) \
  for(x=0;  x<pl; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XL+cl)*h+YT+ct]; \
  for(x=0;  x<pl; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XL+cl)*h+YM+ct]; \
  for(x=0;  x<pl; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XL+cl)*h+YB-cb]; \
  for(x=pl; x<w1; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XM+cl)*h+YT+ct]; \
  for(x=pl; x<w1; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XM+cl)*h+YB-cb]; \
  for(x=w1; x<wb; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XR-cr)*h+YT+ct]; \
  for(x=w1; x<wb; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XR-cr)*h+YM+ct]; \
  for(x=w1; x<wb; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XR-cr)*h+YB-cb];
            // build lookup table for xs and ys if necessary
            if (useLookup) {
                xs = (int*)malloc(wb*sizeof(int)); int h2 = (pt + 1) * 2 * h;
                ys = (int*)malloc(hb*sizeof(int)); int w2 = (pl + 1) * 2 * w;
                if (flag == 2) {
                    for (x = 0; x < wb; x++) { z = (x - pl + w2) % (w * 2); xs[x] = z < w ? z : w * 2 - z - 1; }
                    for (y = 0; y < hb; y++) { z = (y - pt + h2) % (h * 2); ys[y] = z < h ? z : h * 2 - z - 1; }
                }
                else if (flag == 3) {
                    for (x = 0; x < wb; x++) xs[x] = (x - pl + w2) % w;
                    for (y = 0; y < hb; y++) ys[y] = (y - pt + h2) % h;
                }
            }
            // pad by appropriate value
            for (z = 0; z < d; z++) {
                // copy over A to relevant region in B
                for (x = 0; x < w - cr - cl; x++)
                    memcpy(B + (x + pl)*hb + pt, A + (x + cl)*h + ct, sizeof(T)*(h - ct - cb));
                // set boundaries of B to appropriate values
                if (flag == 0 && val != 0) { // "constant"
                    for (x = 0; x < pl; x++) for (y = 0; y < hb; y++) B[x*hb + y] = val;
                    for (x = pl; x < w1; x++) for (y = 0; y < pt; y++) B[x*hb + y] = val;
                    for (x = pl; x < w1; x++) for (y = h1; y < hb; y++) B[x*hb + y] = val;
                    for (x = w1; x < wb; x++) for (y = 0; y < hb; y++) B[x*hb + y] = val;
                }
                else if (useLookup) { // "lookup"
                    PAD(xs[x], xs[x], xs[x], ys[y], ys[y], ys[y]);
                }
                else if (flag == 1) {  // "replicate"
                    PAD(0, x - pl, w - 1, 0, y - pt, h - 1);
                }
                else if (flag == 2) { // "symmetric"
                    PAD(pl - x - 1, x - pl, w + w1 - 1 - x, pt - y - 1, y - pt, h + h1 - 1 - y);
                }
                else if (flag == 3) { // "circular"
                    PAD(x - pl + w, x - pl, x - pl - w, y - pt + h, y - pt, y - pt - h);
                }
                A += h*w;  B += hb*wb;
            }
            if (useLookup) { free(xs); free(ys); }
#undef PAD
        }

       static HRESULT rgbConvert(const MexImageByte& imgSrc,
            const wchar_t* colorSpace,
            MexImageFloat& imgDst)
       {
           wchar_t* validCS[] = {L"gray", L"rgb", L"luv", L"hsv", L"orig"};
           
           int flag = -1;
           for (int i = 0; i < ARRAYSIZE(validCS); ++i)
           {
               if (wcscmp(colorSpace, validCS[i]) == 0)
               {
                   flag = i;
                   break;
               }
           }    // i

           const int width = imgSrc.width;
           const int height = imgSrc.height;
           const int channels = imgSrc.channels;

           assert(flag >= 0);

           imgDst.Initialize(width, height, channels);

           rgbConvert(&imgSrc.buffer[0], width * height, channels, flag, 1.f / 255, &imgDst.buffer[0]);

           return S_OK;
       }

       template<typename T>
       static HRESULT imResample(MexImage<T>& imgSrc,           
           const double scale,
           MexImage<T>& imgDst)
       {
           if (IS_EQUAL(scale, 1))
           {
               imgDst = imgSrc;

               return S_OK;
           }
           
           const int dstWidth = int(imgSrc.width * scale + 0.5);
           const int dstHeight = int(imgSrc.height * scale + 0.5);
           const int channels = imgSrc.channels;

           imgDst.Initialize(dstWidth, dstHeight, channels);

           bool isBilinear = true;  // TODO: handle nearest neighbor case as in matlab code

           if (isBilinear)
           {
               // TODO: deal with the (T == byte) case in matlab
               resample(&imgSrc.buffer[0], &imgDst.buffer[0], imgSrc.height, dstHeight, imgSrc.width, dstWidth, channels, 1.f);
           }
           else
           {
               // TODO: handle nearest neighbor resampling
           }

           return S_OK;
       }

       static HRESULT convTri(const MexImageFloat& imgSrc,           
           const double r,
           MexImageFloat& imgDst,           
           double s = 1,
           double nomex = 0)
       {
           if (IS_ZERO(r) && IS_EQUAL(s, 1))
           {
               imgDst = imgSrc;

               return S_OK;
           }

           const int m = __min(imgSrc.width, imgSrc.height);

           if (m < 4 || 2 * r + 1 >= m)
           {
               nomex = 1;
           }

           if (IS_ZERO(nomex))
           {
               if (r > 0 && r <= 1 && s <= 2)
               {
                   convConst(imgSrc, L"convTri1", int(12 / (r * (r + 2)) - 2), int(s), imgDst);
               }
               else
               {
                   convConst(imgSrc, L"convTri", int(r), int(s), imgDst);
               }
           }
           else
           {
                // TODO: translate the matlab code in here (currently OK as this branch is never reached)
           }

           return S_OK;
       }

       static HRESULT gradientMag(const MexImageFloat& imgSrc,
           MexImageFloat& M,
           MexImageFloat* pO = nullptr,
           double c = 0,
           double normRad = 0,
           double normConst = 0.005,           
           double full = 0)
       {
           if (IS_ZERO(normRad))
           {
               return S_OK;
           }

           gradientMagImpl(imgSrc, (int)c, (int)full, M, pO);

           MexImageFloat S;
           convTri(M, normRad, S);

           // operate on M
           gradientMagImpl(S, normConst, M);

           return S_OK;
       }       

       static HRESULT gradientHist(const MexImageFloat& M,
           const MexImageFloat& O,
           MexImageFloat& H,
           int binSize = 8,
           int nOrients = 9,
           int softBin = 1,
           int useHog = 0,
           float clipHog = 0.2f,
           bool full = false)
       {
           const int width = M.width;
           const int height = M.height;
           //const int channels = M.channels;

           const int hb = height / binSize;
           const int wb = width / binSize;
           const int nChns = useHog == 0 ? nOrients : (useHog == 1 ? nOrients * 4 : nOrients * 3 + 5);

           H.Initialize(wb, hb, nChns, true);   // must init to 0!           

           if (nOrients == 0) 
           {
               return S_OK;
           }          

           if (useHog == 0) 
           {
               gradHist((float*)&M.buffer[0], (float*)&O.buffer[0], &H.buffer[0], height, width, binSize, nOrients, softBin, full);
           }
           else if (useHog == 1) 
           {
               hog((float*)&M.buffer[0], (float*)&O.buffer[0], &H.buffer[0], height, width, binSize, nOrients, softBin, full, clipHog);
           }
           else 
           {
               fhog((float*)&M.buffer[0], (float*)&O.buffer[0], &H.buffer[0], height, width, binSize, nOrients, softBin, clipHog);
           }

           return S_OK;
       }

       static HRESULT gradient2(const MexImageFloat& imgSrc,
           MexImageFloat& Gx,
           MexImageFloat& Gy)
       {
           const int h = imgSrc.height;
           const int w = imgSrc.width;
           const int d = imgSrc.channels;

           Gx.Initialize(w, h, d);
           Gy.Initialize(w, h, d);

           grad2((float*)&imgSrc.buffer[0], &Gx.buffer[0], &Gy.buffer[0], h, w, d);

           return S_OK;
       }

       static HRESULT edgesNms(MexImageFloat& imgE,
           const MexImageFloat& imgO,
           int r,
           int s,
           float m,
           int nThreads)
       {
           MexImageFloat imgDst;

           imgDst.Initialize(imgE.width, imgE.height, imgE.channels);

           float* E0 = &imgE.buffer[0];
           float* O = (float*)&imgO.buffer[0];

           int h = imgE.height;
           int w = imgE.width;

           float* E = &imgDst.buffer[0];

           // suppress edges where edge is stronger in orthogonal direction
#ifdef USEOMP
           nThreads = nThreads<omp_get_max_threads() ? nThreads : omp_get_max_threads();
#pragma omp parallel for num_threads(nThreads)
#else
           UNREFERENCED_PARAMETER(nThreads);
#endif
           for (int x = 0; x < w; x++) for (int y = 0; y<h; y++) {
               float e = E[x*h + y] = E0[x*h + y]; if (!e) continue; e *= m;
               float coso = cos(O[x*h + y]), sino = sin(O[x*h + y]);
               for (int d = -r; d <= r; d++) if (d) {
                   float e0 = interp(E0, h, w, x + d*coso, y + d*sino);
                   if (e < e0) { E[x*h + y] = 0; break; }
               }
           }

           // suppress noisy edge estimates near boundaries
           s = s > w / 2 ? w / 2 : s; s>h / 2 ? h / 2 : s;
           for (int x = 0; x < s; x++) for (int y = 0; y < h; y++) {
               E[x*h + y] *= x / float(s); E[(w - 1 - x)*h + y] *= x / float(s);
           }
           for (int x = 0; x < w; x++) for (int y = 0; y < s; y++) {
               E[x*h + y] *= y / float(s); E[x*h + (h - 1 - y)] *= y / float(s);
           }

           imgE.buffer.swap(imgDst.buffer);

           return S_OK;
       }

#pragma region --- lookup for edge detection ---
       // construct lookup array for mapping fids to channel indices
       static UINT32* buildLookup(int *dims, int w) {
           int c, r, z, n = w*w*dims[2]; UINT32 *cids = new UINT32[n]; n = 0;
           for (z = 0; z < dims[2]; z++) for (c = 0; c < w; c++) for (r = 0; r < w; r++)
               cids[n++] = z*dims[0] * dims[1] + c*dims[0] + r;
           return cids;
       }

       // construct lookup arrays for mapping fids for self-similarity channel
       static void buildLookupSs(UINT32 *&cids1, UINT32 *&cids2, int *dims, int w, int m) {
           int i, j, z, z1, c, r; int locs[1024];
           int m2 = m*m, n = m2*(m2 - 1) / 2 * dims[2], s = int(w / m / 2.0 + .5);
           cids1 = new UINT32[n]; cids2 = new UINT32[n]; n = 0;
           for (i = 0; i < m; i++) locs[i] = UINT32((i + 1)*(w + 2 * s - 1) / (m + 1.0) - s + .5);
           for (z = 0; z < dims[2]; z++) for (i = 0; i < m2; i++) for (j = i + 1; j < m2; j++) {
               z1 = z*dims[0] * dims[1]; n++;
               r = i%m; c = (i - r) / m; cids1[n - 1] = z1 + locs[c] * dims[0] + locs[r];
               r = j%m; c = (j - r) / m; cids2[n - 1] = z1 + locs[c] * dims[0] + locs[r];
           }
       }
#pragma endregion

    private:
#pragma region --- color conversion ---
        // Copy and normalize only
        template<class iT, class oT> 
        static void normalize(iT *I, oT *J, int n, oT nrm) {
            for (int i = 0; i<n; i++) *(J++) = (oT)*(I++)*nrm;
        }

        template<class iT, class oT>
        static void rgbConvert(iT *I, int n, int d, int flag, oT nrm, oT* J) 
        {        
            int i, n1 = d*(n<1000 ? n / 10 : 100); oT thr = oT(1.001);
            if (flag>1 && nrm == 1)
            {
                for (i = 0; i<n1; i++) 
                {
                    if (I[i]>thr)
                    {
                        //wrError("For floats all values in I must be smaller than 1.");

                        return;
                    }
                }            
            }        
            
            bool useSse = false;    //n % 4 == 0 && typeid(oT) == typeid(float); TODO: translate SSE code
            if (flag == 2 && useSse)
            {
                //for (i = 0; i<d / 3; i++) rgb2luv_sse(I + i*n * 3, (float*)(J + i*n * 3), n, (float)nrm);
            }            
            else if ((flag == 0 && d == 1) || flag == 1)
            {
                normalize(I, J, n*d, nrm);
            }
            else if (flag == 0)
            {
                // TODO:
                //for (i = 0; i<d / 3; i++) rgb2gray(I + i*n * 3, J + i*n * 1, n, nrm);
            }            
            else if (flag == 2)
            {
                for (i = 0; i < d / 3; i++) rgb2luv(I + i*n * 3, J + i*n * 3, n, nrm);
            }
            else if (flag == 3) 
            {
                // TODO:
                //for (i = 0; i<d / 3; i++) rgb2hsv(I + i*n * 3, J + i*n * 3, n, nrm);
            }            
            else
            { 
                //wrError("Unknown flag.");
            }
            
            //return J;
        }

        template<class iT, class oT> 
        static void rgb2luv(iT *I, oT *J, int n, oT nrm) {
            oT minu, minv, un, vn, mr[3], mg[3], mb[3];
            oT *lTable = rgb2luv_setup(nrm, mr, mg, mb, minu, minv, un, vn);
            oT *L = J, *U = L + n, *V = U + n; iT *R = I, *G = R + n, *B = G + n;
            for (int i = 0; i<n; i++) {
                oT r, g, b, x, y, z, l;
                r = (oT)*R++; g = (oT)*G++; b = (oT)*B++;
                x = mr[0] * r + mg[0] * g + mb[0] * b;
                y = mr[1] * r + mg[1] * g + mb[1] * b;
                z = mr[2] * r + mg[2] * g + mb[2] * b;
                l = lTable[(int)(y * 1024)];
                *(L++) = l; z = 1 / (x + 15 * y + 3 * z + (oT)1e-35);
                *(U++) = l * (13 * 4 * x*z - 13 * un) - minu;
                *(V++) = l * (13 * 9 * y*z - 13 * vn) - minv;
            }
        }

        template<class oT> 
        static oT* rgb2luv_setup(oT z, oT *mr, oT *mg, oT *mb,
            oT &minu, oT &minv, oT &un, oT &vn)
        {
            // set constants for conversion
            const oT y0 = (oT)((6.0 / 29)*(6.0 / 29)*(6.0 / 29));
            const oT a = (oT)((29.0 / 3)*(29.0 / 3)*(29.0 / 3));
            un = (oT) 0.197833; vn = (oT) 0.468331;
            mr[0] = (oT) 0.430574*z; mr[1] = (oT) 0.222015*z; mr[2] = (oT) 0.020183*z;
            mg[0] = (oT) 0.341550*z; mg[1] = (oT) 0.706655*z; mg[2] = (oT) 0.129553*z;
            mb[0] = (oT) 0.178325*z; mb[1] = (oT) 0.071330*z; mb[2] = (oT) 0.939180*z;
            oT maxi = (oT) 1.0 / 270; minu = -88 * maxi; minv = -134 * maxi;
            // build (padded) lookup table for y->l conversion assuming y in [0,1]
            static oT lTable[1064]; static bool lInit = false;
            if (lInit) return lTable; oT y, l;
            for (int i = 0; i<1025; i++) {
                y = (oT)(i / 1024.0);
                l = y>y0 ? 116 * (oT)pow((double)y, 1.0 / 3.0) - 16 : y*a;
                lTable[i] = l*maxi;
            }
            for (int i = 1025; i<1064; i++) lTable[i] = lTable[i - 1];
            lInit = true; return lTable;    // TODO: bad practice for returning address of static array
        }
#pragma endregion

#pragma region --- image sampling ---
        // compute interpolation values for single column for resapling
        template<class T> 
        static void resampleCoef(int ha, int hb, int &n, int *&yas,
            int *&ybs, T *&wts, int bd[2], int pad = 0)
        {
            const T s = T(hb) / T(ha), sInv = 1 / s; T wt, wt0 = T(1e-3)*s;
            bool ds = ha>hb; int nMax; bd[0] = bd[1] = 0;
            if (ds) { n = 0; nMax = ha + (pad>2 ? pad : 2)*hb; }
            else { n = nMax = hb; }
            // initialize memory
            wts = (T*)alMalloc(nMax*sizeof(T), 16);
            yas = (int*)alMalloc(nMax*sizeof(int), 16);
            ybs = (int*)alMalloc(nMax*sizeof(int), 16);
            if (ds) for (int yb = 0; yb<hb; yb++) {
                // create coefficients for downsampling
                T ya0f = yb*sInv, ya1f = ya0f + sInv, W = 0;
                int ya0 = int(ceil(ya0f)), ya1 = int(ya1f), n1 = 0;
                for (int ya = ya0 - 1; ya<ya1 + 1; ya++) {
                    wt = s; if (ya == ya0 - 1) wt = (ya0 - ya0f)*s; else if (ya == ya1) wt = (ya1f - ya1)*s;
                    if (wt>wt0 && ya >= 0) { ybs[n] = yb; yas[n] = ya; wts[n] = wt; n++; n1++; W += wt; }
                }
                if (W>1) for (int i = 0; i<n1; i++) wts[n - n1 + i] /= W;
                if (n1>bd[0]) bd[0] = n1;
                while (n1<pad) { ybs[n] = yb; yas[n] = yas[n - 1]; wts[n] = 0; n++; n1++; }
            }
            else for (int yb = 0; yb<hb; yb++) {
                // create coefficients for upsampling
                T yaf = (T(.5) + yb)*sInv - T(.5); int ya = (int)floor(yaf);
                wt = 1; if (ya >= 0 && ya<ha - 1) wt = 1 - (yaf - ya);
                if (ya<0) { ya = 0; bd[0]++; } if (ya >= ha - 1) { ya = ha - 1; bd[1]++; }
                ybs[yb] = yb; yas[yb] = ya; wts[yb] = wt;
            }
        }

        // resample A using bilinear interpolation and and store result in B
        template<class T>
        static void resample(T *A, T *B, int ha, int hb, int wa, int wb, int d, T r) {
            int hn, wn, x, x1, y, z, xa, xb, ya; T *A0, *A1, *A2, *A3, *B0, wt, wt1;
            T *C = (T*)alMalloc((ha + 4)*sizeof(T), 16); for (y = ha; y < ha + 4; y++) C[y] = 0;
            bool sse = (typeid(T) == typeid(float)) && !(size_t(A) & 15) && !(size_t(B) & 15);
            // get coefficients for resampling along w and h
            int *xas, *xbs, *yas, *ybs; T *xwts, *ywts; int xbd[2], ybd[2];
            resampleCoef<T>(wa, wb, wn, xas, xbs, xwts, xbd, 0);
            resampleCoef<T>(ha, hb, hn, yas, ybs, ywts, ybd, 4);
            if (wa == 2 * wb) r /= 2; if (wa == 3 * wb) r /= 3; if (wa == 4 * wb) r /= 4;
            r /= T(1 + 1e-6); for (y = 0; y < hn; y++) ywts[y] *= r;
            // resample each channel in turn
            for (z = 0; z < d; z++) for (x = 0; x < wb; x++) {
                if (x == 0) x1 = 0; xa = xas[x1]; xb = xbs[x1]; wt = xwts[x1]; wt1 = 1 - wt; y = 0;
                A0 = A + z*ha*wa + xa*ha; A1 = A0 + ha, A2 = A1 + ha, A3 = A2 + ha; B0 = B + z*hb*wb + xb*hb;
                // variables for SSE (simple casts to float)
                float *Af0, *Af1, *Af2, *Af3, *Bf0, *Cf, *ywtsf, wtf, wt1f;
                Af0 = (float*)A0; Af1 = (float*)A1; Af2 = (float*)A2; Af3 = (float*)A3;
                Bf0 = (float*)B0; Cf = (float*)C;
                ywtsf = (float*)ywts; wtf = (float)wt; wt1f = (float)wt1;
                // resample along x direction (A -> C)
#define FORs(X) if(sse) for(; y<ha-4; y+=4) STR(Cf[y],X);
#define FORr(X) for(; y<ha; y++) C[y] = X;
                if (wa == 2 * wb) {
                    FORs(ADD(LDu(Af0[y]), LDu(Af1[y])));
                    FORr(A0[y] + A1[y]); x1 += 2;
                }
                else if (wa == 3 * wb) {
                    FORs(ADD(LDu(Af0[y]), LDu(Af1[y]), LDu(Af2[y])));
                    FORr(A0[y] + A1[y] + A2[y]); x1 += 3;
                }
                else if (wa == 4 * wb) {
                    FORs(ADD(LDu(Af0[y]), LDu(Af1[y]), LDu(Af2[y]), LDu(Af3[y])));
                    FORr(A0[y] + A1[y] + A2[y] + A3[y]); x1 += 4;
                }
                else if (wa > wb) {
                    int m = 1; while (x1 + m < wn && xb == xbs[x1 + m]) m++; float wtsf[4];
                    for (int x0 = 0; x0 < (m < 4 ? m : 4); x0++) wtsf[x0] = float(xwts[x1 + x0]);
#define U(x) MUL( LDu(*(Af ## x + y)), SET(wtsf[x]) )
#define V(x) *(A ## x + y) * xwts[x1+x]
                    if (m == 1) { FORs(U(0));                     FORr(V(0)); }
                    if (m == 2) { FORs(ADD(U(0), U(1)));           FORr(V(0) + V(1)); }
                    if (m == 3) { FORs(ADD(U(0), U(1), U(2)));      FORr(V(0) + V(1) + V(2)); }
                    if (m >= 4) { FORs(ADD(U(0), U(1), U(2), U(3))); FORr(V(0) + V(1) + V(2) + V(3)); }
#undef U
#undef V
                    for (int x0 = 4; x0 < m; x0++) {
                        A1 = A0 + x0*ha; wt1 = xwts[x1 + x0]; Af1 = (float*)A1; wt1f = float(wt1); y = 0;
                        FORs(ADD(LD(Cf[y]), MUL(LDu(Af1[y]), SET(wt1f)))); FORr(C[y] + A1[y] * wt1);
                    }
                    x1 += m;
                }
                else {
                    bool xBd = x < xbd[0] || x >= wb - xbd[1]; x1++;
                    if (xBd) memcpy(C, A0, ha*sizeof(T));
                    if (!xBd) FORs(ADD(MUL(LDu(Af0[y]), SET(wtf)), MUL(LDu(Af1[y]), SET(wt1f))));
                    if (!xBd) FORr(A0[y] * wt + A1[y] * wt1);
                }
#undef FORs
#undef FORr
                // resample along y direction (B -> C)
                if (ha == hb * 2) {
                    T r2 = r / 2; int k = ((~((size_t)B0) + 1) & 15) / 4; y = 0;
                    for (; y < k; y++)  B0[y] = (C[2 * y] + C[2 * y + 1])*r2;
                    if (sse) for (; y < hb - 4; y += 4) STR(Bf0[y], MUL((float)r2, _mm_shuffle_ps(ADD(
                        LDu(Cf[2 * y]), LDu(Cf[2 * y + 1])), ADD(LDu(Cf[2 * y + 4]), LDu(Cf[2 * y + 5])), 136)));
                    for (; y < hb; y++) B0[y] = (C[2 * y] + C[2 * y + 1])*r2;
                }
                else if (ha == hb * 3) {
                    for (y = 0; y < hb; y++) B0[y] = (C[3 * y] + C[3 * y + 1] + C[3 * y + 2])*(r / 3);
                }
                else if (ha == hb * 4) {
                    for (y = 0; y < hb; y++) B0[y] = (C[4 * y] + C[4 * y + 1] + C[4 * y + 2] + C[4 * y + 3])*(r / 4);
                }
                else if (ha > hb) {
                    y = 0;
                    //if( sse && ybd[0]<=4 ) for(; y<hb; y++) // Requires SSE4
                    //  STR1(Bf0[y],_mm_dp_ps(LDu(Cf[yas[y*4]]),LDu(ywtsf[y*4]),0xF1));
#define U(o) C[ya+o]*ywts[y*4+o]
                    if (ybd[0] == 2) for (; y < hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1); }
                    if (ybd[0] == 3) for (; y < hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1) + U(2); }
                    if (ybd[0] == 4) for (; y < hb; y++) { ya = yas[y * 4]; B0[y] = U(0) + U(1) + U(2) + U(3); }
                    if (ybd[0] > 4)  for (; y < hn; y++) { B0[ybs[y]] += C[yas[y]] * ywts[y]; }
#undef U
                }
                else {
                    for (y = 0; y < ybd[0]; y++) B0[y] = C[yas[y]] * ywts[y];
                    for (; y < hb - ybd[1]; y++) B0[y] = C[yas[y]] * ywts[y] + C[yas[y] + 1] * (r - ywts[y]);
                    for (; y < hb; y++)        B0[y] = C[yas[y]] * ywts[y];
                }
            }
            alFree(xas); alFree(xbs); alFree(xwts); alFree(C);
            alFree(yas); alFree(ybs); alFree(ywts);
        }
#pragma endregion

#pragma region --- convolution ---
        static void convConst(const MexImageFloat& imgSrc,            
            const wchar_t* convType,
            const int r,
            const int s,
            MexImageFloat& imgDst)
        {
            const int m = __min(imgSrc.width, imgSrc.height);            

            imgDst.Initialize(imgSrc.width / s, imgSrc.height / s, imgSrc.channels);

            int ns[] = { imgSrc.height, imgSrc.width };
            float* A = (float*)&imgSrc.buffer[0];
            float* B = &imgDst.buffer[0];
            const int d = imgDst.channels;
            const float p = (float)r;

            if (!wcscmp(convType, L"convBox")) {
                if (r >= m / 2) mexErrMsgTxt("mask larger than image (r too large)");
                convBox(A, B, ns[0], ns[1], d, r, s);
            }
            else if (!wcscmp(convType, L"convTri")) {
                if (r >= m / 2) mexErrMsgTxt("mask larger than image (r too large)");
                convTri(A, B, ns[0], ns[1], d, r, s);
            }
            else if (!wcscmp(convType, L"conv11")) {
                if (s > 2) mexErrMsgTxt("conv11 can sample by at most s=2");
                conv11(A, B, ns[0], ns[1], d, r, s);
            }
            else if (!wcscmp(convType, L"convTri1")) {
                if (s > 2) mexErrMsgTxt("convTri1 can sample by at most s=2");
                convTri1(A, B, ns[0], ns[1], d, p, s);
            }
            else if (!wcscmp(convType, L"convMax")) {
                if (s > 1) mexErrMsgTxt("convMax cannot sample");
                convMax(A, B, ns[0], ns[1], d, r);
            }
            else 
            {
                mexErrMsgTxt("Invalid type.");
            }
        }

        // convolve one column of I by a 2rx1 ones filter
        static void convBoxY(float *I, float *O, int h, int r, int s) {
            float t; int j, p = r + 1, q = 2 * h - (r + 1), h0 = r + 1, h1 = h - r, h2 = h;
            t = 0; for (j = 0; j <= r; j++) t += I[j]; t = 2 * t - I[r]; j = 0;
            if (s == 1) {
                for (; j < h0; j++) O[j] = t -= I[r - j] - I[r + j];
                for (; j < h1; j++) O[j] = t -= I[j - p] - I[r + j];
                for (; j < h2; j++) O[j] = t -= I[j - p] - I[q - j];
            }
            else {
                int k = (s - 1) / 2; h2 = (h / s)*s; if (h0 > h2) h0 = h2; if (h1 > h2) h1 = h2;
                for (; j < h0; j++) { t -= I[r - j] - I[r + j]; k++; if (k == s) { k = 0; *O++ = t; } }
                for (; j < h1; j++) { t -= I[j - p] - I[r + j]; k++; if (k == s) { k = 0; *O++ = t; } }
                for (; j < h2; j++) { t -= I[j - p] - I[q - j]; k++; if (k == s) { k = 0; *O++ = t; } }
            }
        }

        // convolve I by a 2r+1 x 2r+1 ones filter (uses SSE)
        static void convBox(float *I, float *O, int h, int w, int d, int r, int s) {
            float nrm = 1.0f / ((2 * r + 1)*(2 * r + 1)); int i, j, k = (s - 1) / 2, h0, h1, w0;
            if (h % 4 == 0) h0 = h1 = h; else { h0 = h - (h % 4); h1 = h0 + 4; } w0 = (w / s)*s;
            float *T = (float*)alMalloc(h1*sizeof(float), 16);
            while (d-- > 0) {
                // initialize T
                memset(T, 0, h1*sizeof(float));
                for (i = 0; i <= r; i++) for (j = 0; j < h0; j += 4) INC(T[j], LDu(I[j + i*h]));
                for (j = 0; j < h0; j += 4) STR(T[j], MUL(nrm, SUB(MUL(2, LD(T[j])), LDu(I[j + r*h]))));
                for (i = 0; i <= r; i++) for (j = h0; j < h; j++) T[j] += I[j + i*h];
                for (j = h0; j < h; j++) T[j] = nrm*(2 * T[j] - I[j + r*h]);
                // prepare and convolve each column in turn
                k++; if (k == s) { k = 0; convBoxY(T, O, h, r, s); O += h / s; }
                for (i = 1; i < w0; i++) {
                    float *Il = I + (i - 1 - r)*h; if (i <= r) Il = I + (r - i)*h;
                    float *Ir = I + (i + r)*h; if (i >= w - r) Ir = I + (2 * w - r - i - 1)*h;
                    for (j = 0; j < h0; j += 4) DEC(T[j], MUL(nrm, SUB(LDu(Il[j]), LDu(Ir[j]))));
                    for (j = h0; j < h; j++) T[j] -= nrm*(Il[j] - Ir[j]);
                    k++; if (k == s) { k = 0; convBoxY(T, O, h, r, s); O += h / s; }
                }
                I += w*h;
            }
            alFree(T);
        }

        // convolve one column of I by a [1; 1] filter (uses SSE)
        static void conv11Y(float *I, float *O, int h, int side, int s) {
#define C4(m,o) ADD(LDu(I[m*j-1+o]),LDu(I[m*j+o]))
            int j = 0, k = ((~((size_t)O) + 1) & 15) / 4;
            const int d = (side % 4 >= 2) ? 1 : 0, h2 = (h - d) / 2;
            if (s == 2) {
                for (; j < k; j++) O[j] = I[2 * j + d] + I[2 * j + d + 1];
                for (; j < h2 - 4; j += 4) STR(O[j], _mm_shuffle_ps(C4(2, d + 1), C4(2, d + 5), 136));
                for (; j < h2; j++) O[j] = I[2 * j + d] + I[2 * j + d + 1];
                if (d == 1 && h % 2 == 0) O[j] = 2 * I[2 * j + d];
            }
            else {
                if (d == 0) { O[0] = 2 * I[0]; j++; if (k == 0) k = 4; }
                for (; j < k; j++) O[j] = I[j - 1 + d] + I[j + d];
                for (; j < h - 4 - d; j += 4) STR(O[j], C4(1, d));
                for (; j < h - d; j++) O[j] = I[j - 1 + d] + I[j + d];
                if (d == 1) { O[j] = 2 * I[j]; j++; }
            }
#undef C4
        }

        // convolve I by a [1 1; 1 1] filter (uses SSE)
        static void conv11(float *I, float *O, int h, int w, int d, int side, int s) {
            const float nrm = 0.25f; int i, j;
            float *I0, *I1, *T = (float*)alMalloc(h*sizeof(float), 16);
            for (int d0 = 0; d0 < d; d0++) for (i = s / 2; i < w; i += s) {
                I0 = I1 = I + i*h + d0*h*w; if (side % 2) { if (i < w - 1) I1 += h; }
                else { if (i) I0 -= h; }
                for (j = 0; j < h - 4; j += 4) STR(T[j], MUL(nrm, ADD(LDu(I0[j]), LDu(I1[j]))));
                for (; j < h; j++) T[j] = nrm*(I0[j] + I1[j]);
                conv11Y(T, O, h, side, s); O += h / s;
            }
            alFree(T);
        }

        // convolve one column of I by a 2rx1 triangle filter
        static void convTriY(float *I, float *O, int h, int r, int s) {
            r++; float t, u; int j, r0 = r - 1, r1 = r + 1, r2 = 2 * h - r, h0 = r + 1, h1 = h - r + 1, h2 = h;
            u = t = I[0]; for (j = 1; j < r; j++) u += t += I[j]; u = 2 * u - t; t = 0;
            if (s == 1) {
                O[0] = u; j = 1;
                for (; j < h0; j++) O[j] = u += t += I[r - j] + I[r0 + j] - 2 * I[j - 1];
                for (; j < h1; j++) O[j] = u += t += I[j - r1] + I[r0 + j] - 2 * I[j - 1];
                for (; j < h2; j++) O[j] = u += t += I[j - r1] + I[r2 - j] - 2 * I[j - 1];
            }
            else {
                int k = (s - 1) / 2; h2 = (h / s)*s; if (h0 > h2) h0 = h2; if (h1 > h2) h1 = h2;
                if (++k == s) { k = 0; *O++ = u; } j = 1;
                for (; j < h0; j++) { u += t += I[r - j] + I[r0 + j] - 2 * I[j - 1]; if (++k == s){ k = 0; *O++ = u; } }
                for (; j < h1; j++) { u += t += I[j - r1] + I[r0 + j] - 2 * I[j - 1]; if (++k == s){ k = 0; *O++ = u; } }
                for (; j < h2; j++) { u += t += I[j - r1] + I[r2 - j] - 2 * I[j - 1]; if (++k == s){ k = 0; *O++ = u; } }
            }
        }

        // convolve I by a 2rx1 triangle filter (uses SSE)
        static void convTri(float *I, float *O, int h, int w, int d, int r, int s) {
            r++; float nrm = 1.0f / (r*r*r*r); int i, j, k = (s - 1) / 2, h0, h1, w0;
            if (h % 4 == 0) h0 = h1 = h; else { h0 = h - (h % 4); h1 = h0 + 4; } w0 = (w / s)*s;
            float *T = (float*)alMalloc(2 * h1*sizeof(float), 16), *U = T + h1;
            while (d-- > 0) {
                // initialize T and U
                for (j = 0; j < h0; j += 4) STR(U[j], STR(T[j], LDu(I[j])));
                for (i = 1; i < r; i++) for (j = 0; j < h0; j += 4) INC(U[j], INC(T[j], LDu(I[j + i*h])));
                for (j = 0; j < h0; j += 4) STR(U[j], MUL(nrm, (SUB(MUL(2, LD(U[j])), LD(T[j])))));
                for (j = 0; j < h0; j += 4) STR(T[j], 0);
                for (j = h0; j < h; j++) U[j] = T[j] = I[j];
                for (i = 1; i < r; i++) for (j = h0; j < h; j++) U[j] += T[j] += I[j + i*h];
                for (j = h0; j < h; j++) { U[j] = nrm * (2 * U[j] - T[j]); T[j] = 0; }
                // prepare and convolve each column in turn
                k++; if (k == s) { k = 0; convTriY(U, O, h, r - 1, s); O += h / s; }
                for (i = 1; i < w0; i++) {
                    float *Il = I + (i - 1 - r)*h; if (i <= r) Il = I + (r - i)*h; float *Im = I + (i - 1)*h;
                    float *Ir = I + (i - 1 + r)*h; if (i > w - r) Ir = I + (2 * w - r - i)*h;
                    for (j = 0; j < h0; j += 4) {
                        INC(T[j], ADD(LDu(Il[j]), LDu(Ir[j]), MUL(-2, LDu(Im[j]))));
                        INC(U[j], MUL(nrm, LD(T[j])));
                    }
                    for (j = h0; j < h; j++) U[j] += nrm*(T[j] += Il[j] + Ir[j] - 2 * Im[j]);
                    k++; if (k == s) { k = 0; convTriY(U, O, h, r - 1, s); O += h / s; }
                }
                I += w*h;
            }
            alFree(T);
        }

        // convolve one column of I by a [1 p 1] filter (uses SSE)
        static void convTri1Y(float *I, float *O, int h, float p, int s) {
#define C4(m,o) ADD(ADD(LDu(I[m*j-1+o]),MUL(p,LDu(I[m*j+o]))),LDu(I[m*j+1+o]))
            int j = 0, k = ((~((size_t)O) + 1) & 15) / 4, h2 = (h - 1) / 2;
            if (s == 2) {
                for (; j < k; j++) O[j] = I[2 * j] + p*I[2 * j + 1] + I[2 * j + 2];
                for (; j < h2 - 4; j += 4) STR(O[j], _mm_shuffle_ps(C4(2, 1), C4(2, 5), 136));
                for (; j < h2; j++) O[j] = I[2 * j] + p*I[2 * j + 1] + I[2 * j + 2];
                if (h % 2 == 0) O[j] = I[2 * j] + (1 + p)*I[2 * j + 1];
            }
            else {
                O[j] = (1 + p)*I[j] + I[j + 1]; j++; if (k == 0) k = (h <= 4) ? h - 1 : 4;
                for (; j < k; j++) O[j] = I[j - 1] + p*I[j] + I[j + 1];
                for (; j < h - 4; j += 4) STR(O[j], C4(1, 0));
                for (; j < h - 1; j++) O[j] = I[j - 1] + p*I[j] + I[j + 1];
                O[j] = I[j - 1] + (1 + p)*I[j];
            }
#undef C4
        }

        // convolve I by a [1 p 1] filter (uses SSE)
        static void convTri1(float *I, float *O, int h, int w, int d, float p, int s) {
            const float nrm = 1.0f / ((p + 2)*(p + 2)); int i, j, h0 = h - (h % 4);
            float *Il, *Im, *Ir, *T = (float*)alMalloc(h*sizeof(float), 16);
            for (int d0 = 0; d0 < d; d0++) for (i = s / 2; i < w; i += s) {
                Il = Im = Ir = I + i*h + d0*h*w; if (i>0) Il -= h; if (i < w - 1) Ir += h;
                for (j = 0; j < h0; j += 4)
                    STR(T[j], MUL(nrm, ADD(ADD(LDu(Il[j]), MUL(p, LDu(Im[j]))), LDu(Ir[j]))));
                for (j = h0; j < h; j++) T[j] = nrm*(Il[j] + p*Im[j] + Ir[j]);
                convTri1Y(T, O, h, p, s); O += h / s;
            }
            alFree(T);
        }

        // convolve one column of I by a 2rx1 max filter
        static void convMaxY(float *I, float *O, float *T, int h, int r) {
            int y, y0, y1, yi, m = 2 * r + 1;
#define max1(a,b) a>b ? a : b;
#define maxk(y0,y1) { O[y]=I[y0]; \
    for( yi=y0+1; yi<=y1; yi++ ) { if(I[yi]>O[y]) O[y]=I[yi]; }}
            for (y = 0; y<r; y++) { y1 = y + r; if (y1>h - 1) y1 = h - 1; maxk(0, y1); }
            for (; y <= h - m - r; y += m) {
                T[m - 1] = I[y + r];
                for (yi = 1; yi < m; yi++) T[m - 1 - yi] = max1(T[m - 1 - yi + 1], I[y + r - yi]);
                for (yi = 1; yi < m; yi++) T[m - 1 + yi] = max1(T[m - 1 + yi - 1], I[y + r + yi]);
                for (yi = 0; yi < m; yi++) O[y + yi] = max1(T[yi], T[yi + m - 1]);
            }
            for (; y < h - r; y++) { maxk(y - r, y + r); }
            for (; y < h; y++) { y0 = y - r; if (y0 < 0) y0 = 0; maxk(y0, h - 1); }
#undef maxk
#undef max1
        }

        // convolve I by a 2rx1 max filter
        static void convMax(float *I, float *O, int h, int w, int d, int r) {
            if (r > w - 1) r = w - 1; if (r > h - 1) r = h - 1; int m = 2 * r + 1;
            float *T = (float*)alMalloc(m * 2 * sizeof(float), 16);
            for (int d0 = 0; d0 < d; d0++) for (int x = 0; x < w; x++) {
                float *Oc = O + d0*h*w + h*x, *Ic = I + d0*h*w + h*x;
                convMaxY(Ic, Oc, T, h, r);
            }
            alFree(T);
        }
#pragma endregion

#pragma region --- graidents ---
        // "gradientMag" case
        static void gradientMagImpl(const MexImageFloat& imgSrc,
            const int c,
            const int full,
            MexImageFloat& M, 
            MexImageFloat* pO)
        {
            const int height = imgSrc.height;
            const int width = imgSrc.width;
            const int channels = imgSrc.channels;

            M.Initialize(width, height, 1);

            if (pO != nullptr)
            {
                pO->Initialize(width, height, 1);
            }

            if (c > 0 && c <= channels)
            {
                MexImageFloat imgTemp;
                imgTemp.Initialize(width, height, 1);               

                memcpy(&imgTemp.buffer[0], &imgSrc.buffer[(c - 1) * width * height], sizeof(float) * width * height);

                gradMag(&imgTemp.buffer[0], &M.buffer[0], pO == nullptr ? nullptr : &pO->buffer[0], height, width, imgTemp.channels, full > 0);
            }
            else
            {
                gradMag((float*)&imgSrc.buffer[0], &M.buffer[0], pO == nullptr ? nullptr : &pO->buffer[0], height, width, imgSrc.channels, full > 0);
            }
        }

        // "gradientMagNorm" case
        static void gradientMagImpl(const MexImageFloat& S,
            const double normRad,
            MexImageFloat& M)
        {
            gradMagNorm(&M.buffer[0], (float*)&S.buffer[0], M.height, M.width, (float)normRad);
        }

        // compute x and y gradients for just one column (uses sse)
        static void grad1(float *I, float *Gx, float *Gy, int h, int w, int x) {
            int y, y1; float *Ip, *In, r; __m128 *_Ip, *_In, *_G, _r;
            // compute column of Gx
            Ip = I - h; In = I + h; r = .5f;
            if (x == 0) { r = 1; Ip += h; }
            else if (x == w - 1) { r = 1; In -= h; }
            if (h < 4 || h % 4>0 || (size_t(I) & 15) || (size_t(Gx) & 15)) {
                for (y = 0; y < h; y++) *Gx++ = (*In++ - *Ip++)*r;
            }
            else {
                _G = (__m128*) Gx; _Ip = (__m128*) Ip; _In = (__m128*) In; _r = SET(r);
                for (y = 0; y < h; y += 4) *_G++ = MUL(SUB(*_In++, *_Ip++), _r);
            }
            // compute column of Gy
#define GRADY(r) *Gy++=(*In++-*Ip++)*r;
            Ip = I; In = Ip + 1;
            // GRADY(1); Ip--; for(y=1; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
            y1 = ((~((size_t)Gy) + 1) & 15) / 4; if (y1 == 0) y1 = 4; if (y1 > h - 1) y1 = h - 1;
            GRADY(1); Ip--; for (y = 1; y < y1; y++) GRADY(.5f);
            _r = SET(.5f); _G = (__m128*) Gy;
            for (; y + 4 < h - 1; y += 4, Ip += 4, In += 4, Gy += 4)
                *_G++ = MUL(SUB(LDu(*In), LDu(*Ip)), _r);
            for (; y < h - 1; y++) GRADY(.5f); In--; GRADY(1);
#undef GRADY
        }

        // compute x and y gradients at each location (uses sse)
        static void grad2(float *I, float *Gx, float *Gy, int h, int w, int d) {
            int o, x, c, a = w*h; for (c = 0; c < d; c++) for (x = 0; x < w; x++) {
                o = c*a + x*h; grad1(I + o, Gx + o, Gy + o, h, w, x);
            }
        }

        // build lookup table a[] s.t. a[x*n]~=acos(x) for x in [-1,1]
        static float* acosTable() {
            const int n = 10000, b = 10; int i;
            static float a[n * 2 + b * 2]; static bool init = false;
            float *a1 = a + n + b; if (init) return a1;
            for (i = -n - b; i < -n; i++)   a1[i] = PI;
            for (i = -n; i < n; i++)      a1[i] = float(acos(i / float(n)));
            for (i = n; i < n + b; i++)     a1[i] = 0;
            for (i = -n - b; i<n / 10; i++) if (a1[i] > PI - 1e-6f) a1[i] = PI - 1e-6f;
            init = true; return a1;
        }

        // compute gradient magnitude and orientation at each location (uses sse)
        static void gradMag(float *I, float *M, float *O, int h, int w, int d, bool full) {
            int x, y, y1, c, h4, s; float *Gx, *Gy, *M2; __m128 *_Gx, *_Gy, *_M2, _m;
            float *acost = acosTable(), acMult = 10000.0f;
            // allocate memory for storing one column of output (padded so h4%4==0)
            h4 = (h % 4 == 0) ? h : h - (h % 4) + 4; s = d*h4*sizeof(float);
            M2 = (float*)alMalloc(s, 16); _M2 = (__m128*) M2;
            Gx = (float*)alMalloc(s, 16); _Gx = (__m128*) Gx;
            Gy = (float*)alMalloc(s, 16); _Gy = (__m128*) Gy;
            // compute gradient magnitude and orientation for each column
            for (x = 0; x < w; x++) {
                // compute gradients (Gx, Gy) with maximum squared magnitude (M2)
                for (c = 0; c < d; c++) {
                    grad1(I + x*h + c*w*h, Gx + c*h4, Gy + c*h4, h, w, x);
                    for (y = 0; y < h4 / 4; y++) {
                        y1 = h4 / 4 * c + y;
                        _M2[y1] = ADD(MUL(_Gx[y1], _Gx[y1]), MUL(_Gy[y1], _Gy[y1]));
                        if (c == 0) continue; _m = CMPGT(_M2[y1], _M2[y]);
                        _M2[y] = OR(AND(_m, _M2[y1]), ANDNOT(_m, _M2[y]));
                        _Gx[y] = OR(AND(_m, _Gx[y1]), ANDNOT(_m, _Gx[y]));
                        _Gy[y] = OR(AND(_m, _Gy[y1]), ANDNOT(_m, _Gy[y]));
                    }
                }
                // compute gradient mangitude (M) and normalize Gx
                for (y = 0; y < h4 / 4; y++) {
                    _m = MIN(RCPSQRT(_M2[y]), SET(1e10f));
                    _M2[y] = RCP(_m);
                    if (O) _Gx[y] = MUL(MUL(_Gx[y], _m), SET(acMult));
                    if (O) _Gx[y] = XOR(_Gx[y], AND(_Gy[y], SET(-0.f)));
                };
                memcpy(M + x*h, M2, h*sizeof(float));
                // compute and store gradient orientation (O) via table lookup
                if (O != 0) for (y = 0; y < h; y++) O[x*h + y] = acost[(int)Gx[y]];
                if (O != 0 && full) {
                    y1 = ((~size_t(O + x*h) + 1) & 15) / 4; y = 0;
                    for (; y < y1; y++) O[y + x*h] += (Gy[y] < 0)*PI;
                    for (; y < h - 4; y += 4) STRu(O[y + x*h],
                        ADD(LDu(O[y + x*h]), AND(CMPLT(LDu(Gy[y]), SET(0.f)), SET(PI))));
                        for (; y < h; y++) O[y + x*h] += (Gy[y] < 0)*PI;
                }
            }
            alFree(Gx); alFree(Gy); alFree(M2);
        }

        // normalize gradient magnitude at each location (uses sse)
        static void gradMagNorm(float *M, float *S, int h, int w, float norm) {
            __m128 *_M, *_S, _norm; int i = 0, n = h*w, n4 = n / 4;
            _S = (__m128*) S; _M = (__m128*) M; _norm = SET(norm);
            bool sse = !(size_t(M) & 15) && !(size_t(S) & 15);
            if (sse) { for (; i < n4; i++) *_M++ = MUL(*_M, RCP(ADD(*_S++, _norm))); i *= 4; }
            for (; i < n; i++) M[i] /= (S[i] + norm);
        }

        // helper for gradHist, quantize O and M into O0, O1 and M0, M1 (uses sse)
        static void gradQuantize(float *O, float *M, int *O0, int *O1, float *M0, float *M1,
            int nb, int n, float norm, int nOrients, bool full, bool interpolate)
        {
            // assumes all *OUTPUT* matrices are 4-byte aligned
            int i, o0, o1; float o, od, m;
            __m128i _o0, _o1, *_O0, *_O1; __m128 _o, _od, _m, *_M0, *_M1;
            // define useful constants
            const float oMult = (float)nOrients / (full ? 2 * PI : PI); const int oMax = nOrients*nb;
            const __m128 _norm = SET(norm), _oMult = SET(oMult), _nbf = SET((float)nb);
            const __m128i _oMax = SET(oMax), _nb = SET(nb);
            // perform the majority of the work with sse
            _O0 = (__m128i*) O0; _O1 = (__m128i*) O1; _M0 = (__m128*) M0; _M1 = (__m128*) M1;
            if (interpolate) for (i = 0; i <= n - 4; i += 4) {
                _o = MUL(LDu(O[i]), _oMult); _o0 = CVT(_o); _od = SUB(_o, CVT(_o0));
                _o0 = CVT(MUL(CVT(_o0), _nbf)); _o0 = AND(CMPGT(_oMax, _o0), _o0); *_O0++ = _o0;
                _o1 = ADD(_o0, _nb); _o1 = AND(CMPGT(_oMax, _o1), _o1); *_O1++ = _o1;
                _m = MUL(LDu(M[i]), _norm); *_M1 = MUL(_od, _m); *_M0++ = SUB(_m, *_M1); _M1++;
            }
            else for (i = 0; i <= n - 4; i += 4) {
                _o = MUL(LDu(O[i]), _oMult); _o0 = CVT(ADD(_o, SET(.5f)));
                _o0 = CVT(MUL(CVT(_o0), _nbf)); _o0 = AND(CMPGT(_oMax, _o0), _o0); *_O0++ = _o0;
                *_M0++ = MUL(LDu(M[i]), _norm); *_M1++ = SET(0.f); *_O1++ = SET(0);
            }
            // compute trailing locations without sse
            if (interpolate) for (i; i < n; i++) {
                o = O[i] * oMult; o0 = (int)o; od = o - o0;
                o0 *= nb; if (o0 >= oMax) o0 = 0; O0[i] = o0;
                o1 = o0 + nb; if (o1 == oMax) o1 = 0; O1[i] = o1;
                m = M[i] * norm; M1[i] = od*m; M0[i] = m - M1[i];
            }
            else for (i; i < n; i++) {
                o = O[i] * oMult; o0 = (int)(o + .5f);
                o0 *= nb; if (o0 >= oMax) o0 = 0; O0[i] = o0;
                M0[i] = M[i] * norm; M1[i] = 0; O1[i] = 0;
            }
        }

        // compute nOrients gradient histograms per bin x bin block of pixels
        static void gradHist(float *M, float *O, float *H, int h, int w,
            int bin, int nOrients, int softBin, bool full)
        {
            const int hb = h / bin, wb = w / bin, h0 = hb*bin, w0 = wb*bin, nb = wb*hb;
            const float s = (float)bin, sInv = 1 / s, sInv2 = 1 / s / s;
            float *H0, *H1, *M0, *M1; int x, y; int *O0, *O1;
            O0 = (int*)alMalloc(h*sizeof(int), 16); M0 = (float*)alMalloc(h*sizeof(float), 16);
            O1 = (int*)alMalloc(h*sizeof(int), 16); M1 = (float*)alMalloc(h*sizeof(float), 16);
            // main loop
            for (x = 0; x < w0; x++) {
                // compute target orientation bins for entire column - very fast
                gradQuantize(O + x*h, M + x*h, O0, O1, M0, M1, nb, h0, sInv2, nOrients, full, softBin >= 0);

                if (softBin < 0 && softBin % 2 == 0) {
                    // no interpolation w.r.t. either orienation or spatial bin
                    H1 = H + (x / bin)*hb;
#define GH H1[O0[y]]+=M0[y]; y++;
                    if (bin == 1)      for (y = 0; y < h0;) { GH; H1++; }
                    else if (bin == 2) for (y = 0; y < h0;) { GH; GH; H1++; }
                    else if (bin == 3) for (y = 0; y < h0;) { GH; GH; GH; H1++; }
                    else if (bin == 4) for (y = 0; y < h0;) { GH; GH; GH; GH; H1++; }
                    else for (y = 0; y < h0;) { for (int y1 = 0; y1 < bin; y1++) { GH; } H1++; }
#undef GH

                }
                else if (softBin % 2 == 0 || bin == 1) {
                    // interpolate w.r.t. orientation only, not spatial bin
                    H1 = H + (x / bin)*hb;
#define GH H1[O0[y]]+=M0[y]; H1[O1[y]]+=M1[y]; y++;
                    if (bin == 1)      for (y = 0; y < h0;) { GH; H1++; }
                    else if (bin == 2) for (y = 0; y < h0;) { GH; GH; H1++; }
                    else if (bin == 3) for (y = 0; y < h0;) { GH; GH; GH; H1++; }
                    else if (bin == 4) for (y = 0; y < h0;) { GH; GH; GH; GH; H1++; }
                    else for (y = 0; y < h0;) { for (int y1 = 0; y1 < bin; y1++) { GH; } H1++; }
#undef GH

                }
                else {
                    // interpolate using trilinear interpolation
                    float ms[4], xyd, xb, yb, xd, yd, init; __m128 _m, _m0, _m1;
                    bool hasLf, hasRt; int xb0, yb0;
                    if (x == 0) { init = (0 + .5f)*sInv - 0.5f; xb = init; }
                    hasLf = xb >= 0; xb0 = hasLf ? (int)xb : -1; hasRt = xb0 < wb - 1;
                    xd = xb - xb0; xb += sInv; yb = init; y = 0;
                    // macros for code conciseness
#define GHinit yd=yb-yb0; yb+=sInv; H0=H+xb0*hb+yb0; xyd=xd*yd; \
        ms[0]=1-xd-yd+xyd; ms[1]=yd-xyd; ms[2]=xd-xyd; ms[3]=xyd;
#define GH(H,ma,mb) H1=H; STRu(*H1,ADD(LDu(*H1),MUL(ma,mb)));
                    // leading rows, no top bin
                    for (; y < bin / 2; y++) {
                        yb0 = -1; GHinit;
                        if (hasLf) { H0[O0[y] + 1] += ms[1] * M0[y]; H0[O1[y] + 1] += ms[1] * M1[y]; }
                        if (hasRt) { H0[O0[y] + hb + 1] += ms[3] * M0[y]; H0[O1[y] + hb + 1] += ms[3] * M1[y]; }
                    }
                    // main rows, has top and bottom bins, use SSE for minor speedup
                    if (softBin < 0) for (;; y++) {
                        yb0 = (int)yb; if (yb0 >= hb - 1) break; GHinit; _m0 = SET(M0[y]);
                        if (hasLf) { _m = SET(0, 0, ms[1], ms[0]); GH(H0 + O0[y], _m, _m0); }
                        if (hasRt) { _m = SET(0, 0, ms[3], ms[2]); GH(H0 + O0[y] + hb, _m, _m0); }
                    }
                    else for (;; y++) {
                        yb0 = (int)yb; if (yb0 >= hb - 1) break; GHinit;
                        _m0 = SET(M0[y]); _m1 = SET(M1[y]);
                        if (hasLf) {
                            _m = SET(0, 0, ms[1], ms[0]);
                            GH(H0 + O0[y], _m, _m0); GH(H0 + O1[y], _m, _m1);
                        }
                        if (hasRt) {
                            _m = SET(0, 0, ms[3], ms[2]);
                            GH(H0 + O0[y] + hb, _m, _m0); GH(H0 + O1[y] + hb, _m, _m1);
                        }
                    }
                    // final rows, no bottom bin
                    for (; y < h0; y++) {
                        yb0 = (int)yb; GHinit;
                        if (hasLf) { H0[O0[y]] += ms[0] * M0[y]; H0[O1[y]] += ms[0] * M1[y]; }
                        if (hasRt) { H0[O0[y] + hb] += ms[2] * M0[y]; H0[O1[y] + hb] += ms[2] * M1[y]; }
                    }
#undef GHinit
#undef GH
                }
            }
            alFree(O0); alFree(O1); alFree(M0); alFree(M1);
            // normalize boundary bins which only get 7/8 of weight of interior bins
            if (softBin % 2 != 0) for (int o = 0; o < nOrients; o++) {
                x = 0; for (y = 0; y < hb; y++) H[o*nb + x*hb + y] *= 8.f / 7.f;
                y = 0; for (x = 0; x < wb; x++) H[o*nb + x*hb + y] *= 8.f / 7.f;
                x = wb - 1; for (y = 0; y < hb; y++) H[o*nb + x*hb + y] *= 8.f / 7.f;
                y = hb - 1; for (x = 0; x < wb; x++) H[o*nb + x*hb + y] *= 8.f / 7.f;
            }
        }

        // HOG helper: compute 2x2 block normalization values (padded by 1 pixel)
        static float* hogNormMatrix(float *H, int nOrients, int hb, int wb, int bin) {
            float *N, *N1, *n; int o, x, y, dx, dy, hb1 = hb + 1, wb1 = wb + 1;
            float eps = 1e-4f / 4 / bin / bin / bin / bin; // precise backward equality
            N = (float*)wrCalloc(hb1*wb1, sizeof(float)); N1 = N + hb1 + 1;
            for (o = 0; o<nOrients; o++) for (x = 0; x<wb; x++) for (y = 0; y<hb; y++)
                N1[x*hb1 + y] += H[o*wb*hb + x*hb + y] * H[o*wb*hb + x*hb + y];
            for (x = 0; x<wb - 1; x++) for (y = 0; y<hb - 1; y++) {
                n = N1 + x*hb1 + y; *n = 1 / float(sqrt(n[0] + n[1] + n[hb1] + n[hb1 + 1] + eps));
            }
            x = 0;     dx = 1; dy = 1; y = 0;                  N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            x = 0;     dx = 1; dy = 0; for (y = 0; y<hb1; y++)  N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            x = 0;     dx = 1; dy = -1; y = hb1 - 1;              N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            x = wb1 - 1; dx = -1; dy = 1; y = 0;                  N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            x = wb1 - 1; dx = -1; dy = 0; for (y = 0; y<hb1; y++) N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            x = wb1 - 1; dx = -1; dy = -1; y = hb1 - 1;              N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            y = 0;     dx = 0; dy = 1; for (x = 0; x<wb1; x++)  N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            y = hb1 - 1; dx = 0; dy = -1; for (x = 0; x<wb1; x++)  N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
            return N;
        }

        // HOG helper: compute HOG or FHOG channels
        static void hogChannels(float *H, const float *R, const float *N,
            int hb, int wb, int nOrients, float clip, int type)
        {
#define GETT(blk) t=R1[y]*N1[y-(blk)]; if(t>clip) t=clip; c++;
            const float r = .2357f; int o, x, y, c; float t;
            const int nb = wb*hb, nbo = nOrients*nb, hb1 = hb + 1;
            for (o = 0; o < nOrients; o++) for (x = 0; x < wb; x++) {
                const float *R1 = R + o*nb + x*hb, *N1 = N + x*hb1 + hb1 + 1;
                float *H1 = (type <= 1) ? (H + o*nb + x*hb) : (H + x*hb);
                if (type == 0) for (y = 0; y < hb; y++) {
                    // store each orientation and normalization (nOrients*4 channels)
                    c = -1; GETT(0); H1[c*nbo + y] = t; GETT(1); H1[c*nbo + y] = t;
                    GETT(hb1); H1[c*nbo + y] = t; GETT(hb1 + 1); H1[c*nbo + y] = t;
                }
                else if (type == 1) for (y = 0; y < hb; y++) {
                    // sum across all normalizations (nOrients channels)
                    c = -1; GETT(0); H1[y] += t*.5f; GETT(1); H1[y] += t*.5f;
                    GETT(hb1); H1[y] += t*.5f; GETT(hb1 + 1); H1[y] += t*.5f;
                }
                else if (type == 2) for (y = 0; y < hb; y++) {
                    // sum across all orientations (4 channels)
                    c = -1; GETT(0); H1[c*nb + y] += t*r; GETT(1); H1[c*nb + y] += t*r;
                    GETT(hb1); H1[c*nb + y] += t*r; GETT(hb1 + 1); H1[c*nb + y] += t*r;
                }
            }
#undef GETT
        }

        // compute HOG features
        static void hog(float *M, float *O, float *H, int h, int w, int binSize,
            int nOrients, int softBin, bool full, float clip)
        {
            float *N, *R; 
            const int hb = h / binSize;
            const int wb = w / binSize;
            //const int nb = hb*wb;
            // compute unnormalized gradient histograms
            R = (float*)wrCalloc(wb*hb*nOrients, sizeof(float));
            gradHist(M, O, R, h, w, binSize, nOrients, softBin, full);
            // compute block normalization values
            N = hogNormMatrix(R, nOrients, hb, wb, binSize);
            // perform four normalizations per spatial block
            hogChannels(H, R, N, hb, wb, nOrients, clip, 0);
            wrFree(N); wrFree(R);
        }

        // compute FHOG features
        static void fhog(float *M, float *O, float *H, int h, int w, int binSize,
            int nOrients, int softBin, float clip)
        {
            const int hb = h / binSize, wb = w / binSize, nb = hb*wb, nbo = nb*nOrients;
            float *N, *R1, *R2; int o, x;
            // compute unnormalized constrast sensitive histograms
            R1 = (float*)wrCalloc(wb*hb*nOrients * 2, sizeof(float));
            gradHist(M, O, R1, h, w, binSize, nOrients * 2, softBin, true);
            // compute unnormalized contrast insensitive histograms
            R2 = (float*)wrCalloc(wb*hb*nOrients, sizeof(float));
            for (o = 0; o < nOrients; o++) for (x = 0; x < nb; x++)
                R2[o*nb + x] = R1[o*nb + x] + R1[(o + nOrients)*nb + x];
            // compute block normalization values
            N = hogNormMatrix(R2, nOrients, hb, wb, binSize);
            // normalized histograms and texture channels
            hogChannels(H + nbo * 0, R1, N, hb, wb, nOrients * 2, clip, 1);
            hogChannels(H + nbo * 2, R2, N, hb, wb, nOrients * 1, clip, 1);
            hogChannels(H + nbo * 3, R1, N, hb, wb, nOrients * 2, clip, 2);
            wrFree(N); wrFree(R1); wrFree(R2);
        }
#pragma endregion

#pragma region --- edge nms ---
        // return I[x,y] via bilinear interpolation
        static float interp(float *I, int h, int w, float x, float y) {
            x = x<0 ? 0 : (x>w - 1.001f ? w - 1.001f : x);
            y = y<0 ? 0 : (y>h - 1.001f ? h - 1.001f : y);
            int x0 = int(x), y0 = int(y), x1 = x0 + 1, y1 = y0 + 1;
            float dx0 = x - x0, dy0 = y - y0, dx1 = 1 - dx0, dy1 = 1 - dy0;
            return I[x0*h + y0] * dx1*dy1 + I[x1*h + y0] * dx0*dy1 +
                I[x0*h + y1] * dx1*dy0 + I[x1*h + y1] * dx0*dy0;
        }
#pragma endregion
    };
}