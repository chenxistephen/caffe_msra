//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      brief descriptor implementation
//
//------------------------------------------------------------------------------
#pragma once

#include "overlapblock.h"

//+-----------------------------------------------------------------------------
// Class: CBriefTable
//------------------------------------------------------------------------------
inline CPoint randomSampleUnitGaussian2D(int iPatchSize, double sigma)
{
    const int P = 10000;

    double u, v, w;
    double rP = 2.0/double(P);
    do {
        u = (rand() % P) * rP - 1.0;
        v = (rand() % P) * rP - 1.0;
        w = u*u + v*v;
    } while (w >= 1.0);
    
    CVec2d sd = CVec2d(u,v) * sigma*sqrt((-2.0 * log(w))/w);

    CPoint s = CPoint(int(sd.x+0.5)+iPatchSize/2, int(sd.y+0.5)+iPatchSize/2);
    s.x = VtMax(0, VtMin((int)s.x, iPatchSize-1));
    s.y = VtMax(0, VtMin((int)s.y, iPatchSize-1));

    return s;
}

template<int DESC_SIZE>
class BriefTable
{
public:
    static const int iDescSize = DESC_SIZE;

    BriefTable(int iPatchSize) { Initialize(iPatchSize); }
    BriefTable() {};
    void Initialize(int iPatchSize);

    unsigned short OffsetA(int i) const
    { return m_table[2*i]; }

    unsigned short OffsetB(int i) const
    { return m_table[2*i+1]; }

private:
    unsigned short m_table[2*DESC_SIZE];
};

template <int DESC_SIZE>
void BriefTable<DESC_SIZE>::Initialize(int iPatchSize)
{
    srand(42);
    double sigma = 0.2 * double(iPatchSize);

    for( int i = 0; i < DESC_SIZE; ) 
    {
        vt::CPoint a = randomSampleUnitGaussian2D(iPatchSize, sigma);
        vt::CPoint b = randomSampleUnitGaussian2D(iPatchSize, sigma);

        if ( abs(a.x-b.x) > 1 || abs(a.y-b.y) > 1 )
        {
            int x0 = a.y * iPatchSize + a.x;
            int x1 = b.y * iPatchSize + b.x;
            int j  = 0;
            for( ; j < i; j++ )
            {
                if( OffsetA(j) == x0 || OffsetB(j) == x1 )
                {
                    break;
                }
            }
            if( j == i )
            {
                m_table[2*i]   = (unsigned short)x0;
                m_table[2*i+1] = (unsigned short)x1;
                i++;
            }
        }
    }

#if 0
    // code to generate table for 'faster Brief descriptor generation' code below
    const int stride = 640;
    for( int i = 0; i < DESC_SIZE; i++) 
    {
        int tableOffsetA = m_table[2*i+0];
        int offsetAx = tableOffsetA%iPatchSize;
        int offsetAy = tableOffsetA/iPatchSize;
        int offsetA = ((offsetAy-(iPatchSize/2))*stride)+offsetAx;

        int tableOffsetB = m_table[2*i+1];
        int offsetBx = tableOffsetB%iPatchSize;
        int offsetBy = tableOffsetB/iPatchSize;
        int offsetB = ((offsetBy-(iPatchSize/2))*stride)+offsetBx;

        if ((i%4)==0) printf("\n");
        printf("{% 5d,% 5d},",offsetA,offsetB); 
    }
    printf("\n");
#endif
}

//+-----------------------------------------------------------------------------
// Class: BriefDesc
//------------------------------------------------------------------------------
template<int DESC_SIZE>
class BriefDesc
{
public:
    static const int iDescSize = DESC_SIZE;

    void Clear()
    { ZeroMemory(m_desc, sizeof(m_desc)); }

    void SetBit( int i )
    { m_desc[i/32] |= 1<<(i&31); }

    const UInt32* GetDescriptor() const
    { return m_desc; }

//protected:
    UInt32 m_desc[(DESC_SIZE+31)/32];
};

//+-----------------------------------------------------------------------------
// Class: BriefAvgDesc
//------------------------------------------------------------------------------
template<int DESC_SIZE>
class BriefAvgDesc: public BriefDesc<DESC_SIZE>
{
public:
    void Clear()
    { 
        m_cnt = 0;
        ZeroMemory(m_avg, sizeof(m_avg));
        BriefDesc<DESC_SIZE>::Clear();  
    }
#if 0
    void Acc(const BriefDesc<DESC_SIZE>& src)
    {
        const UInt32* pSrc = src.GetDescriptor();
        UInt32* pDesc      = m_desc;
        for( int i = 0; i < DESC_SIZE; i+=32, pSrc++, pDesc++ )
        {
            *pDesc = 0;
            for( int j = 0; j < 32; j++ )
            {
                int bit     = ((*pSrc) >> j)&1;
                m_avg[i+j] += (bit^1) - bit;
                *pDesc     |= ((m_avg[i+j]>>15)&1) << j;
            }
        }
    }
#else
    void Acc(const BriefDesc<DESC_SIZE>& src)
    {
        int shift, rnd;
        if( m_cnt < 7 )
        {
            m_cnt++;
            if( m_cnt < 2 )
            {
                shift = rnd = 0;
            }
            else
            {
                shift = (m_cnt<4)? 1: 2; 
                rnd   = 1 << (shift-1);
            }
        }
        else
        {
            shift = 3;
            rnd   = 4;
        }

        const UInt32* pSrc = src.GetDescriptor();
        UInt32* pDesc      = m_desc;
        for( int i = 0; i < DESC_SIZE; i+=32, pSrc++, pDesc++ )
        {
            *pDesc = 0;
            for( int j = 0; j < 32; j++ )
            {
                int bit    = ((*pSrc) >> j)&1;
                int accw   = (1<<shift)-1;
                int newval = ((bit^1) - bit)*16;
                m_avg[i+j] = (short)((m_avg[i+j]*accw + newval + rnd) >> shift);
                *pDesc    |= ((m_avg[i+j]>>15)&1) << j;
            }
        }
    }
#endif

protected:
    int    m_cnt;
    short  m_avg[DESC_SIZE];
};

//+-----------------------------------------------------------------------------
// function: BriefEvaluatePatch
//------------------------------------------------------------------------------
template<int DESC_SIZE, int PATCH_SIZE>
void inline BriefEvaluatePatch(BriefDesc<DESC_SIZE>& desc,
                        const BriefTable<DESC_SIZE>& BT, 
                        const Byte* pPatch,
                        size_t stride)
{
    desc.Clear();    
    for (int i = 0; i < DESC_SIZE; i++) 
    {
        int offsetA = BT.OffsetA(i);
        int offsetAx = offsetA%PATCH_SIZE;
        int offsetAy = offsetA/PATCH_SIZE;
        int offsetB = BT.OffsetB(i);
        int offsetBx = offsetB%PATCH_SIZE;
        int offsetBy = offsetB/PATCH_SIZE;
        if( pPatch[(offsetAy*stride)+offsetAx] > pPatch[(offsetBy*stride)+offsetBx] )
            desc.SetBit(i);
    }
}

template<int DESC_SIZE, int PATCH_SIZE>
void inline BriefEvaluatePatch(BriefDesc<DESC_SIZE>& desc, 
                        const BriefTable<DESC_SIZE>& BT,
                        const HARRIS_FEATURE_POINT& fp, const CByteImg& src)
{
    int patchx = F2I(fp.x)-(PATCH_SIZE/2);
    int patchy = F2I(fp.y)-(PATCH_SIZE/2);
    VT_ASSERT( (patchx>=0) && (patchy>=0) );
    VT_ASSERT( ((patchx+PATCH_SIZE)<=src.Width()) && ((patchy+PATCH_SIZE)<=src.Height()) );

    BriefEvaluatePatch<DESC_SIZE,PATCH_SIZE>(desc, BT, src.Ptr(patchx,patchy), src.StrideBytes());
}

//+-----------------------------------------------------------------------------
//  function BriefFindMatch
// 
//  find the closest brief feature from a list
// 
//------------------------------------------------------------------------------
//
// TODO: perhaps remove iFrameWidth/iFrameHeight params as these are only used to
//       scale maxDist to pixel coordinates
//
template<int iDescSize>
void BriefFindMatch(int& iNearest, int& iNearestDistance,
					const vt::vector<HARRIS_FEATURE_POINT>& vecFPSearch,
                    const vt::vector<int>* pvecFpSearchList,
                    const vt::vector<BriefDesc<iDescSize>>& vecDescSearch, 
                    const CVec2f& fp, const BriefDesc<iDescSize>& desc,
                    int iFrameWidth, int iFrameHeight, float maxDist, float matchTestRatio)
{
    float fFrameDim = (float)VtMax(iFrameWidth, iFrameHeight);
    int   maxDist2  = F2I(fFrameDim*maxDist);
    maxDist2 = maxDist2*maxDist2;

    int iMatch    = -1;
    int iMatch2   = -1;
    int iMinDist  = iDescSize+1;
    int iMinDist2 = iDescSize+1;

	int iCount = int(pvecFpSearchList? pvecFpSearchList->size(): vecFPSearch.size());

    for( int i = 0; i < iCount; i++ )
    {
        int index = pvecFpSearchList? (*pvecFpSearchList)[i]: i;
        const HARRIS_FEATURE_POINT& fpSearch = vecFPSearch[index];

        // only consider nearby features
        if( CVec2f(fpSearch.x-fp.x, fpSearch.y-fp.y).MagnitudeSq() > maxDist2)
        {
            continue;
        }

        int iDist = BriefDistance(vecDescSearch[index], desc);
        if( iDist < iMinDist2 )
        {
            if( iDist < iMinDist )
            {
                iMatch2   = iMatch;
                iMinDist2 = iMinDist;
                iMatch    = index;
                iMinDist  = iDist;
            }
            else
            {
                iMatch2   = index;
                iMinDist2 = iDist;
            }
        }
    }

    // do the ratio test
    if( iMatch2 != -1 && float(iMinDist)/float(iMinDist2) > matchTestRatio )
    {
        iMatch = -1;
    }

    iNearest = iMatch;
    iNearestDistance = iMinDist;
}


template<int iDescSize>
void BriefFindMatch(int& iNearest, int& iNearestDistance,
					const vt::vector<HARRIS_FEATURE_POINT>& vecFPSearch,
                    const COverlappedBlocks<int>& ovlblk,
                    const vt::vector<BriefDesc<iDescSize>>& vecDescSearch, 
                    const CVec2f& fp, const BriefDesc<iDescSize>& desc,
                    int iFrameWidth, int iFrameHeight, float maxDist, float matchTestRatio)
{
	BriefFindMatch(iNearest, iNearestDistance, 
				   vecFPSearch, ovlblk.GetElements((int)fp.x,(int)fp.y), 
				   vecDescSearch, fp, desc, iFrameWidth, iFrameHeight, maxDist, matchTestRatio); 
}

#if 0
void inline DumpPatch(const wchar_t* pname, const HARRIS_FEATURE_POINT& fp, const CByteImg& src,
               int iPatchSize)
{
    vt::CPoint cntr(F2I(fp.x), F2I(fp.y));
    vt::CRect  rctPatch(cntr.x-iPatchSize/2, cntr.y-iPatchSize/2,
                        cntr.x+iPatchSize/2, cntr.y+iPatchSize/2);

    CByteImg imgPatch;
    src.Share(imgPatch, &rctPatch);

    CByteImg img4x;
    VtResizeImage(img4x, vt::CRect(0,0,rctPatch.Width()*4,rctPatch.Height()*4),
                  imgPatch);
    VtSaveImage(pname, img4x);
}
#endif

//+-----------------------------------------------------------------------------
// function: BriefDistance
//------------------------------------------------------------------------------
#define BC_TWO(c)     (0x1u << (c))
#define BC_MASK(c)    (((unsigned int)(-1)) / (BC_TWO(BC_TWO(c)) + 1u))
#define BC_COUNT(x,c) ((x) & BC_MASK(c)) + (((x) >> (BC_TWO(c))) & BC_MASK(c))

inline int bitcount(UInt32 n)
{
    // TODO: implement an SSE4 version of this via the POPCNT instruction
    //       http://msdn.microsoft.com/en-us/library/bb514083.aspx
    n = BC_COUNT(n, 0) ;
    n = BC_COUNT(n, 1) ;
    n = BC_COUNT(n, 2) ;
    n = BC_COUNT(n, 3) ;
    n = BC_COUNT(n, 4) ;
    return n ;
}

template<int DESC_SIZE>
int BriefDistance(const BriefDesc<DESC_SIZE>& a, const BriefDesc<DESC_SIZE>& b)
{
    int dist = 0;
    for( int i = 0; i < (DESC_SIZE+31)/32; i++ )
    {
        dist += bitcount(a.GetDescriptor()[i] ^ b.GetDescriptor()[i]);
        // TODO: add an early out max distance test
    }
    return dist;
}

//------------------------------------------------------------------------------
//
// faster Brief descriptor generation - functions with compiled-in constants to
// locations of Brief sample points for specific descriptor bit count, patch
// size, and source image stride
//
//------------------------------------------------------------------------------

// table of Brief descriptor offsets relative to patch x,y coord (0,patchSize/2)
//
// offsets are relative to middle scanline of the patch since processors
// (ARM at least) support a limited-range signed absolute offset for a
// single instruction byte load; for ARM, the offsets in this table must
// be in the range +/- 4096 to do the fastest byte load
//
// note that this is about as large a stride as can be optimally supported
// for a patch size of 24 and the +/- 4096 offset range since 12*320 = 3840,
// which is only just under 4096; the code will function for a larger offset
// but will generate slower code for positions that are outside this range
//
// descriptor: 128bits; patch size: 24x24; image stride: 640
const int32_t BriefTableOffsets_d128_p24_s640[128][2] = {
{-3824, 1935},{-3187, 3848},{-1912,    9},{ 1292, 1295},
{ 3860,-3194},{ 2578, 4490},{-2548,-1278},{-1907,-3824},
{ 4488,-3828},{ -628,-1265},{-2541,-6390},{-1908, 1285},
{ 1286,  650},{-1267,   16},{   13,-3833},{-2552, 3214},
{ -623,-1279},{  652,-4468},{ 7059,-1912},{ 2571, 2565},
{-2549, 5778},{    6,   12},{-1905,-1271},{    9, 1294},
{ -621,   13},{ 1930,   17},{-5748, 2573},{  656, 4489},
{-4464,   10},{-1270,-1908},{ 1289,-2551},{-2551,-1913},
{   12, 3215},{-3828,-2556},{-1911, 1933},{-3826, 2572},
{-1271,-1267},{ -629, 5768},{-1906,  653},{ 1938,-2552},
{-1268,-1904},{-3825, 4497},{ 5137, 3219},{  651, -624},
{-1266,  649},{-3188,   14},{ 1299,-1909},{ -636,-1260},
{-3830,    3},{-1904,   11},{ 4485,-1268},{ 3843,    6},
{-3180,-3825},{ 1926,    8},{-1910, -623},{ -625, 1928},
{-1263,-1273},{    5, 1299},{-1901, 1937},{ 1294, -628},
{-4465, 4482},{ 3852,-4464},{ 1293,-3835},{    7,-1269},
{  661, 3861},{   21,  645},{ 7052, 5130},{-1273,-2544},
{   10,  656},{-3191, 1293},{ 7045, 1929},{  654,-6385},
{ 3212, 5774},{-1914, 1932},{ 5126, 4493},{ 1932, 1924},
{ -626,-3832},{-3186,-2541},{ 5131,-2550},{ 3220,-6382},
{ 3211, -630},{ 1940, 2577},{-5112,-2547},{ 3855,-2548},
{-1909, -622},{ -627, 5132},{ 3216,  654},{-1897, 2570},
{  653,-1902},{ 1936, 3851},{ -634, 4488},{ -631, 3212},
{ 1298, 3852},{ 4497, 2564},{  650,   19},{   15,-1272},
{ 2561,-1914},{ 1928,-1907},{ 2572,   22},{ 1934, 3853},
{-6384,-2546},{ 2568, 1931},{-5111,-1270},{ 3849, 1291},
{-1262,-2545},{-3834, 3863},{ 1296,-1275},{ 1291, -631},
{ 1290, 5137},{    8,-1261},{ 3856,  652},{ 3210,-3819},
{-3829,  658},{   19, -626},{-1272,-1264},{ 2569,-1906},
{-4463, 1301},{-1913,  659},{ 4492, 4485},{ -630,  648},
{ 3218,-1911},{ 1287,  655},{ 3214, 1292},{ 3213, 5767},
{ -624,-2542},{ 1295,-1266},{-1257,-4463},{ -633, -625},
};
// descriptor: 128bits; patch size: 24x24; image stride: 320
const int32_t BriefTableOffsets_d128_p24_s320[128][2] = {
{-1904,  975},{-1587, 1928},{ -952,    9},{  652,  655},
{ 1940,-1594},{ 1298, 2250},{-1268, -638},{ -947,-1904},
{ 2248,-1908},{ -308, -625},{-1261,-3190},{ -948,  645},
{  646,  330},{ -627,   16},{   13,-1913},{-1272, 1614},
{ -303, -639},{  332,-2228},{ 3539, -952},{ 1291, 1285},
{-1269, 2898},{    6,   12},{ -945, -631},{    9,  654},
{ -301,   13},{  970,   17},{-2868, 1293},{  336, 2249},
{-2224,   10},{ -630, -948},{  649,-1271},{-1271, -953},
{   12, 1615},{-1908,-1276},{ -951,  973},{-1906, 1292},
{ -631, -627},{ -309, 2888},{ -946,  333},{  978,-1272},
{ -628, -944},{-1905, 2257},{ 2577, 1619},{  331, -304},
{ -626,  329},{-1588,   14},{  659, -949},{ -316, -620},
{-1910,    3},{ -944,   11},{ 2245, -628},{ 1923,    6},
{-1580,-1905},{  966,    8},{ -950, -303},{ -305,  968},
{ -623, -633},{    5,  659},{ -941,  977},{  654, -308},
{-2225, 2242},{ 1932,-2224},{  653,-1915},{    7, -629},
{  341, 1941},{   21,  325},{ 3532, 2570},{ -633,-1264},
{   10,  336},{-1591,  653},{ 3525,  969},{  334,-3185},
{ 1612, 2894},{ -954,  972},{ 2566, 2253},{  972,  964},
{ -306,-1912},{-1586,-1261},{ 2571,-1270},{ 1620,-3182},
{ 1611, -310},{  980, 1297},{-2552,-1267},{ 1935,-1268},
{ -949, -302},{ -307, 2572},{ 1616,  334},{ -937, 1290},
{  333, -942},{  976, 1931},{ -314, 2248},{ -311, 1612},
{  658, 1932},{ 2257, 1284},{  330,   19},{   15, -632},
{ 1281, -954},{  968, -947},{ 1292,   22},{  974, 1933},
{-3184,-1266},{ 1288,  971},{-2551, -630},{ 1929,  651},
{ -622,-1265},{-1914, 1943},{  656, -635},{  651, -311},
{  650, 2577},{    8, -621},{ 1936,  332},{ 1610,-1899},
{-1909,  338},{   19, -306},{ -632, -624},{ 1289, -946},
{-2223,  661},{ -953,  339},{ 2252, 2245},{ -310,  328},
{ 1618, -951},{  647,  335},{ 1614,  652},{ 1613, 2887},
{ -304,-1262},{  655, -626},{ -617,-2223},{ -313, -305},
};
// descriptor: 128bits; patch size: 24x24; image stride: 180
const int32_t BriefTableOffsets_d128_p24_s180[128][2] = {
{-1064,  555},{ -887, 1088},{ -532,    9},{  372,  375},
{ 1100, -894},{  738, 1270},{ -708, -358},{ -527,-1064},
{ 1268,-1068},{ -168, -345},{ -701,-1790},{ -528,  365},
{  366,  190},{ -347,   16},{   13,-1073},{ -712,  914},
{ -163, -359},{  192,-1248},{ 1999, -532},{  731,  725},
{ -709, 1638},{    6,   12},{ -525, -351},{    9,  374},
{ -161,   13},{  550,   17},{-1608,  733},{  196, 1269},
{-1244,   10},{ -350, -528},{  369, -711},{ -711, -533},
{   12,  915},{-1068, -716},{ -531,  553},{-1066,  732},
{ -351, -347},{ -169, 1628},{ -526,  193},{  558, -712},
{ -348, -524},{-1065, 1277},{ 1457,  919},{  191, -164},
{ -346,  189},{ -888,   14},{  379, -529},{ -176, -340},
{-1070,    3},{ -524,   11},{ 1265, -348},{ 1083,    6},
{ -880,-1065},{  546,    8},{ -530, -163},{ -165,  548},
{ -343, -353},{    5,  379},{ -521,  557},{  374, -168},
{-1245, 1262},{ 1092,-1244},{  373,-1075},{    7, -349},
{  201, 1101},{   21,  185},{ 1992, 1450},{ -353, -704},
{   10,  196},{ -891,  373},{ 1985,  549},{  194,-1785},
{  912, 1634},{ -534,  552},{ 1446, 1273},{  552,  544},
{ -166,-1072},{ -886, -701},{ 1451, -710},{  920,-1782},
{  911, -170},{  560,  737},{-1432, -707},{ 1095, -708},
{ -529, -162},{ -167, 1452},{  916,  194},{ -517,  730},
{  193, -522},{  556, 1091},{ -174, 1268},{ -171,  912},
{  378, 1092},{ 1277,  724},{  190,   19},{   15, -352},
{  721, -534},{  548, -527},{  732,   22},{  554, 1093},
{-1784, -706},{  728,  551},{-1431, -350},{ 1089,  371},
{ -342, -705},{-1074, 1103},{  376, -355},{  371, -171},
{  370, 1457},{    8, -341},{ 1096,  192},{  910,-1059},
{-1069,  198},{   19, -166},{ -352, -344},{  729, -526},
{-1243,  381},{ -533,  199},{ 1272, 1265},{ -170,  188},
{  918, -531},{  367,  195},{  914,  372},{  913, 1627},
{ -164, -702},{  375, -346},{ -337,-1243},{ -173, -165},
};

// compute a single bit of the descriptor
#if 0
// with branches
#define CBDbit(_t_,_w_,_i_) \
    if (patch[BriefTableOffsets_##_t_[(_w_*32)+_i_][0]] > \
        patch[BriefTableOffsets_##_t_[(_w_*32)+_i_][1]]) { ret |= (1<<_i_); };
#else  
// with no branches - faster on current x86 processors
#define CBDbit(_t_,_w_,_i_) \
    { union { uint32_t ui; int32_t ii; } i; \
    i.ii = patch[BriefTableOffsets_##_t_[(_w_*32)+_i_][0]] - \
           patch[BriefTableOffsets_##_t_[(_w_*32)+_i_][1]]; \
    i.ui >>= 31; i.ui <<=_i_; ret |= i.ui; }
#endif

// compute a 32 bit word of the multi-word descriptor 
#define CBDdword(_t_,_w_) \
{ \
    register uint32_t ret = 0x0; \
    CBDbit(_t_,_w_, 0); CBDbit(_t_,_w_, 1); CBDbit(_t_,_w_, 2); CBDbit(_t_,_w_, 3); \
    CBDbit(_t_,_w_, 4); CBDbit(_t_,_w_, 5); CBDbit(_t_,_w_, 6); CBDbit(_t_,_w_, 7); \
    CBDbit(_t_,_w_, 8); CBDbit(_t_,_w_, 9); CBDbit(_t_,_w_,10); CBDbit(_t_,_w_,11); \
    CBDbit(_t_,_w_,12); CBDbit(_t_,_w_,13); CBDbit(_t_,_w_,14); CBDbit(_t_,_w_,15); \
    CBDbit(_t_,_w_,16); CBDbit(_t_,_w_,17); CBDbit(_t_,_w_,18); CBDbit(_t_,_w_,19); \
    CBDbit(_t_,_w_,20); CBDbit(_t_,_w_,21); CBDbit(_t_,_w_,22); CBDbit(_t_,_w_,23); \
    CBDbit(_t_,_w_,24); CBDbit(_t_,_w_,25); CBDbit(_t_,_w_,26); CBDbit(_t_,_w_,27); \
    CBDbit(_t_,_w_,28); CBDbit(_t_,_w_,29); CBDbit(_t_,_w_,30); CBDbit(_t_,_w_,31); \
    desc.m_desc[_w_] = ret; \
} 

// compute brief descriptor for descriptor size of 128, patch size
// of 24x24, and img strides of 640,320,180
inline void ComputeBriefDescriptor_d128_p24_s640(BriefDesc<128>& desc, HARRIS_FEATURE_POINT& fp, 
    const CByteImg& img)
{
    const BYTE* patch = img.Ptr(F2I(fp.x)-12,F2I(fp.y));
    CBDdword(d128_p24_s640,0);
    CBDdword(d128_p24_s640,1);
    CBDdword(d128_p24_s640,2);
    CBDdword(d128_p24_s640,3);
}
inline void ComputeBriefDescriptor_d128_p24_s320(BriefDesc<128>& desc, HARRIS_FEATURE_POINT& fp, 
    const CByteImg& img)
{
    const BYTE* patch = img.Ptr(F2I(fp.x)-12,F2I(fp.y));
    CBDdword(d128_p24_s320,0);
    CBDdword(d128_p24_s320,1);
    CBDdword(d128_p24_s320,2);
    CBDdword(d128_p24_s320,3);
}
inline void ComputeBriefDescriptor_d128_p24_s180(BriefDesc<128>& desc, HARRIS_FEATURE_POINT& fp, 
    const CByteImg& img)
{
    const BYTE* patch = img.Ptr(F2I(fp.x)-12,F2I(fp.y));
    CBDdword(d128_p24_s180,0);
    CBDdword(d128_p24_s180,1);
    CBDdword(d128_p24_s180,2);
    CBDdword(d128_p24_s180,3);
}

//------------------------------------------------------------------------------

// end
