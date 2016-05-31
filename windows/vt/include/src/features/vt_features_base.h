//+-----------------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description: base interfaces for the components in this library 
//
//------------------------------------------------------------------------------
#pragma once

//------------------------------------------------------------------------------
// how stabilization components report the current range of results that they
// have buffered up 
//------------------------------------------------------------------------------
struct BUFFER_RANGE
{
	int first_frame;
	int frame_count;

	bool InRange(int frame)
	{ return frame >= first_frame && frame < first_frame+frame_count; }

	BUFFER_RANGE() : first_frame(-1), frame_count(0)
	{}
};

//+-----------------------------------------------------------------------------
//
// Class: IDelayResult
//
//------------------------------------------------------------------------------
class IDelayResult
{
public:
    virtual ~IDelayResult() = 0
	{} 

	virtual int GetMaxDelay() = 0;

	virtual BUFFER_RANGE GetResultsRange() = 0;
};

//+-----------------------------------------------------------------------------
//
// Interface: IFeatureWarpCompute
//
// Synposis: base interface for components that read a set of feature
//           matches and produce a warp transform to 'stabilize' these features.
//
//------------------------------------------------------------------------------
class IFeatureWarpCompute: public IDelayResult
{
public:
	virtual ~IFeatureWarpCompute() = 0
    {}

    virtual HRESULT Begin() = 0;

    virtual HRESULT PushFeatures(const vt::PointMatch* pMatches,
                                 int iMatchCount) = 0;

    virtual HRESULT GetResult(vt::CMtx3x3f& xfrm, int frameNumber) = 0;

    virtual HRESULT End() = 0;
};

//+-----------------------------------------------------------------------------
//
// Class: CRollingBuffer
//
// Synposis: container that holds a rolling set of results
//
//------------------------------------------------------------------------------
template <class T>
class CRollingBuffer
{
public:
    CRollingBuffer() : m_base(0), m_count(0)
    {}

    HRESULT resize(size_t size)
    {
        m_base  = 0;
        m_count = 0;
        return m_data.resize(size);
    }

    void clear()
    {
        m_base  = 0;
        m_count = 0;
        m_data.clear();
    }

    size_t size(void) const
    { return m_data.size(); }

    void advance()
    { 
        m_count++;
        m_base++;
        if( m_base >= (int)m_data.size() )
        {
            m_base = 0;
        }
    }

    void advance(const T& v)
    { 
        VT_ASSERT( m_data.size() > 0 );
        m_data[m_base] = v;
        advance();
    }

    int get_first_id() const
    { return m_count-get_available_count(); }

    int get_last_id() const
    { return m_count-1; }

    int get_available_count() const
    { return m_count<(int)size()? m_count: (int)size(); }

    int get_total_count() const
    { return m_count; }

    const T& get(int id) const
    {
        VT_ASSERT( id >= get_first_id() && id <= get_last_id() );
        return m_data[bound_index(m_base-(m_count-id))];
    }

    T& get(int id)
    {
        VT_ASSERT( id >= get_first_id() && id <= get_last_id() );
        return m_data[bound_index(m_base-(m_count-id))];
    }

    const T& operator[](int id) const
    { return get(id); }

    T& operator[](int id)
    { return get(id); }

    const T& last() const
    { return get(get_last_id()); }

    T& last()
    { return get(get_last_id()); }

    void copy_to(T* valarray, int first, int len) const
    {
        for (int i=0; i<len; i++)
        {
            valarray[i] = get(first+i);
        }
    }

    const T& buffer(int idx) const
    { return m_data[bound_index(m_base+idx)]; }

    T& buffer(int idx)
    { return m_data[bound_index(m_base+idx)]; }

private:
    // no copy allowed
    CRollingBuffer(const CRollingBuffer&);
	CRollingBuffer &operator=(const CRollingBuffer&);

protected:
    int bound_index(int idx) const
    {
        idx = (idx < 0)? idx+(int)size(): (idx >= (int)size())? idx-(int)size(): idx;
        VT_ASSERT(idx >= 0 && idx < (int)size());
        return idx;
    }

protected:
    int m_count; 
    int m_base;
    vt::vector<T> m_data;
};

//+-----------------------------------------------------------------------------
//
// Struct: FEAT_SIMILARITY
//
// Synposis: parameters that describe a similarity transform
//
//------------------------------------------------------------------------------
struct FEAT_SIMILARITY
{
    float tx;
    float ty;
    float s;
    float r;

    FEAT_SIMILARITY() : tx(0), ty(0), s(1.f), r(0.f)
    {}
    FEAT_SIMILARITY(float tx_, float ty_, float s_, float r_) : 
        tx(tx_), ty(ty_), s(s_), r(r_)
    {}
};
