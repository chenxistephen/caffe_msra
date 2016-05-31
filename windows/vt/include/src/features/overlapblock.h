//+-----------------------------------------------------------------------
//
//  Copyright (c) Microsoft Corporation.  All rights reserved.
//
//  Description:
//      overlapped block structure
//
//  used to speed up the feature descriptor search by generating a 
//  block-partitioned data structure of indices or pointers to a list
//  of features, where the list consists of all of the features that
//  are either in or within an overlap dimension of that block
//
//  usage intent is that the type T is an index or pointer
//
//------------------------------------------------------------------------
#pragma once

template <class T>
class COverlappedBlocks
{
public:
    COverlappedBlocks()
    {}

	// reset the individual block data.  use this method when re-using 
	// COverlappedBlocks with the same sized images
	void Reset()
	{
		for( int i = 0; i < (int)m_data.size(); i++ )
		{
			for( int j = 0; j < (int)m_data[i].size(); j++ )
			{
				m_data[i][j].resize(0);
			}
		}
	}

    // create a structure with 1<<blockDimPow2 square sized blocks and a
    // 1<<overlapDimPow2 overlap from each block boundary
    HRESULT Create(const vt::CRect& rect, UInt32 blockDimPow2, UInt32 overlapDimPow2)
    {
    	VT_HR_BEGIN()

		m_data.clear();
        m_rect = rect;
        m_blockDP2 = blockDimPow2;
        m_blockOvlDP2 = overlapDimPow2;

        // code below requires that the overlap dimension not be larger
        // than the block dimension (since it only adds to immediately 
        // adjacent blocks), so the block size needs to be clamped to the
        // overlap size
        m_blockDP2 = VtMax(m_blockDP2, m_blockOvlDP2);

        // compute the number of blocks
        UInt32 spaceW = m_rect.Width();
        UInt32 spaceH = m_rect.Height();
        m_blockSize = 1<<m_blockDP2;
        m_blocksW = (spaceW+m_blockSize-1)/m_blockSize;
        m_blocksH = (spaceH+m_blockSize-1)/m_blockSize;

        // allocate block data struct - always padded by one in both dimensions
        VT_HR_EXIT( m_data.resize(m_blocksW+2) );
        for (int i=0; i<(int)m_data.size(); i++) 
        { VT_HR_EXIT( m_data[i].resize(m_blocksH+2) ); }

        // masks for working with addresses
        m_blkM = (1<<m_blockDP2)-1;
        m_ovlM = m_blkM & ~((1<<m_blockOvlDP2)-1);

    	VT_HR_END()
    }

    // add element to blocks that include or are within the overlap of the position
    HRESULT AddElement(float fx, float fy, T data)
    {
    	VT_HR_BEGIN()

        int x = F2I(fx);
        int y = F2I(fy);

        // form address for this block and add it to it's own block
        int bx = (x>>m_blockDP2)+1;
        int by = (y>>m_blockDP2)+1;
		VT_ASSERT( bx >= 1 && bx < m_blocksW+1 && 
			       by >= 1 && by < m_blocksH+1 );
        VT_HR_EXIT( m_data[bx][by].push_back(data) );

        // form boundary booleans
        bool bxm = (x & m_ovlM)?false:true; 
        bool bym = (y & m_ovlM)?false:true;
        bool bxx = ((x+(1<<m_blockOvlDP2)) & m_ovlM)?false:true;
        bool byx = ((y+(1<<m_blockOvlDP2)) & m_ovlM)?false:true;

        // conditionally add to four side adjacent blocks
        if (bxm) { m_data[bx-1][by  ].push_back(data); }
        if (bym) { m_data[bx  ][by-1].push_back(data); }
        if (bxx) { m_data[bx+1][by  ].push_back(data); }
        if (byx) { m_data[bx  ][by+1].push_back(data); }

        // conditionally add to four corner adjacent blocks
        if (bxm&bym) { m_data[bx-1][by-1].push_back(data); }
        if (bxm&byx) { m_data[bx-1][by+1].push_back(data); }
        if (bxx&bym) { m_data[bx+1][by-1].push_back(data); }
        if (bxx&byx) { m_data[bx+1][by+1].push_back(data); }

    	VT_HR_END()
    }

    // return list of T's for the block corresponding to the position
    const vt::vector<T>* GetElements(int x, int y) const
    {
        int bx = (x>>m_blockDP2)+1;
        int by = (y>>m_blockDP2)+1;
        return &(m_data[bx][by]);
    }

    // copy the contents of an existing COverlappedBlocks
    HRESULT CopyFrom(const COverlappedBlocks<T>& src)
    {
        VT_HR_BEGIN()

        m_rect        = src.m_rect;
        m_blockDP2    = src.m_blockDP2;   
        m_blockOvlDP2 = src.m_blockOvlDP2;
        m_blocksW     = src.m_blocksW;    
        m_blocksH     = src.m_blocksH;    
        m_blockSize   = src.m_blockSize;  
        m_blkM        = src.m_blkM;       
        m_ovlM        = src.m_ovlM;       

        VT_HR_EXIT( m_data.resize( src.m_data.size() ) );
        for (int i=0; i<(int)src.m_data.size(); i++) 
        { VT_HR_EXIT( m_data[i].resize( src.m_data[i].size() ) ); }

        for( int i = 0; i < (int)src.m_data.size(); i++ )
		{
			for( int j = 0; j < (int)src.m_data[i].size(); j++ )
			{
                for (int k = 0; k < (int)src.m_data[i][j].size(); k++ )
                {
				    VT_HR_EXIT( m_data[i][j].push_back( src.m_data[i][j][k] ) );
                }
			}
		}

        VT_HR_END()
    }

private:
    // no copy allowed
    COverlappedBlocks(const COverlappedBlocks&);
	COverlappedBlocks &operator=(const COverlappedBlocks&);

private:
    vt::CRect m_rect;
    int m_blockDP2;
    int m_blockOvlDP2;

    int m_blocksW, m_blocksH;
    int m_blockSize;

    int m_blkM, m_ovlM;

    vt::vector<vt::vector<vt::vector<T>>> m_data;
};

