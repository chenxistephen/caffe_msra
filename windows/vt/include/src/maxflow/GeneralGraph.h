/**************************************************************************\
*
* Copyright (c) 2004, Microsoft Corp.  All Rights Reserved.
*
* Module Name:
*
*   gridgraph8.h
*
* Abstract:
*
*   Maxflow algorithm for general graphs
*
* Modification History:
*
*   Creation vnk 05/07/04
*
\**************************************************************************/

#pragma once

#include "MaxflowGraph.h"
using namespace GraphcutAPI;

template <class NLinkFlowType, class TotalFlowType> 
  struct GeneralGraph<NLinkFlowType,TotalFlowType>::Node
{
    Edge*           m_firstEdge;

    Node*           m_next; // queue of active nodes
    Edge*           m_parent;
    TLinkFlowType   m_tCap; // if m_tCap > 0, then residual capacity of edge source->node is m_tCap
                            // if m_tCap < 0, then residual capacity of edge node->sink is -m_tCap
    unsigned        m_tree : 1; // 0 = source, 1 = sink (if parent!=PARENT_FREE)
    unsigned        m_isActive : 1;
};

template <class NLinkFlowType, class TotalFlowType> 
  struct GeneralGraph<NLinkFlowType,TotalFlowType>::Edge
{
    Node*           m_to;
    Edge*           m_reverse;
    Edge*           m_next; // next edge with the same origin

    NLinkFlowType   m_nCap;
};

template <class NLinkFlowType, class TotalFlowType> 
  struct GeneralGraph<NLinkFlowType,TotalFlowType>::NodeBlock
{
    NodeBlock()
        : m_dataNum(0)
    {
    }

    static const int nodeBlockSize = 128;
    int         m_dataNum; // number of nodes in m_data
    NodeBlock*  m_next;
    Node      m_data[nodeBlockSize];
};

template <class NLinkFlowType, class TotalFlowType> 
  struct GeneralGraph<NLinkFlowType,TotalFlowType>::EdgeBlock
{
    EdgeBlock()
        : m_dataNum(0)
    {
    }

    static const int edgeBlockSize = 128; // Multiple of 2, since edges allocated in pairs
    int         m_dataNum; // number of edges in m_data
    EdgeBlock*  m_next;
    Edge      m_data[edgeBlockSize];
};


template <class NLinkFlowType, class TotalFlowType> 
  inline typename GeneralGraph<NLinkFlowType,TotalFlowType>::NodeId 
    GeneralGraph<NLinkFlowType,TotalFlowType>::AddNode()
{
    if (!m_nodeBlockFirst || m_nodeBlockFirst->m_dataNum==NodeBlock::nodeBlockSize)
    {
        NodeBlock* block =  new(std::nothrow) NodeBlock;
        if (!block)
        {
            return NULL;
        }
        block -> m_next = m_nodeBlockFirst;
        m_nodeBlockFirst = block;
    }
    Node* p = & m_nodeBlockFirst->m_data[m_nodeBlockFirst->m_dataNum++];
    p -> m_tCap = 0;
    p -> m_firstEdge = NULL;
    return p;
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GeneralGraph<NLinkFlowType,TotalFlowType>::AddTEdge(NodeId p, TLinkFlowType w)
{
    ASSERT(IsAllocated());
    ((Node*)p)->m_tCap = (TLinkFlowType)(((Node*)p)->m_tCap + w);
}

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::AddNEdges(NodeId p, NodeId q, NLinkFlowType w_pq, NLinkFlowType w_qp)
{
    ASSERT(IsAllocated());
    if (!m_edgeBlockFirst || m_edgeBlockFirst->m_dataNum>=EdgeBlock::edgeBlockSize-1)
    {
        EdgeBlock* block =  new(std::nothrow) EdgeBlock;
        if (!block)
        {
            return E_OUTOFMEMORY;
        }
        block -> m_next = m_edgeBlockFirst;
        m_edgeBlockFirst = block;
    }
    Edge* e_pq = & m_edgeBlockFirst->m_data[m_edgeBlockFirst->m_dataNum++];
    Edge* e_qp = & m_edgeBlockFirst->m_data[m_edgeBlockFirst->m_dataNum++];

    e_pq -> m_next = ((Node*)p) -> m_firstEdge;
    e_qp -> m_next = ((Node*)q) -> m_firstEdge;

    ((Node*)p) -> m_firstEdge = e_pq;
    ((Node*)q) -> m_firstEdge = e_qp;

    e_pq -> m_reverse = e_qp;
    e_qp -> m_reverse = e_pq;

    e_pq -> m_to = (Node*)q;
    e_qp -> m_to = (Node*)p;

    e_pq -> m_nCap = w_pq;
    e_qp -> m_nCap = w_qp;

    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
inline int GeneralGraph<NLinkFlowType,TotalFlowType>::GetSegmentation(NodeId p)
{
    ASSERT(IsAllocated());
    return (((Node*)p)->m_parent==PARENT_FREE || ((Node*)p)->m_tree) ? 1 : 0;
}
