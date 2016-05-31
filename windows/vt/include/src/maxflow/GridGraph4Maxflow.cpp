/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include <stdio.h>
#include "GridGraph4.h"

#include "instances.inl"

#define INST_PARENT_EDGES(N, T) \
    const unsigned GridGraph4<N, T>::PARENT_TERMINAL = (unsigned)GridGraph4<N, T>::NeighborhoodSize;\
    const unsigned GridGraph4<N, T>::PARENT_ORPHAN   = (unsigned)GridGraph8<N, T>::NeighborhoodSize+1;\
    const unsigned GridGraph4<N, T>::PARENT_FREE     = (unsigned)GridGraph8<N, T>::NeighborhoodSize+2;

INST_PARENT_EDGES(int,   int)
INST_PARENT_EDGES(short, int)
INST_PARENT_EDGES(float, float)

#undef INST_PARENT_EDGES

// /*static*/ const unsigned GridGraph4::PARENT_TERMINAL = GridGraph4::NeighborhoodSize;
// /*static*/ const unsigned GridGraph4::PARENT_ORPHAN   = GridGraph4::NeighborhoodSize+1;
// /*static*/ const unsigned GridGraph4::PARENT_FREE     = GridGraph4::NeighborhoodSize+2;

#define TREE_CAP(p, q, e_pq, e_qp) ((m_nodes[p].m_tree==0) ? m_nodes[p].m_nCap[e_pq] : m_nodes[q].m_nCap[e_qp])

// Invariants:
// 1. m_nodes[p].m_next is valid <=> p is in the queue
// 2. p is active => p is in the queue
// 3. p is free => p is not active

///////////////////////////////////////////////////////////////////////////////////////
// Queue of active nodes
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::ActiveQueueReset()
{
    m_activeQueueTop = m_activeQueueRear = -1;
#   ifndef _NDEBUG
    m_activeQueueCount = 0;
#   endif
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::ActiveQueueAdd(NodeId p)
{
    if (IS_VALID_ID(m_nodes[p].m_next))
    {
        return; // p is already in the queue
    }
    m_nodes[p].m_next = p; // mark it as being in the queue
    if (IS_VALID_ID(m_activeQueueRear))
    {
        m_nodes[m_activeQueueRear].m_next = p;
    }
    else
    {
        m_activeQueueTop = p;
    }
    m_activeQueueRear = p;
#   ifndef NDEBUG
    ++m_activeQueueCount;
#   endif
}

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4<NLinkFlowType,TotalFlowType>::NodeId GridGraph4<NLinkFlowType,TotalFlowType>::ActiveQueueGetTop()
{
    ASSERT(m_activeQueueCount == 0 || (m_activeQueueCount > 0 && IS_VALID_ID(m_activeQueueTop)));
    return m_activeQueueTop;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::ActiveQueueRemoveTop()
{
    if (m_activeQueueTop == m_activeQueueRear)
    {
        ASSERT(m_activeQueueCount == 1);
        ActiveQueueReset();
    }
    else
    {
        ASSERT(m_activeQueueCount > 1);
        NodeId p = m_activeQueueTop;
        m_activeQueueTop = m_nodes[p].m_next;
        m_nodes[p].m_next = -1;
#       ifndef _NDEBUG
        --m_activeQueueCount;
#       endif
    }
}

//////////////////////////////////////////////////////////////////
//                          NodeQueue                           //
//////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::NodeQueue()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::~NodeQueue()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::Allocate(int size)
{
    ASSERT(size > 0);

    DeAllocate();

    m_start =  new(std::nothrow) NodeId[size];
    if (!m_start)
    {
        return E_OUTOFMEMORY;
    }
    m_end = m_start + size;
    Reset();
    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::Reset()
{
    ASSERT(m_start != NULL);
    m_top = m_rear = m_start;
}

// Note: /analyze warning here is a false positive so disabling it. 
#pragma warning( push )
#pragma warning ( disable : 6011 )

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::EnQueue(NodeId p)
{
    ASSERT(m_start != NULL);
    ASSERT(m_rear >= m_start && m_rear<m_end);

    *m_rear = p;

    m_rear ++;
    if (m_rear == m_end)
    {
        m_rear = m_start;
    }
    if (m_rear == m_top)
    {
        INT_PTR size = m_end - m_start;
        INT_PTR sizeNew = 2*size;
        NodeId* startNew =  new(std::nothrow) NodeId[sizeNew];
        if (startNew == NULL)
        {
            return E_OUTOFMEMORY;
        }
        memcpy(startNew, m_top, (m_end - m_top)*sizeof(NodeId));
        memcpy(startNew + (m_end - m_top), m_start, (m_top - m_start)*sizeof(NodeId));
        delete [] m_start;
        m_start = startNew;
        m_end = m_start + sizeNew;
        m_top = m_start;
        m_rear = m_start + size;
    }

    return S_OK;
}

#pragma warning(pop)

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4<NLinkFlowType,TotalFlowType>::NodeId GridGraph4<NLinkFlowType,TotalFlowType>::NodeQueue::DeQueue()
{
    ASSERT(m_start != NULL);

    if (m_top == m_rear)
    {
        return -1;
    }
    NodeId p = *m_top;
    m_top ++;
    if (m_top == m_end)
    {
        m_top = m_start;
    }
    return p;
}

//////////////////////////////////////////////////////////////////
//                          NodeStack                           //
//////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::NodeStack()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::~NodeStack()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::Allocate(int size)
{
    ASSERT(size > 0);

    DeAllocate();

    m_start =  new(std::nothrow) NodeId[size];
    if (!m_start)
    {
        return E_OUTOFMEMORY;
    }
    m_end = m_start + size;
    Reset();
    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::Reset()
{
    ASSERT(m_start != NULL);
    m_current = m_start;
}
// Note: /analyze warning here is a false positive so disabling it. 
#pragma warning( push )
#pragma warning ( disable : 6011 )

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::Push(NodeId p)
{
    ASSERT(m_start != NULL);
    ASSERT(m_current >= m_start && m_current<=m_end);

    if (m_current == m_end)
    {
        INT_PTR size = m_end - m_start;
        INT_PTR sizeNew = 2*size;
        NodeId* startNew =  new(std::nothrow) NodeId[sizeNew];
        if (startNew == NULL)
        {
            return E_OUTOFMEMORY;
        }
        memcpy(startNew, m_start, size*sizeof(NodeId));
        delete [] m_start;
        m_start = startNew;
        m_end = m_start + sizeNew;
        m_current = m_start + size;
    }

    *m_current ++ = p;

    return S_OK;
}

#pragma warning( pop )

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4<NLinkFlowType,TotalFlowType>::NodeId GridGraph4<NLinkFlowType,TotalFlowType>::NodeStack::Pop()
{
    ASSERT(m_start != NULL);

    if (m_current == m_start)
    {
        return -1;
    }
    m_current --;
    return *m_current;
}

///////////////////////////////////////////////////////////////////////////////////////
// Growth
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  bool GridGraph4<NLinkFlowType,TotalFlowType>::Grow(NodeId& source_node, NodeId& sink_node, int& source_edge)
{
    NodeId p;
    NodeId q;
    int e;

    while ( 1 )
    {
        p = ActiveQueueGetTop();
        if (!IS_VALID_ID(p))
        {
            return false; // no more augmenting paths - done
        }
        if (m_nodes[p].m_isActive)
        {
            if (m_nodes[p].m_tree == 0) // p is in the source set
            {
                for (e=0; e<NeighborhoodSize; e++)
                {
                    if (m_nodes[p].m_nCap[e] == 0)
                    {
                        continue;
                    }
                    q = GetNeib(p, e);

                    if (m_nodes[q].m_parent == PARENT_FREE)
                    {
                        // add q to the source search m_tree as active node
                        m_nodes[q].m_tree = 0;
                        ActiveQueueAdd(q);
                        m_nodes[q].m_isActive = 1;
                        m_nodes[q].m_parent = GetReverseEdge(e);
                    }
                    else if (m_nodes[q].m_tree == 1)
                    {
                        // return path
                        source_node = p;
                        sink_node = q;
                        source_edge = e;
                        return true;
                    }
                }
            }
            else // p is in the sink set
            {
                for (e=0; e<NeighborhoodSize; e++)
                {
                    q = GetAndCheckNeib(p, e);
                    if (!IS_VALID_ID(q))
                    {
                        continue;
                    }
                    if (m_nodes[q].m_nCap[GetReverseEdge(e)] == 0)
                    {
                        continue;
                    }

                    if (m_nodes[q].m_parent == PARENT_FREE)
                    {
                        // add q to the search m_tree as active node
                        m_nodes[q].m_tree = 1;
                        ActiveQueueAdd(q);
                        m_nodes[q].m_isActive = 1;
                        m_nodes[q].m_parent = GetReverseEdge(e);
                    }
                    else if (m_nodes[q].m_tree == 0)
                    {
                        // return path
                        source_node = q;
                        sink_node = p;
                        source_edge = GetReverseEdge(e);
                        return true;
                    }
                }
            }
        }
        ActiveQueueRemoveTop();
    }
}

///////////////////////////////////////////////////////////////////////////////////////
// Augmentation
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::Augment(NodeId source_node, NodeId sink_node, int boundaryEdge)
{
    HRESULT hr;
    NodeId p;
    int e;
    TLinkFlowType delta;

    // compute bottleneck capacity
    p = source_node;
    delta = m_nodes[p].m_nCap[boundaryEdge];
    while ( 1 )
    {
        e = m_nodes[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > m_nodes[p].m_tCap)
            {
                delta = m_nodes[p].m_tCap;
            }
            break;
        }
        p = GetNeib(p, e);
        e = GetReverseEdge(e);
        if (e < NeighborhoodSize && delta > m_nodes[p].m_nCap[e])
        {
            delta = m_nodes[p].m_nCap[e];
        }
    }
    p = sink_node;
    while ( 1 )
    {
        e = m_nodes[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > -m_nodes[p].m_tCap)
            {
                delta = - m_nodes[p].m_tCap;
            }
            break;
        }
        if (e < NeighborhoodSize && delta > m_nodes[p].m_nCap[e])
        {
            delta = m_nodes[p].m_nCap[e];
        }
        p = GetNeib(p, e);
    }

    // augment
    m_nodes[source_node].m_nCap[boundaryEdge] = (NLinkFlowType)(m_nodes[source_node].m_nCap[boundaryEdge] - delta);
    m_nodes[sink_node].m_nCap[GetReverseEdge(boundaryEdge)] = (NLinkFlowType)(m_nodes[sink_node].m_nCap[GetReverseEdge(boundaryEdge)] + delta);
    p = source_node;
    while ( 1 )
    {
        e = m_nodes[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            m_nodes[p].m_tCap = (TLinkFlowType)(m_nodes[p].m_tCap - delta);
            if (m_nodes[p].m_tCap == 0)
            {
                m_nodes[p].m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
        if( e < NeighborhoodSize )
          m_nodes[p].m_nCap[e] = (NLinkFlowType)(m_nodes[p].m_nCap[e] + delta);
        NodeId p_prev = p;
        p = GetNeib(p, e);
        e = GetReverseEdge(e);
        if( e < NeighborhoodSize )
          m_nodes[p].m_nCap[e] = (NLinkFlowType)(m_nodes[p].m_nCap[e] - delta);
        if (e < NeighborhoodSize && m_nodes[p].m_nCap[e] == 0)
        {
            m_nodes[p_prev].m_parent = PARENT_ORPHAN;
            hr = m_orphanStack.Push(p_prev);
            if (hr != S_OK)
            {
                return hr;
            }
        }
    }
    p = sink_node;
    while ( 1 )
    {
        e = m_nodes[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            m_nodes[p].m_tCap = (TLinkFlowType)(m_nodes[p].m_tCap + delta);
            if (-m_nodes[p].m_tCap == 0)
            {
                m_nodes[p].m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
        if( e < NeighborhoodSize )
          m_nodes[p].m_nCap[e] = (NLinkFlowType)(m_nodes[p].m_nCap[e] - delta);
        if (e < NeighborhoodSize && m_nodes[p].m_nCap[e] == 0)
        {
            m_nodes[p].m_parent = PARENT_ORPHAN;
            hr = m_orphanStack.Push(p);
            if (hr != S_OK)
            {
                return hr;
            }
        }
        p = GetNeib(p, e);
        m_nodes[p].m_nCap[GetReverseEdge(e)] = (NLinkFlowType)(m_nodes[p].m_nCap[GetReverseEdge(e)] + delta);
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Adoption
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::Adopt(NodeId p)
{
    HRESULT hr;
    NodeId q;
    NodeId r;
    int e, parent_new = PARENT_ORPHAN;
    int d, d_new = 0;
    unsigned int m_tree = m_nodes[p].m_tree;

    // try to find new valid parent
    for (e=0; e<NeighborhoodSize; e++)
    {
        q = GetAndCheckNeib(p, e);
        if (!IS_VALID_ID(q))
        {
            continue;
        }
        if (m_nodes[q].m_parent==PARENT_FREE || m_nodes[q].m_tree!=m_tree)
        {
            continue;
        }
        if (TREE_CAP(q, p, GetReverseEdge(e), e) == 0)
        {
            continue;
        }

        // check origin
        r = q;
        for (d=0; ; d++)
        {
            if (m_nodes[r].m_parent >= NeighborhoodSize)
            {
                break;
            }
            r = GetNeib(r, m_nodes[r].m_parent);
        }
        if (m_nodes[r].m_parent == PARENT_TERMINAL)
        {
            if (parent_new==PARENT_ORPHAN || d_new>d) // new parent is found
            {
                d_new = d;
                parent_new = e;
            }
        }
    }

    if (parent_new != PARENT_ORPHAN)
    {
        m_nodes[p].m_parent = parent_new;
    }
    else // new parent is not found
    {
        for (e=0; e<NeighborhoodSize; e++)
        {
            q = GetAndCheckNeib(p, e);
            if (!IS_VALID_ID(q))
            {
                continue;
            }
            if (m_nodes[q].m_parent==PARENT_FREE || m_nodes[q].m_tree!=m_tree)
            {
                continue;
            }
            if (TREE_CAP(q, p, GetReverseEdge(e), e) > 0)
            {
                ActiveQueueAdd(q);
                m_nodes[q].m_isActive = 1;
            }
            if (m_nodes[q].m_parent == (unsigned int)GetReverseEdge(e))
            {
                m_nodes[q].m_parent = PARENT_ORPHAN;
                hr = m_orphanQueue.EnQueue(q);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
        }
        m_nodes[p].m_parent = PARENT_FREE;
        m_nodes[p].m_isActive = 0;
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::Maxflow(int flag)
{
    HRESULT hr;
    NodeId p;

    // compute initial flow using dynamic programming
    switch (flag)
    {
        case 0:
            hr = MultipleTreesDynamicProgrammingMaxflowHorz();
            if (hr != S_OK)
            {
                return hr;
            }
            break;
        case 1:
            hr = DynamicProgrammingMaxflowVert();
            if (hr != S_OK)
            {
                return hr;
            }
            hr = DynamicProgrammingMaxflowHorz();
            if (hr != S_OK)
            {
                return hr;
            }
            break;
        case 2:
            hr = DynamicProgrammingMaxflowHorz();
            if (hr != S_OK)
            {
                return hr;
            }
            break;
        default:
            break;
    }

	// now main algorithm

    // initialization
    ActiveQueueReset();
    m_orphanQueue.Reset();
    m_orphanStack.Reset();

    for (p=0; p<m_nodeNum; p++)
    {
        m_nodes[p].m_next = -1;
        if (m_nodes[p].m_tCap != 0)
        {
            ActiveQueueAdd(p);
            m_nodes[p].m_parent = PARENT_TERMINAL;
            m_nodes[p].m_tree = (m_nodes[p].m_tCap > 0) ? 0 : 1;
            m_nodes[p].m_isActive = 1;
        }
        else
        {
            m_nodes[p].m_parent = PARENT_FREE;
            m_nodes[p].m_isActive = 0;
        }
    }

    // main loop
    for( int iteration = 0; ; ++iteration )
    {
        if(m_abortCallback && (iteration % 16 == 0) && m_abortCallback(iteration, m_callbackData))
          return S_FALSE;

        NodeId source_node;
        NodeId sink_node;
        int boundaryEdge;

        // growth
        if (!Grow(source_node, sink_node, boundaryEdge))
        {
            break; // done
        }

        // augmentation
        hr = Augment(source_node, sink_node, boundaryEdge);
        if (hr != S_OK)
        {
            return hr;
        }

        // adoption
        while ( 1 )
        {
            NodeId orphan = m_orphanQueue.DeQueue();
            if (!IS_VALID_ID(orphan))
            {
                orphan = m_orphanStack.Pop();
                if (!IS_VALID_ID(orphan))
                {
                    break;
                }
            }

            hr = Adopt(orphan);
            if (hr != S_OK)
            {
                return hr;
            }
        }
    }

    return S_OK;
}

