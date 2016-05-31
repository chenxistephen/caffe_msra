/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include <stdio.h>
#include "GridGraph4Plus.h"

#include "instances.inl"

#define INST_PARENT_EDGES(N, T) \
    const unsigned GridGraph4Plus<N, T>::PARENT_TERMINAL  = (unsigned)GridGraph4Plus<N, T>::NeighborhoodSize;\
    const unsigned GridGraph4Plus<N, T>::PARENT_ORPHAN    = (unsigned)GridGraph4Plus<N, T>::NeighborhoodSize+1;\
    const unsigned GridGraph4Plus<N, T>::PARENT_FREE      = (unsigned)GridGraph4Plus<N, T>::NeighborhoodSize+2;\
    const unsigned GridGraph4Plus<N, T>::PARENT_EXTRA_ARC = (unsigned)GridGraph4Plus<N, T>::NeighborhoodSize+3;

INST_PARENT_EDGES(int,   int)
INST_PARENT_EDGES(short, int)
INST_PARENT_EDGES(float, float)

#undef INST_PARENT_EDGES

// /*static*/ const unsigned GridGraph4Plus::PARENT_TERMINAL = GridGraph4Plus::NeighborhoodSize;
// /*static*/ const unsigned GridGraph4Plus::PARENT_ORPHAN   = GridGraph4Plus::NeighborhoodSize+1;
// /*static*/ const unsigned GridGraph4Plus::PARENT_FREE     = GridGraph4Plus::NeighborhoodSize+2;

#define TREE_CAP(p, q, e_pq, e_qp) ((m_nodesArray[p].m_tree==0) ? m_nodesArray[p].m_nCap[e_pq] : m_nodesArray[q].m_nCap[e_qp])
#define TREE_CAP_EXTRA(p, q, a_pq, a_qp) ((m_nodesArray[p].m_tree==0) ? (a_pq)->m_nCap : (a_qp)->m_nCap)

// Invariants:
// 1. p->m_next is valid <=> p is in the queue
// 2. p is active => p is in the queue
// 3. p is free => p is not active

///////////////////////////////////////////////////////////////////////////////////////
// Queue of active nodes
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::ActiveQueueReset()
{
    m_activeQueueTop = m_activeQueueRear = -1;
#   ifndef _NDEBUG
    m_activeQueueCount = 0;
#   endif
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::ActiveQueueAdd(NodeId p)
{
    if (IS_VALID_ID(m_nodesArray[p].m_next))
    {
        return; // p is already in the queue
    }
    m_nodesArray[p].m_next = p; // mark it as being in the queue
    if (IS_VALID_ID(m_activeQueueRear))
    {
        m_nodesArray[m_activeQueueRear].m_next = p;
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
  inline typename GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeId GridGraph4Plus<NLinkFlowType,TotalFlowType>::ActiveQueueGetTop()
{
    ASSERT(m_activeQueueCount == 0 || (m_activeQueueCount > 0 && IS_VALID_ID(m_activeQueueTop)));
    return m_activeQueueTop;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::ActiveQueueRemoveTop()
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
        m_activeQueueTop = m_nodesArray[p].m_next;
        m_nodesArray[p].m_next = -1;
#       ifndef _NDEBUG
        --m_activeQueueCount;
#       endif
    }
}

//////////////////////////////////////////////////////////////////
//                          NodeQueue                           //
//////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::NodeQueue()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::~NodeQueue()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::Allocate(int size)
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
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::Reset()
{
    ASSERT(m_start != NULL);
    m_top = m_rear = m_start;
}

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::EnQueue(NodeId p)
{
    ASSERT(m_start != NULL);
    ASSERT(m_rear >= m_start && m_rear<m_end);
	__analysis_assume(m_rear != NULL);

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

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeId GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeQueue::DeQueue()
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
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::NodeStack()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::~NodeStack()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::Allocate(int size)
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
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::Reset()
{
    ASSERT(m_start != NULL);
    m_current = m_start;
}

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::Push(NodeId p)
{
    ASSERT(m_start != NULL);
    ASSERT(m_current >= m_start && m_current<=m_end);
	__analysis_assume(m_current != NULL);

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

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeId GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeStack::Pop()
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
  bool GridGraph4Plus<NLinkFlowType,TotalFlowType>::Grow(NodeId& source_node, NodeId& sink_node, int& source_edge)
{
    NodeId p;
    NodeId q;
	ExtraArc* a;
    int e;

    while ( 1 )
    {
        p = ActiveQueueGetTop();
        if (!IS_VALID_ID(p))
        {
            return false; // no more augmenting paths - done
        }
        if (m_nodesArray[p].m_isActive)
        {
            if (m_nodesArray[p].m_tree == 0) // p is in the source set
            {
                for (e=0; e<NeighborhoodSize; e++)
                {
                    if (m_nodesArray[p].m_nCap[e] == 0)
                    {
                        continue;
                    }
                    q = GetNeib(p, e);

                    if (m_nodesArray[q].m_parent == PARENT_FREE)
                    {
                        // add q to the source search m_tree as active node
                        m_nodesArray[q].m_tree = 0;
                        ActiveQueueAdd(q);
                        m_nodesArray[q].m_isActive = 1;
                        m_nodesArray[q].m_parent = GetReverseEdge(e);
                    }
                    else if (m_nodesArray[q].m_tree == 1)
                    {
                        // return path
                        source_node = p;
                        sink_node = q;
                        source_edge = e;
                        return true;
                    }
                }
				for (e=PARENT_EXTRA_ARC, a=m_extraArcs+m_nodesArray[p].m_extraArcs; a<m_extraArcs+m_nodesArray[p+1].m_extraArcs; e++, a++)
                {
                    if (a->m_nCap == 0)
                    {
                        continue;
                    }
					q = a->m_head;

                    if (m_nodesArray[q].m_parent == PARENT_FREE)
                    {
                        // add q to the source search m_tree as active node
                        m_nodesArray[q].m_tree = 0;
                        ActiveQueueAdd(q);
                        m_nodesArray[q].m_isActive = 1;
						m_nodesArray[q].m_parent = a->m_sister - m_nodesArray[q].m_extraArcs + PARENT_EXTRA_ARC;
                    }
                    else if (m_nodesArray[q].m_tree == 1)
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
                    q = GetNeib(p, e);
                    if (m_nodesArray[q].m_nCap[GetReverseEdge(e)] == 0)
                    {
                        continue;
                    }

                    if (m_nodesArray[q].m_parent == PARENT_FREE)
                    {
                        // add q to the search m_tree as active node
                        m_nodesArray[q].m_tree = 1;
                        ActiveQueueAdd(q);
                        m_nodesArray[q].m_isActive = 1;
                        m_nodesArray[q].m_parent = GetReverseEdge(e);
                    }
                    else if (m_nodesArray[q].m_tree == 0)
                    {
                        // return path
                        source_node = q;
                        sink_node = p;
                        source_edge = GetReverseEdge(e);
                        return true;
                    }
                }
				for (e=PARENT_EXTRA_ARC, a=m_extraArcs+m_nodesArray[p].m_extraArcs; a<m_extraArcs+m_nodesArray[p+1].m_extraArcs; e++, a++)
                {
                    q = a->m_head;
					if (m_extraArcs[a->m_sister].m_nCap == 0)
                    {
                        continue;
                    }

                    if (m_nodesArray[q].m_parent == PARENT_FREE)
                    {
                        // add q to the search m_tree as active node
                        m_nodesArray[q].m_tree = 1;
                        ActiveQueueAdd(q);
                        m_nodesArray[q].m_isActive = 1;
						m_nodesArray[q].m_parent = a->m_sister - m_nodesArray[q].m_extraArcs + PARENT_EXTRA_ARC;
                    }
                    else if (m_nodesArray[q].m_tree == 0)
                    {
                        // return path
                        source_node = q;
                        sink_node = p;
						source_edge = a->m_sister - m_nodesArray[q].m_extraArcs + PARENT_EXTRA_ARC;
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
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::Augment(NodeId source_node, NodeId sink_node, int boundaryEdge)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// I noticed that somebody has added the check "if (e < NeighborhoodSize) ..." to my original code.  //
	// I removed this check since Node::m_parent for all nodes in the path are guaranteed to be          //
	// neither PARENT_ORPHAN nor PARENT_FREE - Vladimir Kolmogorov, 05/2009                              //
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	HRESULT hr;
    NodeId p;
    int e;
    TLinkFlowType delta;

    // compute bottleneck capacity
    p = source_node;
	delta = (boundaryEdge < NeighborhoodSize) ? m_nodesArray[p].m_nCap[boundaryEdge] : m_extraArcs[m_nodesArray[p].m_extraArcs+boundaryEdge-PARENT_EXTRA_ARC].m_nCap;
    while ( 1 )
    {
        e = m_nodesArray[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > m_nodesArray[p].m_tCap)
            {
                delta = m_nodesArray[p].m_tCap;
            }
            break;
        }
		if (e < PARENT_EXTRA_ARC)
		{
			__analysis_assume(e < NeighborhoodSize);
			p = GetNeib(p, e);
			e = GetReverseEdge(e);
			if (delta > m_nodesArray[p].m_nCap[e])
			{
				delta = m_nodesArray[p].m_nCap[e];
			}
		}
		else
		{
			ExtraArc* a = &m_extraArcs[m_nodesArray[p].m_extraArcs+e-PARENT_EXTRA_ARC];
			p = a->m_head;
			if (delta > m_extraArcs[a->m_sister].m_nCap)
			{
				delta = m_extraArcs[a->m_sister].m_nCap;
			}
		}
    }
    p = sink_node;
    while ( 1 )
    {
        e = m_nodesArray[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > -m_nodesArray[p].m_tCap)
            {
                delta = - m_nodesArray[p].m_tCap;
            }
            break;
        }
		if (e < PARENT_EXTRA_ARC)
		{
			__analysis_assume(e < NeighborhoodSize);
			if (delta > m_nodesArray[p].m_nCap[e])
			{
				delta = m_nodesArray[p].m_nCap[e];
			}
			p = GetNeib(p, e);
		}
		else
		{
			ExtraArc* a = &m_extraArcs[m_nodesArray[p].m_extraArcs+e-PARENT_EXTRA_ARC];
			if (delta > a->m_nCap)
			{
				delta = a -> m_nCap;
			}
			p = a->m_head;
		}
    }

    // augment
	if (boundaryEdge < NeighborhoodSize)
	{
	    m_nodesArray[source_node].m_nCap[boundaryEdge] = (NLinkFlowType)(m_nodesArray[source_node].m_nCap[boundaryEdge] - delta);
		m_nodesArray[sink_node].m_nCap[GetReverseEdge(boundaryEdge)] = (NLinkFlowType)(m_nodesArray[sink_node].m_nCap[GetReverseEdge(boundaryEdge)] + delta);
	}
	else
	{
			ExtraArc* a = &m_extraArcs[m_nodesArray[source_node].m_extraArcs+boundaryEdge-PARENT_EXTRA_ARC];
		    a->m_nCap = (NLinkFlowType)(a->m_nCap - delta);
			m_extraArcs[a->m_sister].m_nCap = (NLinkFlowType)(m_extraArcs[a->m_sister].m_nCap + delta);
	}
    p = source_node;
    while ( 1 )
    {
        e = m_nodesArray[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            m_nodesArray[p].m_tCap = (TLinkFlowType)(m_nodesArray[p].m_tCap - delta);
            if (m_nodesArray[p].m_tCap == 0)
            {
                m_nodesArray[p].m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
        NodeId p_prev = p;
        if (e < PARENT_EXTRA_ARC)
		{
			__analysis_assume(e < NeighborhoodSize);
			m_nodesArray[p].m_nCap[e] = (NLinkFlowType)(m_nodesArray[p].m_nCap[e] + delta);
			p = GetNeib(p, e);
			e = GetReverseEdge(e);
			m_nodesArray[p].m_nCap[e] = (NLinkFlowType)(m_nodesArray[p].m_nCap[e] - delta);
			if (m_nodesArray[p].m_nCap[e] != 0) continue;
		}
		else
		{
			ExtraArc* a = &m_extraArcs[m_nodesArray[p].m_extraArcs+e-PARENT_EXTRA_ARC];
			a -> m_nCap = (NLinkFlowType)(a -> m_nCap + delta);
			p = a->m_head;
			a = &m_extraArcs[a->m_sister];
			a -> m_nCap = (NLinkFlowType)(a -> m_nCap - delta);
			if (a->m_nCap != 0) continue;
		}
		m_nodesArray[p_prev].m_parent = PARENT_ORPHAN;
		hr = m_orphanStack.Push(p_prev);
		if (hr != S_OK)
		{
			return hr;
		}
    }
    p = sink_node;
    while ( 1 )
    {
        e = m_nodesArray[p].m_parent;
        if (e == PARENT_TERMINAL)
        {
            m_nodesArray[p].m_tCap = (TLinkFlowType)(m_nodesArray[p].m_tCap + delta);
            if (-m_nodesArray[p].m_tCap == 0)
            {
                m_nodesArray[p].m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
		if (e < PARENT_EXTRA_ARC)
		{
			__analysis_assume(e < NeighborhoodSize);
			m_nodesArray[p].m_nCap[e] = (NLinkFlowType)(m_nodesArray[p].m_nCap[e] - delta);
			if (m_nodesArray[p].m_nCap[e] == 0)
			{
				m_nodesArray[p].m_parent = PARENT_ORPHAN;
				hr = m_orphanStack.Push(p);
				if (hr != S_OK)
				{
					return hr;
				}
			}
			p = GetNeib(p, e);
			m_nodesArray[p].m_nCap[GetReverseEdge(e)] = (NLinkFlowType)(m_nodesArray[p].m_nCap[GetReverseEdge(e)] + delta);
		}
		else
		{
			ExtraArc* a = &m_extraArcs[m_nodesArray[p].m_extraArcs+e-PARENT_EXTRA_ARC];
			a -> m_nCap = (NLinkFlowType)(a -> m_nCap - delta);
			if (a->m_nCap == 0)
			{
				m_nodesArray[p].m_parent = PARENT_ORPHAN;
				hr = m_orphanStack.Push(p);
				if (hr != S_OK)
				{
					return hr;
				}
			}
			p = a->m_head;
			a = &m_extraArcs[a->m_sister];
			a -> m_nCap = (NLinkFlowType)(a -> m_nCap + delta);
		}
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Adoption
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::Adopt(NodeId p)
{
    HRESULT hr;
    NodeId q;
    NodeId r;
	ExtraArc* a;
    int e, parent_new = PARENT_ORPHAN;
    int d, d_new = 0;
    unsigned int m_tree = m_nodesArray[p].m_tree;

    // try to find new valid parent
    for (e=0; e<NeighborhoodSize; e++)
    {
        q = GetNeib(p, e);
        if (m_nodesArray[q].m_parent==PARENT_FREE || m_nodesArray[q].m_tree!=m_tree)
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
			if (m_nodesArray[r].m_parent < NeighborhoodSize)
			{
				r = GetNeib(r, m_nodesArray[r].m_parent);
			}
			else if (m_nodesArray[r].m_parent >= PARENT_EXTRA_ARC)
			{
				r = m_extraArcs[m_nodesArray[r].m_extraArcs + m_nodesArray[r].m_parent - PARENT_EXTRA_ARC].m_head;
			}
			else
            {
                break;
            }
        }
        if (m_nodesArray[r].m_parent == PARENT_TERMINAL)
        {
            if (parent_new==PARENT_ORPHAN || d_new>d) // new parent is found
            {
                d_new = d;
                parent_new = e;
            }
        }
    }
	// now try extra arcs
	for (e=PARENT_EXTRA_ARC, a=m_extraArcs+m_nodesArray[p].m_extraArcs; a<m_extraArcs+m_nodesArray[p+1].m_extraArcs; e++, a++)
	{
		q = a->m_head;
		if (m_nodesArray[q].m_parent==PARENT_FREE || m_nodesArray[q].m_tree!=m_tree)
		{
			continue;
		}
		if (TREE_CAP_EXTRA(q, p, &m_extraArcs[a->m_sister], a) == 0)
		{
			continue;
		}

		// check origin
		r = q;
		for (d=0; ; d++)
		{
			if (m_nodesArray[r].m_parent < NeighborhoodSize)
			{
				r = GetNeib(r, m_nodesArray[r].m_parent);
			}
			else if (m_nodesArray[r].m_parent >= PARENT_EXTRA_ARC)
			{
				r = m_extraArcs[m_nodesArray[r].m_extraArcs + m_nodesArray[r].m_parent - PARENT_EXTRA_ARC].m_head;
			}
			else
			{
				break;
			}
		}
		if (m_nodesArray[r].m_parent == PARENT_TERMINAL)
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
        m_nodesArray[p].m_parent = parent_new;
    }
    else // new parent is not found
    {
        for (e=0; e<NeighborhoodSize; e++)
        {
            q = GetNeib(p, e);
            if (m_nodesArray[q].m_parent==PARENT_FREE || m_nodesArray[q].m_tree!=m_tree)
            {
                continue;
            }
            if (TREE_CAP(q, p, GetReverseEdge(e), e) > 0)
            {
                ActiveQueueAdd(q);
                m_nodesArray[q].m_isActive = 1;
            }
            if (m_nodesArray[q].m_parent == (unsigned int)GetReverseEdge(e))
            {
                m_nodesArray[q].m_parent = PARENT_ORPHAN;
                hr = m_orphanQueue.EnQueue(q);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
        }
		for (e=PARENT_EXTRA_ARC, a=m_extraArcs+m_nodesArray[p].m_extraArcs; a<m_extraArcs+m_nodesArray[p+1].m_extraArcs; e++, a++)
        {
			q = a->m_head;
            if (m_nodesArray[q].m_parent==PARENT_FREE || m_nodesArray[q].m_tree!=m_tree)
            {
                continue;
            }
			if (TREE_CAP_EXTRA(q, p, &m_extraArcs[a->m_sister], a) > 0)
            {
                ActiveQueueAdd(q);
                m_nodesArray[q].m_isActive = 1;
            }
			if (m_nodesArray[q].m_extraArcs + m_nodesArray[q].m_parent - PARENT_EXTRA_ARC == a->m_sister)
            {
                m_nodesArray[q].m_parent = PARENT_ORPHAN;
                hr = m_orphanQueue.EnQueue(q);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
        }
        m_nodesArray[p].m_parent = PARENT_FREE;
        m_nodesArray[p].m_isActive = 0;
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::Maxflow()
{
    HRESULT hr;
    NodeId p;

    // initialization

	if (m_extraEdges0)
	{
		ASSERT(!m_extraArcs);
		hr = InitExtraArcs();
        if (hr != S_OK)
        {
            return hr;
        }
	}

    ActiveQueueReset();
    m_orphanQueue.Reset();
    m_orphanStack.Reset();

    for (p=0; p<m_allocatedNodeNum; p++)
    {
        m_nodesArray[p].m_next = -1;
        if (m_nodesArray[p].m_tCap != 0)
        {
            ActiveQueueAdd(p);
            m_nodesArray[p].m_parent = PARENT_TERMINAL;
            m_nodesArray[p].m_tree = (m_nodesArray[p].m_tCap > 0) ? 0 : 1;
            m_nodesArray[p].m_isActive = 1;
        }
        else
        {
			// dummy nodes will also be here 
            m_nodesArray[p].m_parent = PARENT_FREE;
            m_nodesArray[p].m_isActive = 0;
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

