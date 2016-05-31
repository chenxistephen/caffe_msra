/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include <math.h>
#include <assert.h>
#include "GeneralGraph.h"

#include "instances.inl"

/*static*/
#define INST_PARENT_EDGES(N, T) \
    GeneralGraph<N, T>::Edge* GeneralGraph<N, T>::PARENT_TERMINAL = (GeneralGraph<N, T>::Edge*)0;\
    GeneralGraph<N, T>::Edge* GeneralGraph<N, T>::PARENT_ORPHAN   = (GeneralGraph<N, T>::Edge*)1;\
    GeneralGraph<N, T>::Edge* GeneralGraph<N, T>::PARENT_FREE     = (GeneralGraph<N, T>::Edge*)2;

INST_PARENT_EDGES(int,   int)
INST_PARENT_EDGES(short, int)
INST_PARENT_EDGES(float, float)

#undef INST_PARENT_EDGES

//#define   PARENT_TERMINAL ((Edge*)0)
//#define PARENT_ORPHAN   ((Edge*)1)
//#define PARENT_FREE     ((Edge*)2)

#define TREE_CAP(p, q, e_pq, e_qp) (((p)->m_tree==0) ? (e_pq)->m_nCap : (e_qp)->m_nCap)

// Invariants:
// 1. p->m_next is valid <=> p is in the queue
// 2. p is active => p is in the queue
// 3. p is free => p is not active

///////////////////////////////////////////////////////////////////////////////////////
// Queue of active nodes
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::ActiveQueueReset()
{
    m_activeQueueTop = m_activeQueueRear = NULL;
#   ifndef _NDEBUG
    m_activeQueueCount = 0;
#   endif
}

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::ActiveQueueAdd(Node* p)
{
    if (p->m_next)
    {
        return; // p is already in the queue
    }
    p -> m_next = p; // mark it as being in the queue
    if (m_activeQueueRear)
    {
        m_activeQueueRear -> m_next = p;
    }
    else
    {
        m_activeQueueTop = p;
    }
    m_activeQueueRear = p;
#   ifndef _NDEBUG
    ++m_activeQueueCount;
#   endif
}

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GeneralGraph<NLinkFlowType,TotalFlowType>::Node* GeneralGraph<NLinkFlowType,TotalFlowType>::ActiveQueueGetTop()
{
    ASSERT(m_activeQueueCount == 0 || (m_activeQueueCount > 0 && m_activeQueueTop));
    return m_activeQueueTop;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::ActiveQueueRemoveTop()
{
    if (m_activeQueueTop == m_activeQueueRear)
    {
        ASSERT(m_activeQueueCount == 1);
        ActiveQueueReset();
    }
    else
    {
        ASSERT(m_activeQueueCount > 1);
        Node* p = m_activeQueueTop;
        m_activeQueueTop = p -> m_next;
        p -> m_next = NULL;
#       ifndef _NDEBUG
        --m_activeQueueCount;
#       endif
    }
}

//////////////////////////////////////////////////////////////////
//                          NodeQueue                           //
//////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::NodeQueue()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::~NodeQueue()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::Allocate(int size)
{
    ASSERT(size > 0);

    DeAllocate();

    m_start =  new(std::nothrow) Node*[size];
    if (!m_start)
    {
        return E_OUTOFMEMORY;
    }
    m_end = m_start + size;
    Reset();
    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

template <class NLinkFlowType,class TotalFlowType> 
  inline void GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::Reset()
{
    ASSERT(m_start != NULL);
    m_top = m_rear = m_start;
}

// Note: /analyze warning here is a false positive so disabling it. 
#pragma warning( push )
#pragma warning ( disable : 6011 )

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::EnQueue(Node* p)
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
        Node** startNew =  new(std::nothrow) Node*[sizeNew];
        if (startNew == NULL)
        {
            return E_OUTOFMEMORY;
        }
        memcpy(startNew, m_top, (m_end - m_top)*sizeof(Node*));
        memcpy(startNew + (m_end - m_top), m_start, (m_top - m_start)*sizeof(Node*));
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
  inline typename GeneralGraph<NLinkFlowType,TotalFlowType>::Node* GeneralGraph<NLinkFlowType,TotalFlowType>::NodeQueue::DeQueue()
{
    ASSERT(m_start != NULL);

    if (m_top == m_rear)
    {
        return NULL;
    }
    Node* p = *m_top;
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
  GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::NodeStack()
    : m_start(NULL)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::~NodeStack()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::Allocate(int size)
{
    ASSERT(size > 0);

    DeAllocate();

    m_start =  new(std::nothrow) Node*[size];
    if (!m_start)
    {
        return E_OUTOFMEMORY;
    }
    m_end = m_start + size;
    Reset();
    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::DeAllocate()
{
    if (m_start)
    {
        delete [] m_start;
        m_start = NULL;
    }
}

// Note: /analyze warning here is a false positive so disabling it. 
#pragma warning( push )
#pragma warning ( disable : 6011 )

template <class NLinkFlowType,class TotalFlowType> 
  inline void GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::Reset()
{
    ASSERT(m_start != NULL);
    m_current = m_start;
}

template <class NLinkFlowType,class TotalFlowType> 
  inline HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::Push(Node* p)
{
    ASSERT(m_start != NULL);
    ASSERT(m_current >= m_start && m_current<=m_end);

    if (m_current == m_end)
    {
        INT_PTR size = m_end - m_start;
        INT_PTR sizeNew = 2*size;
        Node** startNew =  new(std::nothrow) Node*[sizeNew];
        if (startNew == NULL)
        {
            return E_OUTOFMEMORY;
        }
        memcpy(startNew, m_start, size*sizeof(Node*));
        delete [] m_start;
        m_start = startNew;
        m_end = m_start + sizeNew;
        m_current = m_start + size;
    }

    *m_current ++ = p;

    return S_OK;
}

#pragma warning(pop)

template <class NLinkFlowType,class TotalFlowType> 
  inline typename GeneralGraph<NLinkFlowType,TotalFlowType>::Node* GeneralGraph<NLinkFlowType,TotalFlowType>::NodeStack::Pop()
{
    ASSERT(m_start != NULL);

    if (m_current == m_start)
    {
        return NULL;
    }
    m_current --;
    return *m_current;
}

///////////////////////////////////////////////////////////////////////////////////////
// Growth
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  typename GeneralGraph<NLinkFlowType,TotalFlowType>::Edge* GeneralGraph<NLinkFlowType,TotalFlowType>::Grow()
{
    Node* p;
    Node* q;
    Edge* e;

    while ( 1 )
    {
        p = ActiveQueueGetTop();
        if (p == NULL)
        {
            return false; // no more augmenting paths - done
        }
        if (p->m_isActive)
        {
            if (p->m_tree == 0) // p is in the source set
            {
                for (e=p->m_firstEdge; e; e=e->m_next)
                {
                    if (e->m_nCap == 0)
                    {
                        continue;
                    }
                    q = e -> m_to;

                    if (q->m_parent == PARENT_FREE)
                    {
                        // add q to the source search m_tree as active node
                        q -> m_tree = 0;
                        ActiveQueueAdd(q);
                        q -> m_isActive = 1;
                        q -> m_parent = e -> m_reverse;
                    }
                    else if (q->m_tree == 1)
                    {
                        // return path
                        return e;
                    }
                }
            }
            else // p is in the sink set
            {
                for (e=p->m_firstEdge; e; e=e->m_next)
                {
                    if (e->m_reverse->m_nCap == 0)
                    {
                        continue;
                    }
                    q = e -> m_to;

                    if (q->m_parent == PARENT_FREE)
                    {
                        // add q to the search m_tree as active node
                        q -> m_tree = 1;
                        ActiveQueueAdd(q);
                        q -> m_isActive = 1;
                        q -> m_parent = e -> m_reverse;
                    }
                    else if (q->m_tree == 0)
                    {
                        // return path
                        return e -> m_reverse;
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
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::Augment(Edge* boundaryEdge)
{
    HRESULT hr;
    Node *p;
    Edge* e;
    TLinkFlowType delta;

    // compute bottleneck capacity
    p = boundaryEdge -> m_reverse -> m_to;
    delta = boundaryEdge -> m_nCap;
    while ( 1 )
    {
        e = p -> m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > p->m_tCap)
            {
                delta = p -> m_tCap;
            }
            break;
        }
        if (delta > e->m_reverse->m_nCap)
        {
            delta = e -> m_reverse -> m_nCap;
        }
        p = e -> m_to;
    }
    p = boundaryEdge -> m_to;
    while ( 1 )
    {
        e = p -> m_parent;
        if (e == PARENT_TERMINAL)
        {
            if (delta > -p->m_tCap)
            {
                delta = - p -> m_tCap;
            }
            break;
        }
        if (delta > e->m_nCap)
        {
            delta = e -> m_nCap;
        }
        p = e -> m_to;
    }

    // augment
    boundaryEdge->m_nCap = (NLinkFlowType)(boundaryEdge->m_nCap - delta);
    boundaryEdge->m_reverse->m_nCap = (NLinkFlowType)(boundaryEdge->m_reverse->m_nCap + delta);
    p = boundaryEdge -> m_reverse -> m_to;
    while ( 1 )
    {
        e = p -> m_parent;
        if (e == PARENT_TERMINAL)
        {
            p->m_tCap = (TLinkFlowType)(p->m_tCap - delta);
            if (p->m_tCap == 0)
            {
                p -> m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
        e->m_nCap = (NLinkFlowType)(e->m_nCap + delta);
        e->m_reverse->m_nCap = (NLinkFlowType)(e->m_reverse->m_nCap - delta);
        if (e->m_reverse->m_nCap == 0)
        {
            p -> m_parent = PARENT_ORPHAN;
            hr = m_orphanStack.Push(p);
            if (hr != S_OK)
            {
                return hr;
            }
        }
        p = e -> m_to;
    }
    p = boundaryEdge -> m_to;
    while ( 1 )
    {
        e = p -> m_parent;
        if (e == PARENT_TERMINAL)
        {
            p->m_tCap = (TLinkFlowType)(p->m_tCap + delta);
            if (-p->m_tCap == 0)
            {
                p -> m_parent = PARENT_ORPHAN;
                hr = m_orphanStack.Push(p);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
            break;
        }
        e->m_reverse->m_nCap = (NLinkFlowType)(e->m_reverse->m_nCap + delta);
        e->m_nCap = (NLinkFlowType)(e->m_nCap - delta);
        if (e->m_nCap == 0)
        {
            p -> m_parent = PARENT_ORPHAN;
            hr = m_orphanStack.Push(p);
            if (hr != S_OK)
            {
                return hr;
            }
        }
        p = e -> m_to;
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Adoption
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::Adopt(Node* p)
{
    HRESULT hr;
    Node* q;
    Node* r;
    Edge* e;
    Edge* parent_new = PARENT_ORPHAN;
    int d, d_new = 0;
    unsigned int m_tree = p -> m_tree;

    // try to find new valid parent
    for (e=p->m_firstEdge; e; e=e->m_next)
    {
        q = e -> m_to;
        if (q->m_parent==PARENT_FREE || q->m_tree!=m_tree)
        {
            continue;
        }
        if (TREE_CAP(q, p, e->m_reverse, e) == 0)
        {
            continue;
        }

        // check origin
        r = q;
        for (d=0; ; d++)
        {
            if (r->m_parent == PARENT_TERMINAL || r->m_parent == PARENT_ORPHAN)
            {
                break;
            }
            r = r -> m_parent -> m_to;
        }
        if (r->m_parent == PARENT_TERMINAL)
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
        p -> m_parent = parent_new;
    }
    else // new parent is not found
    {
        for (e=p->m_firstEdge; e; e=e->m_next)
        {
            q = e -> m_to;
            if (q->m_parent==PARENT_FREE || q->m_tree!=m_tree)
            {
                continue;
            }
            if (TREE_CAP(q, p, e->m_reverse, e) > 0)
            {
                ActiveQueueAdd(q);
                q -> m_isActive = 1;
            }
            if (q->m_parent == e->m_reverse)
            {
                q -> m_parent = PARENT_ORPHAN;
                hr = m_orphanQueue.EnQueue(q);
                if (hr != S_OK)
                {
                    return hr;
                }
            }
        }
        p -> m_parent = PARENT_FREE;
        p -> m_isActive = 0;
    }

    return S_OK;
}

///////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::Maxflow()
{
    HRESULT hr;
    NodeBlock* block;
    Node* p;

    // initialization
    ActiveQueueReset();
    m_orphanQueue.Reset();
    m_orphanStack.Reset();

    for (block=m_nodeBlockFirst; block; block=block->m_next)
    {
        for (p=&block->m_data[0]; p<&block->m_data[block->m_dataNum]; p++)
        {
            p -> m_next = NULL;
            if (p->m_tCap != 0)
            {
                ActiveQueueAdd(p);
                p -> m_parent = PARENT_TERMINAL;
                p -> m_tree = (p->m_tCap > 0) ? 0 : 1;
                p -> m_isActive = 1;
            }
            else
            {
                p -> m_parent = PARENT_FREE;
                p -> m_isActive = 0;
            }
        }
    }

    // main loop
    for( int iteration = 1; ; ++iteration )
    {
        if(m_abortCallback && (iteration % 16 == 0) && m_abortCallback(iteration, m_callbackData))
          return S_FALSE;

        // growth
        Edge* boundaryEdge = Grow();
        if (boundaryEdge == NULL)
        {
            break; // done
        }

        // augmentation
        hr = Augment(boundaryEdge);
        if (hr != S_OK)
        {
            return hr;
        }

        // adoption
        while ( 1 )
        {
            Node* orphan = m_orphanQueue.DeQueue();
            if (orphan == NULL)
            {
                orphan = m_orphanStack.Pop();
                if (orphan == NULL)
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
