/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include "GridGraph8.h"

#include "instances.inl"

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph8<NLinkFlowType,TotalFlowType>::GridGraph8(ProgressCallback abortCallback, void* callbackData)
    : m_nodes(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph8<NLinkFlowType,TotalFlowType>::GridGraph8(int sizeX, int sizeY, ProgressCallback abortCallback, 
                       void* callbackData)
    : m_nodes(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
    Allocate(sizeX, sizeY);
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph8<NLinkFlowType,TotalFlowType>::~GridGraph8()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph8<NLinkFlowType,TotalFlowType>::Allocate(int sizeX, int sizeY)
{
    HRESULT hr;

    DeAllocate();

    m_sizeX = sizeX;
    m_sizeY = sizeY;
	m_nodeNum = sizeX*sizeY;

    m_nodes =  new(std::nothrow) Node[m_nodeNum];
    if (m_nodes == NULL)
    {
        return E_OUTOFMEMORY;
    }

    memset(m_nodes, 0, m_nodeNum*sizeof(Node));

    hr = m_orphanQueue.Allocate();
    if (hr != S_OK)
    {
        return hr;
    }
    hr = m_orphanStack.Allocate();
    if (hr != S_OK)
    {
        return hr;
    }

    m_nodeShifts[0] = 1;
    m_nodeShifts[1] = sizeX+1;
    m_nodeShifts[2] = sizeX;
    m_nodeShifts[3] = sizeX-1;
    m_nodeShifts[4] = -1;
    m_nodeShifts[5] = -sizeX-1;
    m_nodeShifts[6] = -sizeX;
    m_nodeShifts[7] = -sizeX+1;

    ActiveQueueReset();

    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph8<NLinkFlowType,TotalFlowType>::DeAllocate()
{
    if (m_nodes)
    {
        delete [] m_nodes;
        m_nodes = NULL;
    }
    m_orphanQueue.DeAllocate();
    m_orphanStack.DeAllocate();
}
