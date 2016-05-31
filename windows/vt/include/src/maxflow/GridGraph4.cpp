/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include "GridGraph4.h"

#include "instances.inl"

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::GridGraph4(ProgressCallback abortCallback, void* callbackData)
    : m_nodes(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::GridGraph4(int sizeX, int sizeY, ProgressCallback abortCallback, 
                       void* callbackData)
    : m_nodes(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
    Allocate(sizeX, sizeY);
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4<NLinkFlowType,TotalFlowType>::~GridGraph4()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4<NLinkFlowType,TotalFlowType>::Allocate(int sizeX, int sizeY)
{
    HRESULT hr;

    DeAllocate();

    m_sizeX = sizeX;
    m_sizeY = sizeY;

    m_nodes =  new(std::nothrow) Node[sizeX*sizeY];
    if (m_nodes == NULL)
    {
        return E_OUTOFMEMORY;
    }

    m_nodeNum = sizeX*sizeY;
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
    m_nodeShifts[1] = sizeX;
    m_nodeShifts[2] = -1;
    m_nodeShifts[3] = -sizeX;

    ActiveQueueReset();

    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4<NLinkFlowType,TotalFlowType>::DeAllocate()
{
    if (m_nodes)
    {
        delete [] m_nodes;
        m_nodes = NULL;
    }
    m_orphanQueue.DeAllocate();
    m_orphanStack.DeAllocate();
}
