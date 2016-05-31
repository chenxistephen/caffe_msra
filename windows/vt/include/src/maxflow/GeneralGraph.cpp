/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include "GeneralGraph.h"

#include "instances.inl"

template <class NLinkFlowType,class TotalFlowType> 
  GeneralGraph<NLinkFlowType,TotalFlowType>::GeneralGraph(ProgressCallback abortCallback, void* callbackData)
	: m_isAllocated(false),
	  m_nodeBlockFirst(NULL),
	  m_edgeBlockFirst(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
	Allocate();
    ActiveQueueReset();
}

template <class NLinkFlowType,class TotalFlowType> 
  GeneralGraph<NLinkFlowType,TotalFlowType>::~GeneralGraph()
{
	DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GeneralGraph<NLinkFlowType,TotalFlowType>::Allocate()
{
	HRESULT hr;

	DeAllocate();

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

	m_isAllocated = true;

	return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GeneralGraph<NLinkFlowType,TotalFlowType>::DeAllocate()
{
	while (m_nodeBlockFirst)
	{
		NodeBlock* block = m_nodeBlockFirst;
		m_nodeBlockFirst = block -> m_next;
		delete block;
	}
	while (m_edgeBlockFirst)
	{
		EdgeBlock* block = m_edgeBlockFirst;
		m_edgeBlockFirst = block -> m_next;
		delete block;
	}
	m_orphanQueue.DeAllocate();
	m_orphanStack.DeAllocate();
	m_isAllocated = false;
}
