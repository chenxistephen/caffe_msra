/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include "GridGraph4Plus.h"

#include "instances.inl"

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::GridGraph4Plus(ProgressCallback abortCallback, void* callbackData)
    : m_nodesArray(NULL),
	  m_nodes(NULL),
	  m_sizeX(NULL),
	  m_sizeY(NULL),
	  m_extraEdges0(NULL),
	  m_extraArcs(NULL),
	  m_nodeShifts(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::GridGraph4Plus(int gridNum, int* sizeX, int* sizeY, int extraEdgeNumMax, 
	                                                          ProgressCallback abortCallback, void* callbackData)
    : m_nodesArray(NULL),
	  m_nodes(NULL),
	  m_sizeX(NULL),
	  m_sizeY(NULL),
	  m_extraEdges0(NULL),
	  m_extraArcs(NULL),
	  m_nodeShifts(NULL),
      m_abortCallback(abortCallback),
      m_callbackData(callbackData)
{
    Allocate(gridNum, sizeX, sizeY, extraEdgeNumMax);
}

template <class NLinkFlowType,class TotalFlowType> 
  GridGraph4Plus<NLinkFlowType,TotalFlowType>::~GridGraph4Plus()
{
    DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::Allocate(int gridNum, int* sizeX, int* sizeY, int extraEdgeNumMax)
{
	GridId t;
    HRESULT hr;

    DeAllocate();

	m_extraEdgeNum = 0; 
	m_extraEdgeNumMax = extraEdgeNumMax;
	m_extraEdges0 = new(std::nothrow) ExtraEdge0[extraEdgeNumMax];
    if (m_extraEdges0 == NULL)
    {
        return E_OUTOFMEMORY;
    }

	m_gridNum = gridNum;
	m_sizeX = new(std::nothrow) int[2*m_gridNum];
    if (m_sizeX == NULL)
    {
        return E_OUTOFMEMORY;
    }
	m_sizeY = m_sizeX + m_gridNum;


	m_allocatedNodeNum = sizeX[0]; // put dummy nodes before the first grid, so that 
	                               // traversing incident edges would not require checking whether the node is at a grid border
	for (t=0; t<m_gridNum; t++)
	{
		if (t > 0) m_allocatedNodeNum += (sizeX[t] > sizeX[t-1]) ? 
			(sizeX[t] - sizeX[t-1]) : (sizeX[t-1] - sizeX[t]); // insert dummy nodes between grids
		m_sizeX[t] = sizeX[t];
		m_sizeY[t] = sizeY[t];
		m_allocatedNodeNum += sizeX[t]*sizeY[t];
	}
	m_allocatedNodeNum += sizeX[t-1];


	m_nodes = new(std::nothrow) NodeId[m_gridNum];
    if (m_nodes == NULL)
    {
        return E_OUTOFMEMORY;
    }
    m_nodesArray = new(std::nothrow) Node[m_allocatedNodeNum];
    if (m_nodesArray == NULL)
    {
        return E_OUTOFMEMORY;
    }
    memset(m_nodesArray, 0, m_allocatedNodeNum*sizeof(Node));



	m_allocatedNodeNum = sizeX[0];
	for (t=0; t<m_gridNum; t++)
	{
		if (t > 0) m_allocatedNodeNum += (sizeX[t] > sizeX[t-1]) ? 
			(sizeX[t] - sizeX[t-1]) : (sizeX[t-1] - sizeX[t]); // insert dummy nodes between grids
		m_nodes[t] = m_allocatedNodeNum;
		for (Node* p = m_nodesArray+m_nodes[t]; p<m_nodesArray+m_nodes[t]+m_sizeX[t]*m_sizeY[t]; p++) p->m_t = t;
		m_allocatedNodeNum += sizeX[t]*sizeY[t];
	}
	m_allocatedNodeNum += sizeX[t-1];


	
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

	m_nodeShifts = new(std::nothrow) int[4*m_gridNum];
    if (m_nodeShifts == NULL)
    {
        return E_OUTOFMEMORY;
    }
	for (t=0; t<m_gridNum; t++)
	{
		m_nodeShifts[4*t+0] = 1;
		m_nodeShifts[4*t+1] = m_sizeX[t];
		m_nodeShifts[4*t+2] = -1;
		m_nodeShifts[4*t+3] = -m_sizeX[t];
	}

    ActiveQueueReset();

    return S_OK;
}

template <class NLinkFlowType,class TotalFlowType> 
  void GridGraph4Plus<NLinkFlowType,TotalFlowType>::DeAllocate()
{
    if (m_nodes)
    {
        delete [] m_nodes;
        m_nodes = NULL;
    }
    if (m_extraEdges0)
    {
        delete [] m_extraEdges0;
        m_extraEdges0 = NULL;
    }
    if (m_extraArcs)
    {
        delete [] m_extraArcs;
        m_extraArcs = NULL;
    }
    if (m_nodesArray)
    {
        delete [] m_nodesArray;
        m_nodesArray = NULL;
    }
    if (m_sizeX)
    {
        delete [] m_sizeX;
        m_sizeX = m_sizeY = NULL;
    }
    if (m_nodeShifts)
    {
        delete [] m_nodeShifts;
        m_nodeShifts = NULL;
    }
    m_orphanQueue.DeAllocate();
    m_orphanStack.DeAllocate();
}

template <class NLinkFlowType,class TotalFlowType> 
  HRESULT GridGraph4Plus<NLinkFlowType,TotalFlowType>::InitExtraArcs()
{
	Node* i;
	ExtraEdge0* e;

	// first and last nodes are dummy nodes, so they must have zero degree
	ASSERT(m_nodesArray[0].m_extraArcs == 0 && m_nodesArray[m_allocatedNodeNum-1].m_extraArcs == 0);

	int _e = 0;
	for (i=m_nodesArray; i<m_nodesArray+m_allocatedNodeNum; i++)
	{
		int degree = i->m_extraArcs;
		i->m_extraArcs = _e;
		_e += degree;
	}

	m_extraArcs = new(std::nothrow) ExtraArc[2*m_extraEdgeNum];
    if (m_extraArcs == NULL)
    {
        return E_OUTOFMEMORY;
    }

	for (e=m_extraEdges0; e<m_extraEdges0+m_extraEdgeNum; e++)
	{
		ExtraArcId _a = m_nodesArray[e->m_i[0]].m_extraArcs ++;
		ExtraArcId _a_rev = m_nodesArray[e->m_i[1]].m_extraArcs ++;
		ExtraArc* a     = &m_extraArcs[_a];
		ExtraArc* a_rev = &m_extraArcs[_a_rev];
		a->m_head = e->m_i[1];
		a->m_nCap = e->m_nCap[0];
		a->m_sister = _a_rev;
		a_rev->m_head = e->m_i[0];
		a_rev->m_nCap = e->m_nCap[1];
		a_rev->m_sister = _a;
	}

	_e = 0;
	for (i=m_nodesArray; i<m_nodesArray+m_allocatedNodeNum-1; i++)
	{
		int _e_prev = i->m_extraArcs;
		i->m_extraArcs = _e;
		_e = _e_prev;
	}

	return S_OK;
}
