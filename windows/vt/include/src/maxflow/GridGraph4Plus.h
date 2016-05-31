/**************************************************************************\
*
* Copyright (c) 2004, Microsoft Corp.  All Rights Reserved.
*
* Module Name:
*
*   GridGraph4Plus.h
*
* Abstract:
*
*   Maxflow algorithm for regular grid graphs with 8-neighborhood system
*   Algorithm is based on the paper
*     "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision",
*     Y. Boykov and V. Kolmogorov, PAMI'04,
*   and my PhD thesis. Additionally, speed-up techniques based on dynamic programming are implemented.
*
*
*
*   Nodes are referenced using struct GridGraph4Plus<NLinkFlowType,TotalFlowType>::NodeId. Nodes are stored in a row-by-row order.
*   Below is an example of fast pass over nodes:
*
*   GridGraph4Plus<int,int>::NodeId p = graph->GetNodeId(0, 0);
*   for (int y=0; y<graph->GetSizeY(); y++)
*   for (int x=0; x<graph->GetSizeX(); x++)
*   {
*       ... (your code; p corresponds to node (x, y))
*       
*       p ++;
*   }
*
*   Edges are referenced using integers from 0 to 4
*   (going in clockwise direction: 0 corresponds to shift (1,0), 1 - to (0,1), and so on)
*
* Modification History:
*
*   Creation vnk 03/24/05
*
\**************************************************************************/

#pragma once

#include "MaxflowGraph.h"
using namespace GraphcutAPI;

const int gridGraph4PlusXShifts[4] = { 1, 0, -1,  0};
const int gridGraph4PlusYShifts[4] = { 0, 1,  0, -1};



template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::AddTEdge(NodeId p, TLinkFlowType w)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    m_nodesArray[p].m_tCap = (TLinkFlowType)(m_nodesArray[p].m_tCap + w);
}
template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::AddNEdge(NodeId p, int edge, NLinkFlowType w)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    m_nodesArray[p].m_nCap[edge] = (NLinkFlowType)(m_nodesArray[p].m_nCap[edge] + w);
}
template <class NLinkFlowType,class TotalFlowType> 
  inline typename GridGraph4Plus<NLinkFlowType,TotalFlowType>::TLinkFlowType GridGraph4Plus<NLinkFlowType,TotalFlowType>::GetTEdge(NodeId p)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    return m_nodesArray[p].m_tCap;
}
template <class NLinkFlowType,class TotalFlowType> 
  inline typename NLinkFlowType GridGraph4Plus<NLinkFlowType,TotalFlowType>::GetNEdge(NodeId p, int edge)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    return m_nodesArray[p].m_nCap[edge];
}
template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::AddNEdges(NodeId p, int edge, NLinkFlowType w, NLinkFlowType w_rev)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    m_nodesArray[p].m_nCap[edge] = (NLinkFlowType)(m_nodesArray[p].m_nCap[edge] + w);
    NLinkFlowType& cap = m_nodesArray[GetNeib(p, edge)].m_nCap[GetReverseEdge(edge)];
    cap = (NLinkFlowType)(cap + w_rev);
}
template <class NLinkFlowType,class TotalFlowType> 
  inline void GridGraph4Plus<NLinkFlowType,TotalFlowType>::AddExtraEdges(NodeId p, NodeId q, NLinkFlowType w, NLinkFlowType w_rev)
{
    ASSERT(IsAllocated());
    ASSERT(m_extraEdgeNum < m_extraEdgeNumMax);
	ASSERT(m_nodesArray[p].m_extraArcs < 256-7);
	ASSERT(m_nodesArray[q].m_extraArcs < 256-7);
	m_extraEdges0[m_extraEdgeNum].m_i[0] = p;
	m_extraEdges0[m_extraEdgeNum].m_i[1] = q;
	m_extraEdges0[m_extraEdgeNum].m_nCap[0] = w;
	m_extraEdges0[m_extraEdgeNum].m_nCap[1] = w_rev;
	m_nodesArray[p].m_extraArcs ++;
	m_nodesArray[q].m_extraArcs ++;
	m_extraEdgeNum ++;
}
template <class NLinkFlowType,class TotalFlowType> 
  inline int GridGraph4Plus<NLinkFlowType,TotalFlowType>::GetSegmentation(NodeId p)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_allocatedNodeNum);
    return (m_nodesArray[p].m_parent==PARENT_FREE || m_nodesArray[p].m_tree) ? 1 : 0;
}




