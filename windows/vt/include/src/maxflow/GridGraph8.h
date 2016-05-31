/**************************************************************************\
*
* Copyright (c) 2004, Microsoft Corp.  All Rights Reserved.
*
* Module Name:
*
*   GridGraph8.h
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
*   Nodes are referenced using struct GridGraph8<NLinkFlowType,TotalFlowType>::NodeId. Nodes are stored in a row-by-row order.
*   Below is an example of fast pass over nodes:
*
*   GridGraph8<int,int>::NodeId p = graph->GetNodeId(0, 0);
*   for (int y=0; y<graph->GetSizeY(); y++)
*   for (int x=0; x<graph->GetSizeX(); x++)
*   {
*       ... (your code; p corresponds to node (x, y))
*       
*       p ++;
*   }
*
*   Edges are referenced using integers from 0 to 7
*   (going in clockwise direction: 0 corresponds to shift (1,0), 1 - to (1,1), and so on)
*
* Modification History:
*
*   Creation vnk 03/23/04
*
\**************************************************************************/

#pragma once

#include "MaxflowGraph.h"
using namespace GraphcutAPI;

const int gridGraph8XShifts[8] = { 1, 1, 0, -1, -1, -1,  0,  1 };
const int gridGraph8YShifts[8] = { 0, 1, 1,  1,  0, -1, -1, -1 };


template <class NLinkFlowType, class TotalFlowType> 
  struct GridGraph8<NLinkFlowType,TotalFlowType>::Square
{
    int     m_x, m_y;       // coordinates of left top corner
    NodeId  m_left_top;     // left top corner
    NodeId  m_left_bottom;  // left bottom corner == GetNeib(m_left_top, 2)
    NodeId  m_right_top;    // right top corner == GetNeib(m_left_top, 0)
    NodeId  m_right_bottom; // right bottom corner == GetNeib(m_left_top, 1)
};

template <class NLinkFlowType, class TotalFlowType> 
  struct GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2
{
    void Clear();
    void Set(TotalFlowType E00, TotalFlowType E01, TotalFlowType E10, TotalFlowType E11);

    // Given graph with the edges shown below, computes virtual graph on nodes p,q.
    //    p0 -- p
    //    | \  /
    //    |  \/   
    //    |  /\
    //    | /  \
    //    q0 -- q
    void Set(
            TotalFlowType w_p0, TotalFlowType w_q0,
            TotalFlowType w_p, TotalFlowType w_q,
            TotalFlowType w_p0_q0, TotalFlowType w_q0_p0,
            TotalFlowType w_p0_p, TotalFlowType w_p_p0,
            TotalFlowType w_p0_q, TotalFlowType w_q_p0,
            TotalFlowType w_q0_p, TotalFlowType w_p_q0,
            TotalFlowType w_q0_q, TotalFlowType w_q_q0
            );

    // construct a graph on two nodes p,q which has the same costs for minimum cuts
    // (up to a constant) as the graph defined by 'this'.
    // w1 is set to the cost of the edge source->p
    // w2 is set to the cost of the edge source->q
    // w12 is set to the cost of the edge p->q
    void ConstructGraph(TotalFlowType& w_p, TotalFlowType& w_q, TotalFlowType& w_pq);

    TotalFlowType   m_E01;
    TotalFlowType   m_E10;
    TotalFlowType   m_E11;
}; // end struct VirtualGraph2


template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::AddTEdge(NodeId p, TLinkFlowType w)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    m_nodes[p].m_tCap = (TLinkFlowType)(m_nodes[p].m_tCap + w);
}
template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::AddNEdge(NodeId p, int edge, NLinkFlowType w)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    ASSERT(GetX(p)+gridGraph8XShifts[edge]>=0 && GetX(p)+gridGraph8XShifts[edge]<m_sizeX);
    ASSERT(GetY(p)+gridGraph8YShifts[edge]>=0 && GetY(p)+gridGraph8YShifts[edge]<m_sizeY);
    m_nodes[p].m_nCap[edge] = (NLinkFlowType)(m_nodes[p].m_nCap[edge] + w);
}
template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::AddNEdges(NodeId p, int edge, NLinkFlowType w, NLinkFlowType w_rev)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    ASSERT(GetX(p)+gridGraph8XShifts[edge]>=0 && GetX(p)+gridGraph8XShifts[edge]<m_sizeX);
    ASSERT(GetY(p)+gridGraph8YShifts[edge]>=0 && GetY(p)+gridGraph8YShifts[edge]<m_sizeY);
    m_nodes[p].m_nCap[edge] = (NLinkFlowType)(m_nodes[p].m_nCap[edge] + w);
    NLinkFlowType& cap = m_nodes[GetNeib(p, edge)].m_nCap[GetReverseEdge(edge)];
    cap = (NLinkFlowType)(cap + w_rev);
}
template <class NLinkFlowType, class TotalFlowType> 
  inline typename GridGraph8<NLinkFlowType,TotalFlowType>::TLinkFlowType GridGraph8<NLinkFlowType,TotalFlowType>::GetTEdge(NodeId p)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    return m_nodes[p].m_tCap;
}
template <class NLinkFlowType, class TotalFlowType> 
  inline typename NLinkFlowType 
    GridGraph8<NLinkFlowType,TotalFlowType>::GetNEdge(NodeId p, int edge)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    ASSERT(GetX(p)+gridGraph8XShifts[edge]>=0 && GetX(p)+gridGraph8XShifts[edge]<m_sizeX);
    ASSERT(GetY(p)+gridGraph8YShifts[edge]>=0 && GetY(p)+gridGraph8YShifts[edge]<m_sizeY);
    return m_nodes[p].m_nCap[edge];
}
template <class NLinkFlowType, class TotalFlowType> 
  inline int GridGraph8<NLinkFlowType,TotalFlowType>::GetSegmentation(NodeId p)
{
    ASSERT(IsAllocated());
    ASSERT(p>=0 && p<m_nodeNum);
    return (m_nodes[p].m_parent==PARENT_FREE || m_nodes[p].m_tree) ? 1 : 0;
}



template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::SetSquare(int x, int y, Square& s)
{
    s.m_x = x;
    s.m_y = y;
    s.m_left_top = GetNodeId(x, y);
    s.m_left_bottom = GetNeib(s.m_left_top,2);
    s.m_right_top = GetNeib(s.m_left_top,0);
    s.m_right_bottom = GetNeib(s.m_left_top,1);

}

template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::MoveSquareRight(Square& s)
{
    ASSERT(GetNodeId(s.m_x, s.m_y) == s.m_left_top);
    ASSERT(s.m_left_bottom==GetNeib(s.m_left_top,2) && s.m_right_top==GetNeib(s.m_left_top,0) && s.m_right_bottom==GetNeib(s.m_left_top,1));
    s.m_x ++;
    s.m_left_top = s.m_right_top;
    s.m_right_top ++;
    s.m_left_bottom = s.m_right_bottom;
    s.m_right_bottom ++;
}

template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::MoveSquareDown(Square& s)
{
    ASSERT(GetNodeId(s.m_x, s.m_y) == s.m_left_top);
    ASSERT(s.m_left_bottom==GetNeib(s.m_left_top,2) && s.m_right_top==GetNeib(s.m_left_top,0) && s.m_right_bottom==GetNeib(s.m_left_top,1));
    s.m_y ++;
    s.m_left_top = s.m_left_bottom;
    s.m_left_bottom = GetNeib(s.m_left_bottom, 2);
    s.m_right_top = s.m_right_bottom;
    s.m_right_bottom = s.m_left_bottom + 1;
}

template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::MoveSquareLeft(Square& s)
{
    ASSERT(GetNodeId(s.m_x, s.m_y) == s.m_left_top);
    ASSERT(s.m_left_bottom==GetNeib(s.m_left_top,2) && s.m_right_top==GetNeib(s.m_left_top,0) && s.m_right_bottom==GetNeib(s.m_left_top,1));
    s.m_x --;
    s.m_right_top = s.m_left_top;
    s.m_left_top --;
    s.m_right_bottom = s.m_left_bottom;
    s.m_left_bottom --;
}

template <class NLinkFlowType, class TotalFlowType> 
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::MoveSquareUp(Square& s)
{
    ASSERT(GetNodeId(s.m_x, s.m_y) == s.m_left_top);
    ASSERT(s.m_left_bottom==GetNeib(s.m_left_top,2) && s.m_right_top==GetNeib(s.m_left_top,0) && s.m_right_bottom==GetNeib(s.m_left_top,1));
    s.m_y --;
    s.m_left_bottom = s.m_left_top;
    s.m_left_top = GetNeib(s.m_left_top, 6);
    s.m_right_bottom = s.m_right_top;
    s.m_right_top = s.m_left_top + 1;
}
