/* Copyright 2004 Microsoft Corporation */

#include "stdafx.h"
#include "GridGraph8.h"

#include "instances.inl"

template <class T>
  inline typename T MIN(typename T a, typename T b)
{
	return (a<b) ? a : b;
}

///////////////////////////////////////////////////////////////////////////////////////////
// VirtualGraph2 functions
///////////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2::Clear()
{
	m_E01 = 0;
	m_E10 = 0;
	m_E11 = 0;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2::Set(TotalFlowType E00, TotalFlowType E01, TotalFlowType E10, TotalFlowType E11)
{
	ASSERT(E00 + E11 <= E01 + E10); // graph must be regular
	m_E01 = E01 - E00;
	m_E10 = E10 - E00;
	m_E11 = E11 - E00;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2::ConstructGraph(TotalFlowType& w_p, TotalFlowType& w_q, TotalFlowType& w_pq)
{
	ASSERT(m_E11 <= m_E01 + m_E10); // graph must be regular
	w_p = m_E10;
	w_q = m_E11 - w_p;
	w_pq = m_E01 - w_q;
}

/*
template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2::Set(
		TotalFlowType w_p0, TotalFlowType w_q0,
		TotalFlowType w_p, TotalFlowType w_q,
		TotalFlowType w_p0_q0, TotalFlowType w_q0_p0,
		TotalFlowType w_p0_p, TotalFlowType w_p_p0,
		TotalFlowType w_p0_q, TotalFlowType w_q_p0,
		TotalFlowType w_q0_p, TotalFlowType w_p_q0,
		TotalFlowType w_q0_q, TotalFlowType w_q_q0
		)
{
	ASSERT(w_p0_q0 >= 0 && w_q0_p0 >= 0);
	ASSERT(w_p0_p >= 0 && w_p_p0 >= 0);
	ASSERT(w_p0_q >= 0 && w_q_p0 >= 0);
	ASSERT(w_q0_p >= 0 && w_p_q0 >= 0);
	ASSERT(w_q0_q >= 0 && w_q_q0 >= 0);

	// push flow source->p->p0->sink
	w_p0 -= w_p0_p;
	w_p_p0 += w_p0_p;
	w_p += w_p0_p;
	// w_p0_p = 0;

	// push flow source->p->q0->sink
	w_q0 -= w_q0_p;
	w_p_q0 += w_q0_p;
	w_p += w_q0_p;
	// w_q0_p = 0;

	// push flow source->q->p0->sink
	w_p0 -= w_p0_q;
	w_q_p0 += w_p0_q;
	w_q += w_p0_q;
	// w_p0_q = 0;

	// push flow source->q->q0->sink
	w_q0 -= w_q0_q;
	w_q_q0 += w_q0_q;
	w_q += w_q0_q;
	// w_q0_q = 0;

	// push flow source->p0->q0->sink
	w_q0 -= w_q0_p0;
	w_p0_q0 += w_q0_p0;
	w_p0 += w_q0_p0;
	// w_q0_p0 = 0;

	TotalFlowType E00 = MIN<TotalFlowType>(
		MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_p_q0 + w_q_q0),
		MIN<TotalFlowType>(0, w_q0 + w_p_q0 + w_q_q0          ) + w_p0 + w_p_p0 + w_q_p0 );

	TotalFlowType E01 = MIN<TotalFlowType>(
		MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_p_q0),
		MIN<TotalFlowType>(0, w_q0 + w_p_q0          ) + w_p0 + w_p_p0);

	TotalFlowType E10 = MIN<TotalFlowType>(
		MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_q_q0),
		MIN<TotalFlowType>(0, w_q0 + w_q_q0          ) + w_p0 + w_q_p0 );

	TotalFlowType E11 = MIN<TotalFlowType>(
		MIN<TotalFlowType>(0, w_q0 + w_p0_q0),
		MIN<TotalFlowType>(0, w_q0          ) + w_p0 );

	Set(E00, E01, E10, E11);

	m_E01 += w_q;
	m_E10 += w_p;
	m_E11 += w_p + w_q;
}
*/
template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::VirtualGraph2::Set(
		TotalFlowType w_p0, TotalFlowType w_q0,
		TotalFlowType w_p, TotalFlowType w_q,
		TotalFlowType w_p0_q0, TotalFlowType w_q0_p0,
		TotalFlowType w_p0_p, TotalFlowType w_p_p0,
		TotalFlowType w_p0_q, TotalFlowType w_q_p0,
		TotalFlowType w_q0_p, TotalFlowType w_p_q0,
		TotalFlowType w_q0_q, TotalFlowType w_q_q0
		)
{
	ASSERT(w_p0_q0 >= 0 && w_q0_p0 >= 0);
	ASSERT(w_p0_p >= 0 && w_p_p0 >= 0);
	ASSERT(w_p0_q >= 0 && w_q_p0 >= 0);
	ASSERT(w_q0_p >= 0 && w_p_q0 >= 0);
	ASSERT(w_q0_q >= 0 && w_q_q0 >= 0);

	// push flow source->p->p0->sink
	w_p0 -= w_p0_p;
	w_p_p0 += w_p0_p;
	w_p += w_p0_p;
	// w_p0_p = 0;

	// push flow source->p->q0->sink
	w_q0 -= w_q0_p;
	w_p_q0 += w_q0_p;
	w_p += w_q0_p;
	// w_q0_p = 0;

	// push flow source->q->p0->sink
	w_p0 -= w_p0_q;
	w_q_p0 += w_p0_q;
	w_q += w_p0_q;
	// w_p0_q = 0;

	// push flow source->q->q0->sink
	w_q0 -= w_q0_q;
	w_q_q0 += w_q0_q;
	w_q += w_q0_q;
	// w_q0_q = 0;

	// E(xp, xq) = maximum flow that can be pushed through the graph
	// when 'p' is connected to terminal 'xp' and 'q' is connected to terminal 'xq'
	// (xp,xq are 0 (source) or 1 (sink)

	if (w_p0 >= w_p0_q0)
	{
		w_q0 += w_p0_q0;

		// now only  edges w_q0, w_p_q0 and w_q_q0 matter - flow trough other edges will not go
		if (w_q0 >= 0)
		{
			m_E01 = w_q;
			m_E10 = w_p;
			m_E11 = w_p + w_q;
			return;
		}

		TotalFlowType E10 = MIN<TotalFlowType>(-w_q0, w_q_q0);
		if (-w_q0 <= w_p_q0)
		{
			m_E01 = 0;
			m_E10 = E10 + w_q0;
			m_E11 = w_q0;
		}
		else
		{
			// E00 = MIN<TotalFlowType>(-w_q0, w_p_q0 + w_q_q0), E01 = w_p_q0, E10 = E10, E11 = 0
			Set(MIN<TotalFlowType>(-w_q0, w_p_q0 + w_q_q0), w_p_q0, E10, 0);
		}
	}
	else
	{
		// w_p0 < w_p0_q0
		if (-w_p0 <= w_q0_p0)
		{
			w_q0 += w_p0;
			if (w_q0 >= 0)
			{
				m_E01 = w_q;
				m_E10 = w_p;
				m_E11 = w_p + w_q;
				return;
			}
			w_p0_q0 -= w_p0;

			// w_p0 is now zero, w_q0 < 0
			// only edges w_q0, w_p_q0, w_q_q0, w_p_p0, w_q_p0, w_p0_q0 matter - flow trough other edges will not go

			TotalFlowType flow01 = w_p_q0 + MIN<TotalFlowType>(w_p0_q0, w_p_p0);
			TotalFlowType flow10 = w_q_q0 + MIN<TotalFlowType>(w_p0_q0, w_q_p0);

			if (-w_q0 <= flow01)
			{
				// E00 = -w_q0, E01 = -w_q0, E10 = MIN<TotalFlowType>(-w_q0, flow10), E11 = 0
				m_E01 =                             w_q;
				m_E10 = w_q0 + MIN<TotalFlowType>(-w_q0, flow10) + w_p;
				m_E11 = w_q0                      + w_p + w_q;
				return;
			}

			if (-w_q0 <= flow10)
			{
				// E00 = -w_q0, E01 = flow01, E10 = -w_q0, E11 = 0
				m_E01 = w_q0 + flow01 + w_q;
				m_E10 =                 w_p;
				m_E11 = w_q0          + w_p + w_q;
				return;
			}

			TotalFlowType flow11 = w_p_q0 + w_q_q0 + MIN<TotalFlowType>(w_p0_q0, w_p_p0 + w_q_p0);

			// E00 = MIN<TotalFlowType>(-w_pq0, flow11), E01 = flow01, E10 = flow10, E11 = 0
			Set(MIN<TotalFlowType>(-w_q0, flow11), flow01, flow10, 0);
		}
		else
		{
			// w_p0 < -w_q0_p0
			w_p0 += w_q0_p0;
			w_p0_q0 += w_q0_p0;
			w_q0 -= w_q0_p0;

			TotalFlowType E00 = MIN<TotalFlowType>(
				MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_p_q0 + w_q_q0),
				MIN<TotalFlowType>(0, w_q0 + w_p_q0 + w_q_q0          ) + w_p0 + w_p_p0 + w_q_p0 );

			TotalFlowType E01 = MIN<TotalFlowType>(
				MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_p_q0),
				MIN<TotalFlowType>(0, w_q0 + w_p_q0          ) + w_p0 + w_p_p0);

			TotalFlowType E10 = MIN<TotalFlowType>(
				MIN<TotalFlowType>(0, w_q0 + w_p0_q0 + w_q_q0),
				MIN<TotalFlowType>(0, w_q0 + w_q_q0          ) + w_p0 + w_q_p0 );

//			TotalFlowType E11 = MIN<TotalFlowType>(
//				MIN<TotalFlowType>(0, w_q0 + w_p0_q0),
//				MIN<TotalFlowType>(0, w_q0          ) + w_p0 );

			TotalFlowType E11 = MIN<TotalFlowType>(0, w_q0) + w_p0;

			Set(E00, E01, E10, E11);
		}

	}

	m_E01 += w_q;
	m_E10 += w_p;
	m_E11 += w_p + w_q;
}

///////////////////////////////////////////////////////////////////////////////////////////
// Functions ComputeVirtualGraphFrom<Side1...SideN>To<Side>()
///////////////////////////////////////////////////////////////////////////////////////////


template <class NLinkFlowType, class TotalFlowType>
  void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromRightToLeft(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG)
{
	// p0=right_top, q0=right_bottom, p=left_top, q=left_bottom
	
	TotalFlowType rightVG_p0, rightVG_q0, rightVG_p0_q0;
	rightVG->ConstructGraph(rightVG_p0, rightVG_q0, rightVG_p0_q0);

	leftVG->Set(
		 m_nodes[s.m_right_top].m_tCap + rightVG_p0,       m_nodes[s.m_right_bottom].m_tCap + rightVG_q0,
		 0,                                        0,
		 m_nodes[s.m_right_top].m_nCap[2] + rightVG_p0_q0, m_nodes[s.m_right_bottom].m_nCap[6],
		 m_nodes[s.m_right_top].m_nCap[4],                 m_nodes[s.m_left_top].m_nCap[0],
		 m_nodes[s.m_right_top].m_nCap[3],                 m_nodes[s.m_left_bottom].m_nCap[7],
		 m_nodes[s.m_right_bottom].m_nCap[5],              m_nodes[s.m_left_top].m_nCap[1],
		 m_nodes[s.m_right_bottom].m_nCap[4],              m_nodes[s.m_left_bottom].m_nCap[0]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromTopToBottom(Square& s, VirtualGraph2* topVG, VirtualGraph2* bottomVG)
{
	// p0=left_top, q0=right_top, p=left_bottom, q=right_bottom
	
	TotalFlowType topVG_p0, topVG_q0, topVG_p0_q0;
	topVG->ConstructGraph(topVG_p0, topVG_q0, topVG_p0_q0);

	bottomVG->Set(
		 m_nodes[s.m_left_top].m_tCap + topVG_p0,       m_nodes[s.m_right_top].m_tCap + topVG_q0,
		 0,                                     0,
		 m_nodes[s.m_left_top].m_nCap[0] + topVG_p0_q0, m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[2],               m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[1],               m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_right_top].m_nCap[3],              m_nodes[s.m_left_bottom].m_nCap[7],
		 m_nodes[s.m_right_top].m_nCap[2],              m_nodes[s.m_right_bottom].m_nCap[6]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromRightLeftToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* bottomVG)
{
	// p0=left_top, q0=right_top, p=left_bottom, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_p, leftVG_p0_p;
	leftVG->ConstructGraph(leftVG_p0, leftVG_p, leftVG_p0_p);

	TotalFlowType rightVG_q0, rightVG_q, rightVG_q0_q;
	rightVG->ConstructGraph(rightVG_q0, rightVG_q, rightVG_q0_q);

	bottomVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0,        m_nodes[s.m_right_top].m_tCap + rightVG_q0,
		                        leftVG_p,                                 rightVG_q,
		 m_nodes[s.m_left_top].m_nCap[0],                 m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_p,   m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[1],                 m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_right_top].m_nCap[3],                m_nodes[s.m_left_bottom].m_nCap[7],
		 m_nodes[s.m_right_top].m_nCap[2] + rightVG_q0_q, m_nodes[s.m_right_bottom].m_nCap[6]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromLeftTopToBottom(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG)
{
	// p0=left_top, q0=right_top, p=left_bottom, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_p, leftVG_p0_p;
	leftVG->ConstructGraph(leftVG_p0, leftVG_p, leftVG_p0_p);

	TotalFlowType topVG_p0, topVG_q0, topVG_p0_q0;
	topVG->ConstructGraph(topVG_p0, topVG_q0, topVG_p0_q0);

	bottomVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0 + topVG_p0, m_nodes[s.m_right_top].m_tCap + topVG_q0,
		                        leftVG_p,             0,
		 m_nodes[s.m_left_top].m_nCap[0] + topVG_p0_q0,       m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_p,       m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[1],                     m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_right_top].m_nCap[3],                    m_nodes[s.m_left_bottom].m_nCap[7],
		 m_nodes[s.m_right_top].m_nCap[2],                    m_nodes[s.m_right_bottom].m_nCap[6]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromRightLeftTopToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG)
{
	// p0=left_top, q0=right_top, p=left_bottom, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_p, leftVG_p0_p;
	leftVG->ConstructGraph(leftVG_p0, leftVG_p, leftVG_p0_p);

	TotalFlowType topVG_p0, topVG_q0, topVG_p0_q0;
	topVG->ConstructGraph(topVG_p0, topVG_q0, topVG_p0_q0);

	TotalFlowType rightVG_q0, rightVG_q, rightVG_q0_q;
	rightVG->ConstructGraph(rightVG_q0, rightVG_q, rightVG_q0_q);

	bottomVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0 + topVG_p0, m_nodes[s.m_right_top].m_tCap + rightVG_q0 + topVG_q0,
		                        leftVG_p,             rightVG_q,
		 m_nodes[s.m_left_top].m_nCap[0] + topVG_p0_q0,       m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_p,       m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[1],                     m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_right_top].m_nCap[3],                    m_nodes[s.m_left_bottom].m_nCap[7],
		 m_nodes[s.m_right_top].m_nCap[2] + rightVG_q0_q,     m_nodes[s.m_right_bottom].m_nCap[6]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromLeftToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* rightVG)
{
	// p0=left_top, q0=left_bottom, p=right_top, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_q0, leftVG_p0_q0;
	leftVG->ConstructGraph(leftVG_p0, leftVG_q0, leftVG_p0_q0);

	rightVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0,       m_nodes[s.m_left_bottom].m_tCap + leftVG_q0,
		 0,                                      0,
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_q0, m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[0],                m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[1],                m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_left_bottom].m_nCap[7],             m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_left_bottom].m_nCap[0],             m_nodes[s.m_right_bottom].m_nCap[4]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromLeftBottomToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* bottomVG, VirtualGraph2* rightVG)
{
	// p0=left_top, q0=left_bottom, p=right_top, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_q0, leftVG_p0_q0;
	leftVG->ConstructGraph(leftVG_p0, leftVG_q0, leftVG_p0_q0);

	TotalFlowType bottomVG_q0, bottomVG_q, bottomVG_q0_q;
	bottomVG->ConstructGraph(bottomVG_q0, bottomVG_q, bottomVG_q0_q);

	rightVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0,           m_nodes[s.m_left_bottom].m_tCap + leftVG_q0 + bottomVG_q0,
		 0,                                          bottomVG_q,
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_q0,     m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[0],                    m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[1],                    m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_left_bottom].m_nCap[7],                 m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_left_bottom].m_nCap[0] + bottomVG_q0_q, m_nodes[s.m_right_bottom].m_nCap[4]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromLeftTopToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* rightVG)
{
	// p0=left_top, q0=left_bottom, p=right_top, q=right_bottom
	
	TotalFlowType leftVG_p0, leftVG_q0, leftVG_p0_q0;
	leftVG->ConstructGraph(leftVG_p0, leftVG_q0, leftVG_p0_q0);

	TotalFlowType topVG_p0, topVG_p, topVG_p0_p;
	topVG->ConstructGraph(topVG_p0, topVG_p, topVG_p0_p);

	rightVG->Set(
		 m_nodes[s.m_left_top].m_tCap + leftVG_p0 + topVG_p0, m_nodes[s.m_left_bottom].m_tCap + leftVG_q0,
		 topVG_p,                                     0,
		 m_nodes[s.m_left_top].m_nCap[2] + leftVG_p0_q0,      m_nodes[s.m_left_bottom].m_nCap[6],
		 m_nodes[s.m_left_top].m_nCap[0] + topVG_p0_p,        m_nodes[s.m_right_top].m_nCap[4],
		 m_nodes[s.m_left_top].m_nCap[1],                     m_nodes[s.m_right_bottom].m_nCap[5],
		 m_nodes[s.m_left_bottom].m_nCap[7],                  m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_left_bottom].m_nCap[0],                  m_nodes[s.m_right_bottom].m_nCap[4]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromBottomToTop(Square& s, VirtualGraph2* bottomVG, VirtualGraph2* topVG)
{
	// p0=left_bottom, q0=right_bottom, p=left_top, q=right_top
	
	TotalFlowType bottomVG_p0, bottomVG_q0, bottomVG_p0_q0;
	bottomVG->ConstructGraph(bottomVG_p0, bottomVG_q0, bottomVG_p0_q0);

	topVG->Set(
		 m_nodes[s.m_left_bottom].m_tCap + bottomVG_p0,       m_nodes[s.m_right_bottom].m_tCap + bottomVG_q0,
		 0,                                           0,
		 m_nodes[s.m_left_bottom].m_nCap[0] + bottomVG_p0_q0, m_nodes[s.m_right_bottom].m_nCap[4],
		 m_nodes[s.m_left_bottom].m_nCap[6],                  m_nodes[s.m_left_top].m_nCap[2],
		 m_nodes[s.m_left_bottom].m_nCap[7],                  m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_right_bottom].m_nCap[5],                 m_nodes[s.m_left_top].m_nCap[1],
		 m_nodes[s.m_right_bottom].m_nCap[6],                 m_nodes[s.m_right_top].m_nCap[2]
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromRightLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG)
{
	// p0=left_bottom, q0=right_bottom, p=left_top, q=right_top
	
	TotalFlowType leftVG_p, leftVG_p0, leftVG_p_p0;
	leftVG->ConstructGraph(leftVG_p, leftVG_p0, leftVG_p_p0);

	TotalFlowType rightVG_q, rightVG_q0, rightVG_q_q0;
	rightVG->ConstructGraph(rightVG_q, rightVG_q0, rightVG_q_q0);

	topVG->Set(
		 m_nodes[s.m_left_bottom].m_tCap + leftVG_p0, m_nodes[s.m_right_bottom].m_tCap + rightVG_q0,
		 leftVG_p,                            rightVG_q,
		 m_nodes[s.m_left_bottom].m_nCap[0],          m_nodes[s.m_right_bottom].m_nCap[4],
		 m_nodes[s.m_left_bottom].m_nCap[6],          m_nodes[s.m_left_top].m_nCap[2] + leftVG_p_p0,
		 m_nodes[s.m_left_bottom].m_nCap[7],          m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_right_bottom].m_nCap[5],         m_nodes[s.m_left_top].m_nCap[1],
		 m_nodes[s.m_right_bottom].m_nCap[6],         m_nodes[s.m_right_top].m_nCap[2] + rightVG_q_q0
		);
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeVirtualGraphFromRightBottomLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* bottomVG, VirtualGraph2* leftVG, VirtualGraph2* topVG)
{
	// p0=left_bottom, q0=right_bottom, p=left_top, q=right_top
	
	TotalFlowType bottomVG_p0, bottomVG_q0, bottomVG_p0_q0;
	bottomVG->ConstructGraph(bottomVG_p0, bottomVG_q0, bottomVG_p0_q0);

	TotalFlowType leftVG_p, leftVG_p0, leftVG_p_p0;
	leftVG->ConstructGraph(leftVG_p, leftVG_p0, leftVG_p_p0);

	TotalFlowType rightVG_q, rightVG_q0, rightVG_q_q0;
	rightVG->ConstructGraph(rightVG_q, rightVG_q0, rightVG_q_q0);

	topVG->Set(
		 m_nodes[s.m_left_bottom].m_tCap + bottomVG_p0 + leftVG_p0, m_nodes[s.m_right_bottom].m_tCap + bottomVG_q0 + rightVG_q0,
		 leftVG_p,                                          rightVG_q,
		 m_nodes[s.m_left_bottom].m_nCap[0] + bottomVG_p0_q0,       m_nodes[s.m_right_bottom].m_nCap[4],
		 m_nodes[s.m_left_bottom].m_nCap[6],                        m_nodes[s.m_left_top].m_nCap[2] + leftVG_p_p0,
		 m_nodes[s.m_left_bottom].m_nCap[7],                        m_nodes[s.m_right_top].m_nCap[3],
		 m_nodes[s.m_right_bottom].m_nCap[5],                       m_nodes[s.m_left_top].m_nCap[1],
		 m_nodes[s.m_right_bottom].m_nCap[6],                       m_nodes[s.m_right_top].m_nCap[2] + rightVG_q_q0
		);
}

///////////////////////////////////////////////////////////////////////////////////////////
// Functions ComputeFlow<Side>()
///////////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowLeft(Square& s, VirtualGraph2 *leftVG, TotalFlowType* left_top_flow, TotalFlowType *left_bottom_flow)
{
	TotalFlowType w_p, w_q, w_pq, flow_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, t_old, delta;

	leftVG->ConstructGraph(w_p, w_q, w_pq);
	flow_pq = 0; // flow through edge with capacity w_pq

	t_p = m_nodes[s.m_left_top].m_tCap + w_p;
	t_q = m_nodes[s.m_left_bottom].m_tCap + w_q;

	// Push flow through through n-links within square s
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	// Calculate the amount of flow that must be pushed through w0_pq (but do not push it).
	if (t_p>0 && t_q<0)
	{
		t_old = t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[2]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[2] -= delta; m_nodes[s.m_left_bottom].m_nCap[6] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[0]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[3]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[0] -= delta; m_nodes[s.m_right_top].m_nCap[4] += delta;
			m_nodes[s.m_right_top].m_nCap[3] -= delta; m_nodes[s.m_left_bottom].m_nCap[7] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[1]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[4]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[1] -= delta; m_nodes[s.m_right_bottom].m_nCap[5] += delta;
			m_nodes[s.m_right_bottom].m_nCap[4] -= delta; m_nodes[s.m_left_bottom].m_nCap[0] += delta;
		}

		flow_pq = MIN<TotalFlowType>(t, w_pq);

		delta = t_old - t;
		t_p -= delta;
		t_q += delta;
	}
	else if (t_p<0 && t_q>0)
	{
		t_old = t = MIN<TotalFlowType>(-t_p, t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[6]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[2] += delta; m_nodes[s.m_left_bottom].m_nCap[6] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[4]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_bottom].m_nCap[7]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[0] += delta; m_nodes[s.m_right_top].m_nCap[4] -= delta;
			m_nodes[s.m_right_top].m_nCap[3] += delta; m_nodes[s.m_left_bottom].m_nCap[7] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_bottom].m_nCap[5]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_bottom].m_nCap[0]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[1] += delta; m_nodes[s.m_right_bottom].m_nCap[5] -= delta;
			m_nodes[s.m_right_bottom].m_nCap[4] += delta; m_nodes[s.m_left_bottom].m_nCap[0] -= delta;
		}
		delta = t_old - t;
		t_p += delta;
		t_q -= delta;
	}

	// now separate flow that goes through leftVG and flow that goes trough real t-links of p,q
	t_p -= w_p;
	t_q -= w_q;

	w_p -= flow_pq;
	w_q += flow_pq;

	if (w_p<0 && t_p>0)
	{
		*left_top_flow = MIN<TotalFlowType>(-w_p, t_p);
	}
	else if (w_p>0 && t_p<0)
	{
		*left_top_flow = -MIN<TotalFlowType>(w_p, -t_p);
	}
	else
	{
		*left_top_flow = 0;
	}
	m_nodes[s.m_left_top].m_tCap = t_p - *left_top_flow;

	if (w_q<0 && t_q>0)
	{
		*left_bottom_flow = MIN<TotalFlowType>(-w_q, t_q);
	}
	else if (w_q>0 && t_q<0)
	{
		*left_bottom_flow = -MIN<TotalFlowType>(w_q, -t_q);
	}
	else
	{
		*left_bottom_flow = 0;
	}
	m_nodes[s.m_left_bottom].m_tCap = t_q - *left_bottom_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowRight(Square& s, VirtualGraph2 *rightVG, TotalFlowType* right_top_flow, TotalFlowType *right_bottom_flow)
{
	TotalFlowType w_p, w_q, w_pq, flow_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, t_old, delta;

	rightVG->ConstructGraph(w_p, w_q, w_pq);
	flow_pq = 0; // flow through edge with capacity w_pq

	t_p = m_nodes[s.m_right_top].m_tCap + w_p;
	t_q = m_nodes[s.m_right_bottom].m_tCap + w_q;

	// Push flow through through n-links within square s
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	// Calculate the amount of flow that must be pushed through w0_pq (but do not push it).
	if (t_p>0 && t_q<0)
	{
		t_old = t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[2]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[2] -= delta; m_nodes[s.m_right_bottom].m_nCap[6] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[4]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_top].m_nCap[1]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[4] -= delta; m_nodes[s.m_left_top].m_nCap[0] += delta;
			m_nodes[s.m_left_top].m_nCap[1] -= delta; m_nodes[s.m_right_bottom].m_nCap[5] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[3]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_bottom].m_nCap[0]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[3] -= delta; m_nodes[s.m_left_bottom].m_nCap[7] += delta;
			m_nodes[s.m_left_bottom].m_nCap[0] -= delta; m_nodes[s.m_right_bottom].m_nCap[4] += delta;
		}

		flow_pq = MIN<TotalFlowType>(t, w_pq);

		delta = t_old - t;
		t_p -= delta;
		t_q += delta;
	}
	else if (t_p<0 && t_q>0)
	{
		t_old = t = MIN<TotalFlowType>(-t_p, t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_bottom].m_nCap[6]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[2] += delta; m_nodes[s.m_right_bottom].m_nCap[6] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[0]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[5]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[4] += delta; m_nodes[s.m_left_top].m_nCap[0] -= delta;
			m_nodes[s.m_left_top].m_nCap[1] += delta; m_nodes[s.m_right_bottom].m_nCap[5] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[7]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[4]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_right_top].m_nCap[3] += delta; m_nodes[s.m_left_bottom].m_nCap[7] -= delta;
			m_nodes[s.m_left_bottom].m_nCap[0] += delta; m_nodes[s.m_right_bottom].m_nCap[4] -= delta;
		}
		delta = t_old - t;
		t_p += delta;
		t_q -= delta;
	}

	// now separate flow that goes through rightVG and flow that goes trough real t-links of p,q
	t_p -= w_p;
	t_q -= w_q;

	w_p -= flow_pq;
	w_q += flow_pq;

	if (w_p<0 && t_p>0)
	{
		*right_top_flow = MIN<TotalFlowType>(-w_p, t_p);
	}
	else if (w_p>0 && t_p<0)
	{
		*right_top_flow = -MIN<TotalFlowType>(w_p, -t_p);
	}
	else
	{
		*right_top_flow = 0;
	}
	m_nodes[s.m_right_top].m_tCap = t_p - *right_top_flow;

	if (w_q<0 && t_q>0)
	{
		*right_bottom_flow = MIN<TotalFlowType>(-w_q, t_q);
	}
	else if (w_q>0 && t_q<0)
	{
		*right_bottom_flow = -MIN<TotalFlowType>(w_q, -t_q);
	}
	else
	{
		*right_bottom_flow = 0;
	}
	m_nodes[s.m_right_bottom].m_tCap = t_q - *right_bottom_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowTop(Square& s, VirtualGraph2 *topVG, TotalFlowType *left_top_flow, TotalFlowType *right_top_flow)
{
	TotalFlowType w_p, w_q, w_pq, flow_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, t_old, delta;

	topVG->ConstructGraph(w_p, w_q, w_pq);
	flow_pq = 0; // flow through edge with capacity w_pq

	t_p = m_nodes[s.m_left_top].m_tCap + w_p;
	t_q = m_nodes[s.m_right_top].m_tCap + w_q;

	// Push flow through through n-links within square s
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	// Calculate the amount of flow that must be pushed through w0_pq (but do not push it).
	if (t_p>0 && t_q<0)
	{
		t_old = t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[0]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[0] -= delta; m_nodes[s.m_right_top].m_nCap[4] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[2]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_bottom].m_nCap[7]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[2] -= delta; m_nodes[s.m_left_bottom].m_nCap[6] += delta;
			m_nodes[s.m_left_bottom].m_nCap[7] -= delta; m_nodes[s.m_right_top].m_nCap[3] += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[1]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[6]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[1] -= delta; m_nodes[s.m_right_bottom].m_nCap[5] += delta;
			m_nodes[s.m_right_bottom].m_nCap[6] -= delta; m_nodes[s.m_right_top].m_nCap[2] += delta;
		}

		flow_pq = MIN<TotalFlowType>(t, w_pq);

		delta = t_old - t;
		t_p -= delta;
		t_q += delta;
	}
	else if (t_p<0 && t_q>0)
	{
		t_old = t = MIN<TotalFlowType>(-t_p, t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[4]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[0] += delta; m_nodes[s.m_right_top].m_nCap[4] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[6]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[3]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[2] += delta; m_nodes[s.m_left_bottom].m_nCap[6] -= delta;
			m_nodes[s.m_left_bottom].m_nCap[7] += delta; m_nodes[s.m_right_top].m_nCap[3] -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_bottom].m_nCap[5]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[2]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_nCap[1] += delta; m_nodes[s.m_right_bottom].m_nCap[5] -= delta;
			m_nodes[s.m_right_bottom].m_nCap[6] += delta; m_nodes[s.m_right_top].m_nCap[2] -= delta;
		}
		delta = t_old - t;
		t_p += delta;
		t_q -= delta;
	}

	// now separate flow that goes through topVG and flow that goes trough real t-links of p,q
	t_p -= w_p;
	t_q -= w_q;

	w_p -= flow_pq;
	w_q += flow_pq;

	if (w_p<0 && t_p>0)
	{
		*left_top_flow = MIN<TotalFlowType>(-w_p, t_p);
	}
	else if (w_p>0 && t_p<0)
	{
		*left_top_flow = -MIN<TotalFlowType>(w_p, -t_p);
	}
	else
	{
		*left_top_flow = 0;
	}
	m_nodes[s.m_left_top].m_tCap = t_p - *left_top_flow;

	if (w_q<0 && t_q>0)
	{
		*right_top_flow = MIN<TotalFlowType>(-w_q, t_q);
	}
	else if (w_q>0 && t_q<0)
	{
		*right_top_flow = -MIN<TotalFlowType>(w_q, -t_q);
	}
	else
	{
		*right_top_flow = 0;
	}
	m_nodes[s.m_right_top].m_tCap = t_q - *right_top_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowBottomInsideSquare(Square& s, VirtualGraph2 *bottomVG)
{
	TotalFlowType w_p, w_q, w_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, delta;

	bottomVG->ConstructGraph(w_p, w_q, w_pq);

	t_p = m_nodes[s.m_left_bottom].m_tCap + w_p;
	t_q = m_nodes[s.m_right_bottom].m_tCap + w_q;

	// Push flow through through n-links within square s
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	if (t_p>0 && t_q<0)
	{
		t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[6]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_top].m_nCap[1]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_bottom].m_tCap -= delta;
			m_nodes[s.m_left_bottom].m_nCap[6] -= delta; m_nodes[s.m_left_top].m_nCap[2] += delta;
			m_nodes[s.m_left_top].m_nCap[1] -= delta; m_nodes[s.m_right_bottom].m_nCap[5] += delta;
			m_nodes[s.m_right_bottom].m_tCap += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[7]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[2]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_bottom].m_tCap -= delta;
			m_nodes[s.m_left_bottom].m_nCap[7] -= delta; m_nodes[s.m_right_top].m_nCap[3] += delta;
			m_nodes[s.m_right_top].m_nCap[2] -= delta; m_nodes[s.m_right_bottom].m_nCap[6] += delta;
			m_nodes[s.m_right_bottom].m_tCap += delta;
		}
	}
	else if (t_p<0 && t_q>0)
	{
		t = MIN<TotalFlowType>(-t_p, t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[2]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[5]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_bottom].m_tCap += delta;
			m_nodes[s.m_left_bottom].m_nCap[6] += delta; m_nodes[s.m_left_top].m_nCap[2] -= delta;
			m_nodes[s.m_left_top].m_nCap[1] += delta; m_nodes[s.m_right_bottom].m_nCap[5] -= delta;
			m_nodes[s.m_right_bottom].m_tCap -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_top].m_nCap[3]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[6]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_bottom].m_tCap += delta;
			m_nodes[s.m_left_bottom].m_nCap[7] += delta; m_nodes[s.m_right_top].m_nCap[3] -= delta;
			m_nodes[s.m_right_top].m_nCap[2] += delta; m_nodes[s.m_right_bottom].m_nCap[6] -= delta;
			m_nodes[s.m_right_bottom].m_tCap -= delta;
		}
	}
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowTopInsideSquare(Square& s, VirtualGraph2 *topVG)
{
	TotalFlowType w_p, w_q, w_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, delta;

	topVG->ConstructGraph(w_p, w_q, w_pq);

	t_p = m_nodes[s.m_left_top].m_tCap + w_p;
	t_q = m_nodes[s.m_right_top].m_tCap + w_q;

	// Push flow through through n-links within square s
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	if (t_p>0 && t_q<0)
	{
		t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[2]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_left_bottom].m_nCap[7]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_tCap -= delta;
			m_nodes[s.m_left_top].m_nCap[2] -= delta; m_nodes[s.m_left_bottom].m_nCap[6] += delta;
			m_nodes[s.m_left_bottom].m_nCap[7] -= delta; m_nodes[s.m_right_top].m_nCap[3] += delta;
			m_nodes[s.m_right_top].m_tCap += delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_top].m_nCap[1]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_bottom].m_nCap[6]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_tCap -= delta;
			m_nodes[s.m_left_top].m_nCap[1] -= delta; m_nodes[s.m_right_bottom].m_nCap[5] += delta;
			m_nodes[s.m_right_bottom].m_nCap[6] -= delta; m_nodes[s.m_right_top].m_nCap[2] += delta;
			m_nodes[s.m_right_top].m_tCap += delta;
		}
	}
	else if (t_p<0 && t_q>0)
	{
		t = MIN<TotalFlowType>(-t_p, t_q);
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_left_bottom].m_nCap[6]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[3]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_tCap += delta;
			m_nodes[s.m_left_top].m_nCap[2] += delta; m_nodes[s.m_left_bottom].m_nCap[6] -= delta;
			m_nodes[s.m_left_bottom].m_nCap[7] += delta; m_nodes[s.m_right_top].m_nCap[3] -= delta;
			m_nodes[s.m_right_top].m_tCap -= delta;
		}
		delta = MIN<TotalFlowType>(t, m_nodes[s.m_right_bottom].m_nCap[5]);
		delta = MIN<TotalFlowType>(delta, m_nodes[s.m_right_top].m_nCap[2]);
		if (delta > 0)
		{
			t -= delta;
			m_nodes[s.m_left_top].m_tCap += delta;
			m_nodes[s.m_left_top].m_nCap[1] += delta; m_nodes[s.m_right_bottom].m_nCap[5] -= delta;
			m_nodes[s.m_right_bottom].m_nCap[6] += delta; m_nodes[s.m_right_top].m_nCap[2] -= delta;
			m_nodes[s.m_right_top].m_tCap -= delta;
		}
	}
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeFlowForTwoVirtualGraphs(NodeId p, NodeId q, VirtualGraph2* vg0, VirtualGraph2* vg, TotalFlowType* p0_flow, TotalFlowType* q0_flow)
{
	TotalFlowType w0_p, w0_q, w0_pq, flow0_pq;
	TotalFlowType w_p, w_q, w_pq;
	TotalFlowType t_p, t_q;
	TotalFlowType t, t_old, delta;

	vg0->ConstructGraph(w0_p, w0_q, w0_pq);
	vg->ConstructGraph(w_p, w_q, w_pq);
	flow0_pq = 0; // flow through edge with capacity w0_pq

	t_p = m_nodes[p].m_tCap + w_p + w0_p;
	t_q = m_nodes[q].m_tCap + w_q + w0_q;

	// Push flow through through edge w_pq
	// assuming that t-link capacities of nodes p,q are t_p, t_q.
	// Calculate the amount of flow that must be pushed through w0_pq (but do not push it).
	if (t_p>0 && t_q<0)
	{
		t_old = t = MIN<TotalFlowType>(t_p, -t_q);
		delta = MIN<TotalFlowType>(t, w_pq);
		t -= delta;

		flow0_pq = MIN<TotalFlowType>(t, w0_pq);

		delta = t_old - t;
		t_p -= delta;
		t_q += delta;
	}

	// now separate flow that goes through vg0 and flow that goes trough real t-links of p,q
	t_p -= w0_p;
	t_q -= w0_q;

	w0_p -= flow0_pq;
	w0_q += flow0_pq;

	if (w0_p<0 && t_p>0)
	{
		*p0_flow = MIN<TotalFlowType>(-w0_p, t_p);
		m_nodes[p].m_tCap -= *p0_flow;
	}
	else if (w0_p>0 && t_p<0)
	{
		*p0_flow = -MIN<TotalFlowType>(w0_p, -t_p);
		m_nodes[p].m_tCap -= *p0_flow;
	}
	else
	{
		*p0_flow = 0;
	}

	if (w0_q<0 && t_q>0)
	{
		*q0_flow = MIN<TotalFlowType>(-w0_q, t_q);
		m_nodes[q].m_tCap -= *q0_flow;
	}
	else if (w0_q>0 && t_q<0)
	{
		*q0_flow = -MIN<TotalFlowType>(w0_q, -t_q);
		m_nodes[q].m_tCap -= *q0_flow;
	}
	else
	{
		*q0_flow = 0;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
// Functions PushFlow<Side>()
///////////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::PushFlowRight(Square& s, TotalFlowType right_top_flow, TotalFlowType right_bottom_flow)
{
	TotalFlowType delta;

	if (right_top_flow > 0)
	{
		delta = MIN<TotalFlowType>(right_top_flow, m_nodes[s.m_right_top].m_nCap[4]);
		right_top_flow -= delta;
		m_nodes[s.m_right_top].m_nCap[4] -= delta; m_nodes[s.m_left_top].m_nCap[0] += delta;
		m_nodes[s.m_left_top].m_tCap += delta;
		ASSERT(right_top_flow <= m_nodes[s.m_right_top].m_nCap[3]);
	}
	else if (right_top_flow < 0)
	{
		delta = MIN<TotalFlowType>(-right_top_flow, m_nodes[s.m_left_top].m_nCap[0]);
		right_top_flow += delta;
		m_nodes[s.m_right_top].m_nCap[4] += delta; m_nodes[s.m_left_top].m_nCap[0] -= delta;
		m_nodes[s.m_left_top].m_tCap -= delta;
		ASSERT(-right_top_flow <= m_nodes[s.m_left_bottom].m_nCap[7]);
	}
	m_nodes[s.m_right_top].m_nCap[3] -= right_top_flow; m_nodes[s.m_left_bottom].m_nCap[7] += right_top_flow;
	m_nodes[s.m_left_bottom].m_tCap += right_top_flow;

	if (right_bottom_flow > 0)
	{
		delta = MIN<TotalFlowType>(right_bottom_flow, m_nodes[s.m_right_bottom].m_nCap[5]);
		right_bottom_flow -= delta;
		m_nodes[s.m_right_bottom].m_nCap[5] -= delta; m_nodes[s.m_left_top].m_nCap[1] += delta;
		m_nodes[s.m_left_top].m_tCap += delta;
		ASSERT(right_bottom_flow <= m_nodes[s.m_right_bottom].m_nCap[4]);
	}
	else if (right_bottom_flow < 0)
	{
		delta = MIN<TotalFlowType>(-right_bottom_flow, m_nodes[s.m_left_top].m_nCap[1]);
		right_bottom_flow += delta;
		m_nodes[s.m_right_bottom].m_nCap[5] += delta; m_nodes[s.m_left_top].m_nCap[1] -= delta;
		m_nodes[s.m_left_top].m_tCap -= delta;
		ASSERT(-right_bottom_flow <= m_nodes[s.m_left_bottom].m_nCap[0]);
	}
	m_nodes[s.m_right_bottom].m_nCap[4] -= right_bottom_flow; m_nodes[s.m_left_bottom].m_nCap[0] += right_bottom_flow;
	m_nodes[s.m_left_bottom].m_tCap += right_bottom_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::PushFlowBottom(Square& s, TotalFlowType left_bottom_flow, TotalFlowType right_bottom_flow)
{
	TotalFlowType delta;

	if (left_bottom_flow > 0)
	{
		delta = MIN<TotalFlowType>(left_bottom_flow, m_nodes[s.m_left_bottom].m_nCap[6]);
		left_bottom_flow -= delta;
		m_nodes[s.m_left_bottom].m_nCap[6] -= delta; m_nodes[s.m_left_top].m_nCap[2] += delta;
		m_nodes[s.m_left_top].m_tCap += delta;
		ASSERT(left_bottom_flow <= m_nodes[s.m_left_bottom].m_nCap[7]);
	}
	else if (left_bottom_flow < 0)
	{
		delta = MIN<TotalFlowType>(-left_bottom_flow, m_nodes[s.m_left_top].m_nCap[2]);
		left_bottom_flow += delta;
		m_nodes[s.m_left_bottom].m_nCap[6] += delta; m_nodes[s.m_left_top].m_nCap[2] -= delta;
		m_nodes[s.m_left_top].m_tCap -= delta;
		ASSERT(-left_bottom_flow <= m_nodes[s.m_right_top].m_nCap[3]);
	}
	m_nodes[s.m_left_bottom].m_nCap[7] -= left_bottom_flow; m_nodes[s.m_right_top].m_nCap[3] += left_bottom_flow;
	m_nodes[s.m_right_top].m_tCap += left_bottom_flow;

	if (right_bottom_flow > 0)
	{
		delta = MIN<TotalFlowType>(right_bottom_flow, m_nodes[s.m_right_bottom].m_nCap[5]);
		right_bottom_flow -= delta;
		m_nodes[s.m_right_bottom].m_nCap[5] -= delta; m_nodes[s.m_left_top].m_nCap[1] += delta;
		m_nodes[s.m_left_top].m_tCap += delta;
		ASSERT(right_bottom_flow <= m_nodes[s.m_right_bottom].m_nCap[6]);
	}
	else if (right_bottom_flow < 0)
	{
		delta = MIN<TotalFlowType>(-right_bottom_flow, m_nodes[s.m_left_top].m_nCap[1]);
		right_bottom_flow += delta;
		m_nodes[s.m_right_bottom].m_nCap[5] += delta; m_nodes[s.m_left_top].m_nCap[1] -= delta;
		m_nodes[s.m_left_top].m_tCap -= delta;
		ASSERT(-right_bottom_flow <= m_nodes[s.m_right_top].m_nCap[2]);
	}
	m_nodes[s.m_right_bottom].m_nCap[6] -= right_bottom_flow; m_nodes[s.m_right_top].m_nCap[2] += right_bottom_flow;
	m_nodes[s.m_right_top].m_tCap += right_bottom_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::PushFlowLeft(Square& s, TotalFlowType left_top_flow, TotalFlowType left_bottom_flow)
{
	TotalFlowType delta;

	if (left_top_flow > 0)
	{
		delta = MIN<TotalFlowType>(left_top_flow, m_nodes[s.m_left_top].m_nCap[0]);
		left_top_flow -= delta;
		m_nodes[s.m_left_top].m_nCap[0] -= delta; m_nodes[s.m_right_top].m_nCap[4] += delta;
		m_nodes[s.m_right_top].m_tCap += delta;
		ASSERT(left_top_flow <= m_nodes[s.m_left_top].m_nCap[1]);
	}
	else if (left_top_flow < 0)
	{
		delta = MIN<TotalFlowType>(-left_top_flow, m_nodes[s.m_right_top].m_nCap[4]);
		left_top_flow += delta;
		m_nodes[s.m_left_top].m_nCap[0] += delta; m_nodes[s.m_right_top].m_nCap[4] -= delta;
		m_nodes[s.m_right_top].m_tCap -= delta;
		ASSERT(-left_top_flow <= m_nodes[s.m_right_bottom].m_nCap[5]);
	}
	m_nodes[s.m_left_top].m_nCap[1] -= left_top_flow; m_nodes[s.m_right_bottom].m_nCap[5] += left_top_flow;
	m_nodes[s.m_right_bottom].m_tCap += left_top_flow;

	if (left_bottom_flow > 0)
	{
		delta = MIN<TotalFlowType>(left_bottom_flow, m_nodes[s.m_left_bottom].m_nCap[7]);
		left_bottom_flow -= delta;
		m_nodes[s.m_left_bottom].m_nCap[7] -= delta; m_nodes[s.m_right_top].m_nCap[3] += delta;
		m_nodes[s.m_right_top].m_tCap += delta;
		ASSERT(left_bottom_flow <= m_nodes[s.m_left_bottom].m_nCap[0]);
	}
	else if (left_bottom_flow < 0)
	{
		delta = MIN<TotalFlowType>(-left_bottom_flow, m_nodes[s.m_right_top].m_nCap[3]);
		left_bottom_flow += delta;
		m_nodes[s.m_left_bottom].m_nCap[7] += delta; m_nodes[s.m_right_top].m_nCap[3] -= delta;
		m_nodes[s.m_right_top].m_tCap -= delta;
		ASSERT(-left_bottom_flow <= m_nodes[s.m_right_bottom].m_nCap[4]);
	}
	m_nodes[s.m_left_bottom].m_nCap[0] -= left_bottom_flow; m_nodes[s.m_right_bottom].m_nCap[4] += left_bottom_flow;
	m_nodes[s.m_right_bottom].m_tCap += left_bottom_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::PushFlowTop(Square& s, TotalFlowType left_top_flow, TotalFlowType right_top_flow)
{
	TotalFlowType delta;

	if (left_top_flow > 0)
	{
		delta = MIN<TotalFlowType>(left_top_flow, m_nodes[s.m_left_top].m_nCap[2]);
		left_top_flow -= delta;
		m_nodes[s.m_left_top].m_nCap[2] -= delta; m_nodes[s.m_left_bottom].m_nCap[6] += delta;
		m_nodes[s.m_left_bottom].m_tCap += delta;
		ASSERT(left_top_flow <= m_nodes[s.m_left_top].m_nCap[1]);
	}
	else if (left_top_flow < 0)
	{
		delta = MIN<TotalFlowType>(-left_top_flow, m_nodes[s.m_left_bottom].m_nCap[6]);
		left_top_flow += delta;
		m_nodes[s.m_left_top].m_nCap[2] += delta; m_nodes[s.m_left_bottom].m_nCap[6] -= delta;
		m_nodes[s.m_left_bottom].m_tCap -= delta;
		ASSERT(-left_top_flow <= m_nodes[s.m_right_bottom].m_nCap[5]);
	}
	m_nodes[s.m_left_top].m_nCap[1] -= left_top_flow; m_nodes[s.m_right_bottom].m_nCap[5] += left_top_flow;
	m_nodes[s.m_right_bottom].m_tCap += left_top_flow;

	if (right_top_flow > 0)
	{
		delta = MIN<TotalFlowType>(right_top_flow, m_nodes[s.m_right_top].m_nCap[3]);
		right_top_flow -= delta;
		m_nodes[s.m_right_top].m_nCap[3] -= delta; m_nodes[s.m_left_bottom].m_nCap[7] += delta;
		m_nodes[s.m_left_bottom].m_tCap += delta;
		ASSERT(right_top_flow <= m_nodes[s.m_right_top].m_nCap[2]);
	}
	else if (right_top_flow < 0)
	{
		delta = MIN<TotalFlowType>(-right_top_flow, m_nodes[s.m_left_bottom].m_nCap[7]);
		right_top_flow += delta;
		m_nodes[s.m_right_top].m_nCap[3] += delta; m_nodes[s.m_left_bottom].m_nCap[7] -= delta;
		m_nodes[s.m_left_bottom].m_tCap -= delta;
		ASSERT(-right_top_flow <= m_nodes[s.m_right_bottom].m_nCap[6]);
	}
	m_nodes[s.m_right_top].m_nCap[2] -= right_top_flow; m_nodes[s.m_right_bottom].m_nCap[6] += right_top_flow;
	m_nodes[s.m_right_bottom].m_tCap += right_top_flow;
}

template <class NLinkFlowType, class TotalFlowType>
  inline void GridGraph8<NLinkFlowType,TotalFlowType>::ComputeMaximumFlow2(NodeId p, NodeId q, int e_pq, int e_qp)
{
	ASSERT(q==GetNeib(p, e_pq) && p==GetNeib(q, e_qp));

	TotalFlowType delta;

	if (m_nodes[p].m_tCap>0 && m_nodes[q].m_tCap<0)
	{
		delta = MIN<TotalFlowType>(m_nodes[p].m_tCap, -m_nodes[q].m_tCap);
		delta = MIN<TotalFlowType>(delta, m_nodes[p].m_nCap[e_pq]);
		m_nodes[p].m_tCap -= delta;
		m_nodes[p].m_nCap[e_pq] -= delta;
		m_nodes[q].m_nCap[e_qp] += delta;
		m_nodes[q].m_tCap += delta;
	}
	else if (m_nodes[p].m_tCap<0 && m_nodes[q].m_tCap>0)
	{
		delta = MIN<TotalFlowType>(-m_nodes[p].m_tCap, m_nodes[q].m_tCap);
		delta = MIN<TotalFlowType>(delta, m_nodes[q].m_nCap[e_qp]);
		m_nodes[q].m_tCap -= delta;
		m_nodes[q].m_nCap[e_qp] -= delta;
		m_nodes[p].m_nCap[e_pq] += delta;
		m_nodes[p].m_tCap += delta;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
//                             Main DP functions                                         //
///////////////////////////////////////////////////////////////////////////////////////////

template <class NLinkFlowType, class TotalFlowType>
  HRESULT GridGraph8<NLinkFlowType,TotalFlowType>::DynamicProgrammingMaxflowHorz()
{
	ASSERT(IsAllocated());
	if (m_sizeX < 2)
	{
		return S_OK;
	}

	int y;
	Square s;
	TotalFlowType top_flow, bottom_flow;
	VirtualGraph2* buf;
	VirtualGraph2* leftVG;
	VirtualGraph2* rightVG;

	buf =  new(std::nothrow) VirtualGraph2[m_sizeX];
	if (buf == NULL)
	{
		return E_OUTOFMEMORY;
	}

	for (y=0; y<m_sizeY-1; y++)
	{
		// forward pass
		leftVG = buf;
		rightVG = leftVG + 1;
		leftVG->Clear();
		for (SetSquare(0, y, s); ; MoveSquareRight(s))
		{
			ComputeVirtualGraphFromLeftToRight(s, leftVG, rightVG);

			if (s.m_x == m_sizeX-2)
			{
				break;
			}

			leftVG = rightVG;
			rightVG ++;
		}

		// backward pass
		VirtualGraph2 rightBorderVG;
		rightBorderVG.Set(0, m_nodes[s.m_right_top].m_nCap[2], m_nodes[s.m_right_bottom].m_nCap[6], 0);
		ComputeFlowForTwoVirtualGraphs(s.m_right_top, s.m_right_bottom, rightVG, &rightBorderVG, &top_flow, &bottom_flow);
		ComputeMaximumFlow2(s.m_right_top, s.m_right_bottom, 2, 6);

		while ( 1 )
		{
			PushFlowRight(s, top_flow, bottom_flow);
			ComputeFlowLeft(s, leftVG, &top_flow, &bottom_flow);

			if (s.m_x == 0)
			{
				ASSERT(top_flow == 0 && bottom_flow == 0);
				break;
			}

			MoveSquareLeft(s);
			leftVG --;
		}
	}

	delete [] buf;

	return S_OK;
}

template <class NLinkFlowType, class TotalFlowType>
  HRESULT GridGraph8<NLinkFlowType,TotalFlowType>::DynamicProgrammingMaxflowVert()
{
	ASSERT(IsAllocated());
	if (m_sizeY < 2)
	{
		return S_OK;
	}

	int x;
	Square s;
	TotalFlowType left_flow, right_flow;
	VirtualGraph2* buf;
	VirtualGraph2* topVG;
	VirtualGraph2* bottomVG;

	buf =  new(std::nothrow) VirtualGraph2[m_sizeY];
	if (buf == NULL)
	{
		return E_OUTOFMEMORY;
	}

	for (x=0; x<m_sizeX-1; x++)
	{
		// forward pass
		topVG = buf;
		bottomVG = topVG + 1;
		topVG->Clear();
		for (SetSquare(x, 0, s); ; MoveSquareDown(s))
		{
			ComputeVirtualGraphFromTopToBottom(s, topVG, bottomVG);

			if (s.m_y == m_sizeY-2)
			{
				break;
			}

			topVG = bottomVG;
			bottomVG ++;
		}

		// backward pass
		VirtualGraph2 bottomBorderVG;
		bottomBorderVG.Set(0, m_nodes[s.m_left_bottom].m_nCap[0], m_nodes[s.m_right_bottom].m_nCap[4], 0);
		ComputeFlowForTwoVirtualGraphs(s.m_left_bottom, s.m_right_bottom, bottomVG, &bottomBorderVG, &left_flow, &right_flow);
		ComputeMaximumFlow2(s.m_left_bottom, s.m_right_bottom, 0, 4);

		while ( 1 )
		{
			PushFlowBottom(s, left_flow, right_flow);
			ComputeFlowTop(s, topVG, &left_flow, &right_flow);

			if (s.m_y == 0)
			{
				ASSERT(left_flow == 0 && right_flow == 0);
				break;
			}

			MoveSquareUp(s);
			topVG --;
		}
	}
	delete [] buf;

	return S_OK;
}

template <class NLinkFlowType, class TotalFlowType>
  HRESULT GridGraph8<NLinkFlowType,TotalFlowType>::MultipleTreesDynamicProgrammingMaxflowHorz()
{
	ASSERT(IsAllocated());
	if (m_sizeX < 2 || m_sizeY < 4)
	{
		return S_OK;
	}

	int y;
	Square s;
	TotalFlowType left_flow, right_flow, top_flow, bottom_flow;
	VirtualGraph2* bufHorz;
	VirtualGraph2* bufVert;
	VirtualGraph2* leftVG;
	VirtualGraph2* rightVG;
	VirtualGraph2* topVG;
	VirtualGraph2* bottomVG;

	bufHorz =  new(std::nothrow) VirtualGraph2[m_sizeX*(m_sizeY/2)];
	if (bufHorz == NULL)
	{
		return E_OUTOFMEMORY;
	}
	bufVert =  new(std::nothrow) VirtualGraph2[m_sizeY];
	if (bufVert == NULL)
	{
		delete [] bufHorz;
		return E_OUTOFMEMORY;
	}

	// pass 1 - initializing bufHorz and bufVert
	topVG = bufVert;
	bottomVG = topVG + 1;
	topVG->Clear();
	for (y=0; ; y+=2)
	{
		// forward pass
		leftVG = bufHorz + m_sizeX*(y/2);
		rightVG = leftVG + 1;
		leftVG->Clear();
		for (SetSquare(0, y, s); s.m_x<m_sizeX-2; MoveSquareRight(s))
		{
			ComputeVirtualGraphFromLeftToRight(s, leftVG, rightVG);

			leftVG = rightVG;
			rightVG ++;
		}

		rightVG->Clear();

		if (y >= m_sizeY - 3) // it is the last line
		{
			break;
		}

		ComputeVirtualGraphFromLeftTopToBottom(s, leftVG, topVG, bottomVG);
		topVG = bottomVG;
		bottomVG ++;

		MoveSquareDown(s);
		ComputeVirtualGraphFromTopToBottom(s, topVG, bottomVG);
		topVG = bottomVG;
		bottomVG ++;
	}
	bottomVG->Clear();

	// pass 2
	while ( 1 )
	{
		bool lastPass;
	
		if (s.m_x == 0)
		{
			lastPass = true; // do last pass on a colum [0;1] x [0; m_sizeY)
		}
		else
		{
			lastPass = false;
		}

		if (((m_sizeX - s.m_x) & 1) == 0)
		{
			// move up
			do
			{
				VirtualGraph2 tmpVG;
				ComputeVirtualGraphFromRightLeftToTop(s, rightVG, leftVG, &tmpVG);
				tmpVG.m_E01 += m_nodes[s.m_left_top].m_nCap[0];
				tmpVG.m_E10 += m_nodes[s.m_right_top].m_nCap[4];
				ComputeFlowForTwoVirtualGraphs(s.m_left_top, s.m_right_top, topVG, &tmpVG, &left_flow, &right_flow);

				MoveSquareUp(s); // moving up
				topVG --;        //

				PushFlowBottom(s, left_flow, right_flow);
				ComputeFlowTopInsideSquare(s, topVG);

				if (lastPass)
				{
					tmpVG.m_E01 -= m_nodes[s.m_left_bottom].m_nCap[0];
					tmpVG.m_E10 -= m_nodes[s.m_right_bottom].m_nCap[4];
					ComputeVirtualGraphFromBottomToTop(s, &tmpVG, topVG);

					MoveSquareUp(s);      //
					bottomVG = topVG;     //
					topVG --;             // moving up
					leftVG -= m_sizeX;    //
					rightVG = leftVG + 1; //
				}
				else
				{
					MoveSquareDown(s); // moving down
					topVG ++;          //

					ComputeVirtualGraphFromRightToLeft(s, rightVG, leftVG);

					MoveSquareLeft(s); //
					rightVG = leftVG;  // moving left
					leftVG --;         //

					ComputeVirtualGraphFromRightBottomLeftToTop(s, rightVG, bottomVG, leftVG, topVG);

					MoveSquareUp(s);  //
					bottomVG = topVG; // moving up
					topVG --;         //

					ComputeVirtualGraphFromBottomToTop(s, bottomVG, topVG);

					MoveSquareUp(s);      //
					bottomVG = topVG;     //
					topVG --;             // moving up
					leftVG -= m_sizeX;    //
					rightVG = leftVG + 1; //

					ComputeVirtualGraphFromLeftBottomToRight(s, leftVG, bottomVG, rightVG);

					MoveSquareRight(s); //
					leftVG = rightVG;   // moving right
					rightVG ++;         //
				}
			} while (s.m_y > 1);
		}
		else
		{
			// move down
			do
			{
				VirtualGraph2 tmpVG;
				ComputeVirtualGraphFromRightLeftToBottom(s, rightVG, leftVG, &tmpVG);
				tmpVG.m_E01 += m_nodes[s.m_left_bottom].m_nCap[0];
				tmpVG.m_E10 += m_nodes[s.m_right_bottom].m_nCap[4];
				ComputeFlowForTwoVirtualGraphs(s.m_left_bottom, s.m_right_bottom, bottomVG, &tmpVG, &left_flow, &right_flow);

				MoveSquareDown(s);  // moving down
				bottomVG ++;        //

				PushFlowTop(s, left_flow, right_flow);
				ComputeFlowBottomInsideSquare(s, bottomVG);

				if (lastPass)
				{
					tmpVG.m_E01 -= m_nodes[s.m_left_top].m_nCap[0];
					tmpVG.m_E10 -= m_nodes[s.m_right_top].m_nCap[4];
					ComputeVirtualGraphFromTopToBottom(s, &tmpVG, bottomVG);

					MoveSquareDown(s);    //
					topVG = bottomVG;     //
					bottomVG ++;          // moving down
					leftVG += m_sizeX;    //
					rightVG = leftVG + 1; //
				}
				else
				{
					MoveSquareUp(s); // moving up
					bottomVG --;     //

					ComputeVirtualGraphFromRightToLeft(s, rightVG, leftVG);

					MoveSquareLeft(s); //
					rightVG = leftVG;  // moving left
					leftVG --;         //

					ComputeVirtualGraphFromRightLeftTopToBottom(s, rightVG, leftVG, topVG, bottomVG);

					MoveSquareDown(s); //
					topVG = bottomVG;  // moving down
					bottomVG ++;       //

					ComputeVirtualGraphFromTopToBottom(s, topVG, bottomVG);

					MoveSquareDown(s);    //
					topVG = bottomVG;     //
					bottomVG ++;          // moving down
					leftVG += m_sizeX;    //
					rightVG = leftVG + 1; //

					ComputeVirtualGraphFromLeftTopToRight(s, leftVG, topVG, rightVG);

					MoveSquareRight(s); //
					leftVG = rightVG;   // moving right
					rightVG ++;         //
				}
			} while (s.m_y < m_sizeY-3);
		}

		if (lastPass)
		{
			break;
		}

		ComputeVirtualGraphFromRightToLeft(s, rightVG, leftVG);

		MoveSquareLeft(s); //
		rightVG = leftVG;  // moving left
		leftVG --;         //
	}

	// pass 3
	for (y=0; y<m_sizeY-1; y+=2)
	{
		leftVG = bufHorz + m_sizeX*(y/2);
		rightVG = leftVG + 1;
		leftVG->Clear();
		SetSquare(0, y, s);

		ComputeVirtualGraphFromRightToLeft(s, rightVG, leftVG);

		VirtualGraph2 leftBorderVG;
		leftBorderVG.Set(0, m_nodes[s.m_left_top].m_nCap[2], m_nodes[s.m_left_bottom].m_nCap[6], 0);
		ComputeFlowForTwoVirtualGraphs(s.m_left_top, s.m_left_bottom, leftVG, &leftBorderVG, &top_flow, &bottom_flow);
		ComputeMaximumFlow2(s.m_left_top, s.m_left_bottom, 2, 6);

		while ( 1 )
		{
			PushFlowLeft(s, top_flow, bottom_flow);
			ComputeFlowRight(s, rightVG, &top_flow, &bottom_flow);

			if (s.m_x == m_sizeX-2)
			{
				ASSERT(top_flow == 0 && bottom_flow == 0);
				break;
			}

			MoveSquareRight(s);
			rightVG ++;
		}
	}

	delete [] bufHorz;
	delete [] bufVert;

	return S_OK;
}
