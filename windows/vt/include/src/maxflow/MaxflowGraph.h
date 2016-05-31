/**************************************************************************\
*
* Copyright (c) 2004, Microsoft Corp.  All Rights Reserved.
*
* Module Name:
*
*   MaxflowGraph.h
*
* Abstract:
*
*   MSRC Image Editing public header file for max flow graph classes
*
\**************************************************************************/

#pragma once

#ifdef CUTOUT_EXPORTS
#define CUTOUT_API
#else
#define CUTOUT_API
#endif


#include <new>
#include <assert.h>
#ifndef ASSERT
#define ASSERT assert
#endif

namespace GraphcutAPI
{

// Callback function specifically for these classes - combines the "upper"
// abort and progress callbacks, but avoids pushing any notion of time down
// to the graph classes. I deliberately don't use the existing objects
// because (a) I don't want dependencies "up" and (b) I want to emphasise
// the notion that the count (~time) might not monotonically decrease.
// Usage: the graph core routines will call this, passing an estimated count
//        of number of <somethings> done so far - doesn't matter what these
//        are, but the intention is that they monotonically increase. The
//        second argument is whatever was handed in at construction time.
//        The routine is expected to return true if processing should be
//        interrupted.
typedef bool (*ProgressCallback)( int numSoFar, void* callbackData );

/****************************************************************************
*
*   Maxflow algorithm for general graphs
*
****************************************************************************/

template <class NLinkFlowType, class TotalFlowType> 
  class CUTOUT_API GeneralGraph
{
    struct Node; // forward declaration
public:

    typedef NLinkFlowType T_NLinkFlowType;
    typedef TotalFlowType TLinkFlowType; // capacities of t-links. Should be the same as TotalFlowType

    typedef Node* NodeId;

    // Default constructor. Creates a graph with no nodes and no edges.
    // After constructor has been called, IsAllocated() can be used to check whether allocation was successfull
    GeneralGraph(
                ProgressCallback abortCallback = NULL,
                void* callbackData = NULL
                );
    // Destructor
    ~GeneralGraph();
    // Creates a graph with no nodes and no edges.
    // If it is already allocated, it is deallocated first
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Allocate();
    bool IsAllocated()
    {
        return m_isAllocated;
    }
    void DeAllocate();

    // Adds a node to the graph
    // Returns NULL if not enough memory, and a valid NodeId otherwise
    NodeId AddNode();

    // Adds an edge source->p with weight w (if w>0) or p->sink with weight -w (if w<0)
    void AddTEdge(NodeId p, TLinkFlowType w);
    // Adds edges p->q and q->p with weights w_pq and w_qp, respectively
    HRESULT AddNEdges(NodeId p, NodeId q, NLinkFlowType w_pq, NLinkFlowType w_qp);
    // Computes maxflow. Can be called several times. Returns S_OK if success and E_OUTOFMEMORY otherwise, S_FALSE on abort
    HRESULT Maxflow();
    // After Maxflow() has been called, GetSegmentation(p) can be used to determine the label of p
    // Returns 0 if p is with the source, 1 otherwise
    int GetSegmentation(NodeId p);

/**************************************************************************\
* Implementation
\**************************************************************************/
private:

    /////////////////////////////////////////////////////////////////////////
    //                             Structures                              //
    /////////////////////////////////////////////////////////////////////////

    struct Edge;
    struct NodeBlock;
    struct EdgeBlock;

    struct NodeQueue
    {
        NodeQueue();
        ~NodeQueue();
        HRESULT Allocate(int size = 128); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the queue

        HRESULT EnQueue(Node* p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        Node* DeQueue();

    private:
        // queue contains at most (m_end - m_start - 1) elements.
        // If m_top == m_rear, then queue is empty
        Node**      m_start; // array of size m_end - m_start
        Node**      m_end; 
        Node**      m_top; // pointer to top element in the queue
        Node**      m_rear; // pointer to next element after the last element in the queue
    };

    struct NodeStack
    {
        NodeStack();
        ~NodeStack();
        HRESULT Allocate(int size = 128); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the stack

        HRESULT Push(Node* p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        Node* Pop();

    private:
        Node**      m_start; // array of size m_end - m_start
        Node**      m_end; 
        Node**      m_current; // pointer to first element which is NOT in the stack
    };


    /////////////////////////////////////////////////////////////////////////
    //                           Data members                              //
    /////////////////////////////////////////////////////////////////////////

    static Edge* PARENT_TERMINAL;
    static Edge* PARENT_ORPHAN;
    static Edge* PARENT_FREE;

    bool        m_isAllocated;

    NodeBlock*  m_nodeBlockFirst;
    EdgeBlock*  m_edgeBlockFirst;

    // Queue of active nodes. m_ActiveQueueTop is the first node in the list,
    // m_ActiveQueueRear is the last. Node p is in the queue iff p->m_next != NULL.
    // The last node (m_ActiveQueueRear) points to itself.
    // Note: some nodes in the queue may be passive (Node::m_isActive = 0); however,
    // if a node is active, then it is in the queue
    Node*       m_activeQueueTop;
    Node*       m_activeQueueRear;
    int         m_activeQueueCount; // only used in debug mode

    // Queue of orphans.
    NodeQueue   m_orphanQueue;

    // Stack of orphans.
    NodeStack   m_orphanStack;

    /////////////////////////////////////////////////////////////////////////
    //                          Member functions                           //
    /////////////////////////////////////////////////////////////////////////

    // Returns path (defined by an edge from source tree to the sink tree) if path is found and NULL if there are no more augmenting paths
    Edge* Grow();
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Augment(Edge* boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Adopt(Node *p);

    void ActiveQueueReset(); // empties the queue. In addition, all pointers Node::m_next must be set to NULL
    void ActiveQueueAdd(Node* p); // add p to the rear of the queue
    Node* ActiveQueueGetTop();
    void ActiveQueueRemoveTop();

    ProgressCallback m_abortCallback;
    void* m_callbackData;

    // Prevent copy & assign
private:
    GeneralGraph( const GeneralGraph& );
    GeneralGraph& operator=( const GeneralGraph& );
}; // end class GeneralGraph

/****************************************************************************
*
*   Maxflow algorithm for regular grid graphs with 8-neighborhood system
*   Algorithm is based on the paper
*     "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision",
*     Y. Boykov and V. Kolmogorov, PAMI'04,
*   and my PhD thesis. Additionally, speed-up techniques based on dynamic programming are implemented.
*
*
*
*   Nodes are referenced using struct GridGraph8::NodeId. Nodes are stored in a row-by-row order.
*   Below is an example of fast pass over nodes:
*
*   GridGraph8::NodeId p = graph->GetNodeId(0, 0);
*   for (int y=0; y<graph->GetSizeY(); y++)
*   for (int x=0; x<graph->GetSizeX(); x++)
*   {
*       ... (your code; p corresponds to node (x, y))
*       
*       p ++;
*   }
*
*   Edges are referenced using integers from 0 to 7
*   (going in clockwise direction: 0 corresponds to shift (1,0), 1 to (1,1), and so on)
*
****************************************************************************/

template <class NLinkFlowType, class TotalFlowType> 
  class CUTOUT_API GridGraph8
{
    struct Node; // Forward decl
public:

    static const int NeighborhoodSize = 8;

    typedef NLinkFlowType T_NLinkFlowType;
    typedef TotalFlowType TLinkFlowType; // capacities of t-links. Should be the same as TotalFlowType

    typedef __int32 NodeId;

    NodeId GetNodeId(int x, int y)
    {
        ASSERT(x>=0 && x<m_sizeX && y>=0 && y<m_sizeY);
        return x + y*m_sizeX;
    };

    int GetX(NodeId p) // very slow
    {
        ASSERT(p>=0 && p<m_nodeNum);
        return (int)(p % m_sizeX);
    }

    int GetY(NodeId p) // very slow
    {
        ASSERT(p>=0 && p<m_nodeNum);
        return (int)(p / m_sizeX);
    }

    // Default constructor
    GridGraph8(ProgressCallback abortCallback = NULL,
               void* callbackData = NULL
               );
    // Constructor. All edge weights are set to 0
    GridGraph8(int sizeX, int sizeY,
                ProgressCallback abortCallback = NULL,
                void* callbackData = NULL
              );
    // Destructor
    ~GridGraph8();
    // Initializes graph of size sizeX x sizeY. All edge weights are set to zero.
    // If graph is already allocated, it is deallocated first.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise
    HRESULT Allocate(int sizeX, int sizeY);
    bool IsAllocated()
    {
        return (m_nodes != NULL);
    }
    void DeAllocate();

    // Get grid dimensions
    int GetSizeX()
    {
        ASSERT(IsAllocated());
        return m_sizeX;
    }
    int GetSizeY()
    {
        ASSERT(IsAllocated());
        return m_sizeY;
    }
    // Adds an edge source->p with weight w (if w>0) or p->sink with weight -w (if w<0)
    void AddTEdge(NodeId p, TLinkFlowType w);
    // Adds an edge p->q with weight w
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdge(NodeId p, int edge, NLinkFlowType w);
    // Adds edges p->q and q->p with weights w and w_rev, respectively
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdges(NodeId p, int edge, NLinkFlowType w, NLinkFlowType w_rev);
    // returns the weight of TEdge at node p
    TLinkFlowType GetTEdge(NodeId p);
    // returns the weight of NEdge node p and specified edge 
    NLinkFlowType GetNEdge(NodeId p, int edge);
    // Computes maxflow. Can be called several times.
    // Optional parameter 'flag' (flag=0,1,2 or 3) is an estimation of how much flow can be pushed through the graph:
    // 0 means large flow, 3 small flow.
    // This flags helps to choose the fastest algorithm: for flag=0,1 or 2 initial flow is computed
    // via dynamic programming.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise, S_FALSE if interrupted
    HRESULT Maxflow(int flag = 0);
    // After Maxflow() has been called, GetSegmentation(p) can be used to determine the label of p
    // Returns 0 if p is with the source, 1 otherwise
    int GetSegmentation(NodeId p);


/**************************************************************************\
* Implementation
\**************************************************************************/
private:

    struct Node
    {
        NodeId          m_next; // queue of active nodes
        NLinkFlowType   m_nCap[NeighborhoodSize];
        TLinkFlowType   m_tCap; // if m_tCap > 0, then residual capacity of edge source->node is m_tCap
                                // if m_tCap < 0, then residual capacity of edge node->sink is -m_tCap
        unsigned        m_parent : 4;
        unsigned        m_tree : 1; // 0 = source, 1 = sink (if parent!=PARENT_FREE)
        unsigned        m_isActive : 1;
    };

    static const unsigned PARENT_TERMINAL;
    static const unsigned PARENT_ORPHAN;
    static const unsigned PARENT_FREE;

    /////////////////////////////////////////////////////////////////////////
    //                             Structures                              //
    /////////////////////////////////////////////////////////////////////////

    struct NodeQueue
    {
        NodeQueue();
        ~NodeQueue();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the queue

        HRESULT EnQueue(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId DeQueue();

    private:
        // queue constains at most (m_end - m_start - 1) elements.
        // If m_top == m_rear, then queue is empty
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_top; // pointer to top element in the queue
        NodeId*      m_rear; // pointer to next element after the last element in the queue
    };

    struct NodeStack
    {
        NodeStack();
        ~NodeStack();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the stack

        HRESULT Push(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId Pop();

    private:
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_current; // pointer to first element which is NOT in the stack
    };


    /////////////////////////////////////////////////////////////////////////
    //                           Data members                              //
    /////////////////////////////////////////////////////////////////////////

    int         m_nodeShifts[NeighborhoodSize]; // for functions GetNeib() and GetAndCheckNeib()

    int         m_sizeX, m_sizeY; // image dimensions
    int         m_nodeNum; // = m_sizeX*m_sizeY

    Node*       m_nodes; // array of size m_nodeNum

    // Queue of active nodes. m_ActiveQueueTop is the first node in the list,
    // m_ActiveQueueRear is the last. Node p is in the queue iff p->m_next != NULL.
    // The last node (m_ActiveQueueRear) points to itself.
    // Note: some nodes in the queue may be passive (Node::m_isActive = 0); however,
    // if a node is active, then it is in the queue
    NodeId      m_activeQueueTop;
    NodeId      m_activeQueueRear;
    int         m_activeQueueCount; // Only used in debug mode

    // Queue of orphans.
    NodeQueue   m_orphanQueue;

    // Stack of orphans.
    NodeStack   m_orphanStack;

    /////////////////////////////////////////////////////////////////////////
    //                          Member functions                           //
    /////////////////////////////////////////////////////////////////////////

	bool IS_VALID_ID(NodeId p) { return (p >= 0); }

    int GetReverseEdge(int e)
    {
        return e ^ 4;
    };

    NodeId GetNeib(NodeId p, int e)
    {
        return (NodeId)(p + m_nodeShifts[e]);
    };

    NodeId GetAndCheckNeib(NodeId p, int e)
    {
        NodeId q = (NodeId)(p + m_nodeShifts[e]);
        return (q<m_nodeNum) ? q : -1;
    };

    ///////////////////////////////////////////////////////////////////////////////////
    //                        Maxflow via dynamic programming                        //
    ///////////////////////////////////////////////////////////////////////////////////

    struct Square;

    void SetSquare(int x, int y, Square& s);
    void MoveSquareRight(Square& s);
    void MoveSquareDown(Square& s);
    void MoveSquareLeft(Square& s);
    void MoveSquareUp(Square& s);

    // Structure VirtualGraph2 defines a graph with two main nodes p, q (and perhaps with other implicit nodes).
    // The graph is defined not by a set of edges, but rather by costs of minimum cuts for every possible
    // configuration of nodes p1, p2:
    //      m_E01 is the cost of the minimum cut when node p is with the source and node q is with the sink
    //      m_E10 is the cost of the minimum cut when node p is with the sink and node q is with the source
    //      m_E11 is the cost of the minimum cut when nodes p, q are with the sink
    // m_E01, m_E10 and m_E11 are normalized (by adding a constant) so that the cost of the minimum cut
    // for configuration '00' is zero.
    struct VirtualGraph2;

    // Functions below take a subgraph in square 's', merge it with given input virtual graphs,
    // and compute the virtual graph on the output side.
    // By convention, virtual graph with main nodes p and q does not include edges
    // from p,q to the terminals and edge p--q.
    // For example, ComputeVirtualGraphFromLeftToRight() computes virtual graph with main nodes
    // s.m_right_top and s.m_right_bottom summarizing graph containing the following:
    //   - graph leftVG with main nodes s.m_left_top and s.m_left_bottom
    //   - all edges withing square s (including edges to terminals), except
    //     for edges s.m_right_top--terminal, s.m_right_bottom--terminal and s.m_right_top--s.m_right_bottom.
    void ComputeVirtualGraphFromRightToLeft(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG);

    void ComputeVirtualGraphFromTopToBottom(Square& s, VirtualGraph2* topVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromRightLeftToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromLeftTopToBottom(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromRightLeftTopToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG);

    void ComputeVirtualGraphFromLeftToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* rightVG);
    void ComputeVirtualGraphFromLeftBottomToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* bottomVG, VirtualGraph2* rightVG);
    void ComputeVirtualGraphFromLeftTopToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* rightVG);

    void ComputeVirtualGraphFromBottomToTop(Square& s, VirtualGraph2* bottomVG, VirtualGraph2* topVG);
    void ComputeVirtualGraphFromRightLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG);
    void ComputeVirtualGraphFromRightBottomLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* bottomVG, VirtualGraph2* leftVG, VirtualGraph2* topVG);

    // Functions below compute maximum flow in a small graph.
    // The last word of function name is square side which defines two main nodes p,q
    // (e.g. for side 'left' p=m_left_top, q=m_left_bottom).
    // The graph is as follows: it is virtual graph sideVG with main nodes p,q
    // plus all edges within square 's' involving at least one of the nodes p,q.
    // Residual capacities of edges within square 's' are modified accordingly.
    // Note that some flow must also be pushed through graph defined by sideVG. Since it is not
    // given by residual capacities, flow that must be pushed through this graph is
    // returned in the last two arguments. For example, for ComputeFlowLeft()
    // left_top_flow is set to the amount of flow that must be pushed from node s.m_left_top,
    // and left_bottom_flow is set to the amount of flow that must be pushed from node s.m_left_bottom.
    // This operation must be done by the caller.
    // 
    // Note that the flow returned in the last two arguments can be either positive or negative.
    // Positive flow means that it must be pushed from the node, negative - that the [negated]
    // flow must be pushed to the node
    void ComputeFlowLeft(Square& s, VirtualGraph2* leftVG, TotalFlowType* left_top_flow, TotalFlowType* left_bottom_flow);
    void ComputeFlowRight(Square& s, VirtualGraph2* rightVG, TotalFlowType* right_top_flow, TotalFlowType* right_bottom_flow);
    void ComputeFlowTop(Square& s, VirtualGraph2* topVG, TotalFlowType* left_top_flow, TotalFlowType* right_top_flow);

    // Functions below are similar to functions computeFlow<Side>() except
    // that they are not concerned how the flow goes in the graph
    //   sideVG + t-links of p,q + edge between p,q
    // (where p,q are nodes of the side 'Side')
    // The remaining flow is sent through t-links of p,q
    void ComputeFlowBottomInsideSquare(Square& s, VirtualGraph2* bottomVG);
    void ComputeFlowTopInsideSquare(Square& s, VirtualGraph2* topVG);

    // The function takes two virtual graphs vg0 and vg on nodes p,q and t-links of nodes p,q,
    // merges them together. Flow that must be pushed through vg0 is returned in the last
    // two arguments (similar to ComputeFlow<Side>() functions). T-link capacity of node p
    // is reduced by *p0_flow and t-link capacity of node q is reduced by *q0_flow.
    void ComputeFlowForTwoVirtualGraphs(NodeId p, NodeId q, VirtualGraph2* vg0, VirtualGraph2* vg, TotalFlowType* p0_flow, TotalFlowType* q0_flow);

    // Functions below push given amount of flow from nodes of the given side of
    // square 's' (side=last word in the function name) through edges within 's'
    // going to the other side.
    // For example, PushFlowRight() pushes 'right_top_flow' amount of flow from s.m_right_top
    // through edges s.m_right_top--s.m_left_top and s.m_right_top--s.m_left_bottom,
    // and 'right_bottom_flow' amount of flow from s.m_right_bottom
    // through edges s.m_right_bottom--s.m_left_top and s.m_right_bottom--s.m_left_bottom.
    // Residual capacities of edges must be large enough to hold the amount of flow needed.
    void PushFlowLeft(Square& s, TotalFlowType left_top_flow, TotalFlowType left_bottom_flow);
    void PushFlowBottom(Square& s, TotalFlowType left_bottom_flow, TotalFlowType right_bottom_flow);
    void PushFlowRight(Square& s, TotalFlowType right_top_flow, TotalFlowType right_bottom_flow);
    void PushFlowTop(Square& s, TotalFlowType left_top_flow, TotalFlowType right_top_flow);

    // computes maximum flow on a graph with two nodes - p and q
    void ComputeMaximumFlow2(NodeId p, NodeId q, int e_pq, int e_qp);

    // main DP functions

    // Do DP on subgraphs of size 'm_sizeX' x 2
    HRESULT DynamicProgrammingMaxflowHorz();
    // Do DP on subgraphs of size 2 x 'm_sizeY'
    HRESULT DynamicProgrammingMaxflowVert();

    HRESULT MultipleTreesDynamicProgrammingMaxflowHorz();

    ///////////////////////////////////////////////////////////////////////////////////
    //                                Main algorithm                                 //
    ///////////////////////////////////////////////////////////////////////////////////

    // Returns true if path is found and false if there are no more augmenting paths
    // If path is found, then source_node is set to be a node in the source set,
    // sink_node - node in the sink set, and boundaryEdge - edge from source_node to sink_node
    bool Grow(NodeId& source_node, NodeId& sink_node, int& boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Augment(NodeId source_node, NodeId sink_node, int boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Adopt(NodeId p);

    void ActiveQueueReset(); // empties the queue. In addition, all pointers Node::m_next must be set to NULL
    void ActiveQueueAdd(NodeId p); // add p to the rear of the queue
    NodeId ActiveQueueGetTop();
    void ActiveQueueRemoveTop();

    ProgressCallback m_abortCallback;
    void* m_callbackData;

    // Prevent copy & assign
private:
    GridGraph8( const GridGraph8& );
    GridGraph8& operator=( const GridGraph8& );
}; // end class GridGraph8

/****************************************************************************
*
*   Maxflow algorithm for regular grid graphs with 4-neighborhood system
*   Algorithm is based on the paper
*     "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision",
*     Y. Boykov and V. Kolmogorov, PAMI'04,
*   and my PhD thesis. Additionally, speed-up techniques based on dynamic programming are implemented.
*
*
*
*   Nodes are referenced using struct GridGraph4::NodeId. Nodes are stored in a row-by-row order.
*   Below is an example of fast pass over nodes:
*
*   GridGraph4::NodeId p = graph->GetNodeId(0, 0);
*   for (int y=0; y<graph->GetSizeY(); y++)
*   for (int x=0; x<graph->GetSizeX(); x++)
*   {
*       ... (your code; p corresponds to node (x, y))
*       
*       p ++;
*   }
*
*   Edges are referenced using integers from 0 to 4
*   (going in clockwise direction: 0 corresponds to shift (1,0), 1 to (0,1), and so on)
*
****************************************************************************/

template <class NLinkFlowType, class TotalFlowType> 
  class CUTOUT_API GridGraph4
{
    struct Node; // Forward decl
public:

    static const int NeighborhoodSize = 4;

    typedef NLinkFlowType T_NLinkFlowType;
    typedef TotalFlowType TLinkFlowType; // capacities of t-links. Should be the same as TotalFlowType

    typedef __int32 NodeId;

    NodeId GetNodeId(int x, int y)
    {
        ASSERT(x>=0 && x<m_sizeX && y>=0 && y<m_sizeY);
        return x + y*m_sizeX;
    };

    int GetX(NodeId p) // very slow
    {
        ASSERT(p>=0 && p<m_nodeNum);
        return (int)(p % m_sizeX);
    }

    int GetY(NodeId p) // very slow
    {
        ASSERT(p>=0 && p<m_nodeNum);
        return (int)(p / m_sizeX);
    }

    // Default constructor
    GridGraph4(ProgressCallback abortCallback = NULL,
               void* callbackData = NULL
               );
    // Constructor. All edge weights are set to 0
    GridGraph4(int sizeX, int sizeY,
                ProgressCallback abortCallback = NULL,
                void* callbackData = NULL
              );
    // Destructor
    ~GridGraph4();
    // Initializes graph of size sizeX x sizeY. All edge weights are set to zero.
    // If graph is already allocated, it is deallocated first.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise
    HRESULT Allocate(int sizeX, int sizeY);
    bool IsAllocated()
    {
        return (m_nodes != NULL);
    }
    void DeAllocate();

    // Get grid dimensions
    int GetSizeX()
    {
        ASSERT(IsAllocated());
        return m_sizeX;
    }
    int GetSizeY()
    {
        ASSERT(IsAllocated());
        return m_sizeY;
    }
    // Adds an edge source->p with weight w (if w>0) or p->sink with weight -w (if w<0)
    void AddTEdge(NodeId p, TLinkFlowType w);
    // Adds an edge p->q with weight w
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdge(NodeId p, int edge, NLinkFlowType w);
    // Adds edges p->q and q->p with weights w and w_rev, respectively
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdges(NodeId p, int edge, NLinkFlowType w, NLinkFlowType w_rev);
    // returns the weight of TEdge at node p
    TLinkFlowType GetTEdge(NodeId p);
    // returns the weight of NEdge node p and specified edge 
    NLinkFlowType GetNEdge(NodeId p, int edge);
    // Computes maxflow. Can be called several times.
    // Optional parameter 'flag' (flag=0,1,2 or 3) is an estimation of how much flow can be pushed through the graph:
    // 0 means large flow, 3 small flow.
    // This flags helps to choose the fastest algorithm: for flag=0,1 or 2 initial flow is computed
    // via dynamic programming.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise, S_FALSE if interrupted
    HRESULT Maxflow(int flag = 0);
    // After Maxflow() has been called, GetSegmentation(p) can be used to determine the label of p
    // Returns 0 if p is with the source, 1 otherwise
    int GetSegmentation(NodeId p);


/**************************************************************************\
* Implementation
\**************************************************************************/
private:

    struct Node
    {
        NodeId          m_next; // queue of active nodes
        NLinkFlowType   m_nCap[NeighborhoodSize];
        TLinkFlowType   m_tCap; // if m_tCap > 0, then residual capacity of edge source->node is m_tCap
                                // if m_tCap < 0, then residual capacity of edge node->sink is -m_tCap
        unsigned        m_parent : 4;
        unsigned        m_tree : 1; // 0 = source, 1 = sink (if parent!=PARENT_FREE)
        unsigned        m_isActive : 1;
    };

    static const unsigned PARENT_TERMINAL;
    static const unsigned PARENT_ORPHAN;
    static const unsigned PARENT_FREE;

    /////////////////////////////////////////////////////////////////////////
    //                             Structures                              //
    /////////////////////////////////////////////////////////////////////////

    struct NodeQueue
    {
        NodeQueue();
        ~NodeQueue();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the queue

        HRESULT EnQueue(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId DeQueue();

    private:
        // queue constains at most (m_end - m_start - 1) elements.
        // If m_top == m_rear, then queue is empty
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_top; // pointer to top element in the queue
        NodeId*      m_rear; // pointer to next element after the last element in the queue
    };

    struct NodeStack
    {
        NodeStack();
        ~NodeStack();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the stack

        HRESULT Push(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId Pop();

    private:
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_current; // pointer to first element which is NOT in the stack
    };


    /////////////////////////////////////////////////////////////////////////
    //                           Data members                              //
    /////////////////////////////////////////////////////////////////////////

    int         m_nodeShifts[NeighborhoodSize]; // for functions GetNeib() and GetAndCheckNeib()

    int         m_sizeX, m_sizeY; // image dimensions
	int         m_nodeNum; // = m_sizeX*m_sizeY

    Node*       m_nodes; // array of size m_nodeNum

    // Queue of active nodes. m_ActiveQueueTop is the first node in the list,
    // m_ActiveQueueRear is the last. Node p is in the queue iff p->m_next != NULL.
    // The last node (m_ActiveQueueRear) points to itself.
    // Note: some nodes in the queue may be passive (Node::m_isActive = 0); however,
    // if a node is active, then it is in the queue
    NodeId      m_activeQueueTop;
    NodeId      m_activeQueueRear;
    int         m_activeQueueCount; // Only used in debug mode

    // Queue of orphans.
    NodeQueue   m_orphanQueue;

    // Stack of orphans.
    NodeStack   m_orphanStack;

    /////////////////////////////////////////////////////////////////////////
    //                          Member functions                           //
    /////////////////////////////////////////////////////////////////////////

	bool IS_VALID_ID(NodeId p) { return (p >= 0); }

    int GetReverseEdge(int e)
    {
        return e ^ 2;
    };

    NodeId GetNeib(NodeId p, int e)
    {
        return (NodeId)(p + m_nodeShifts[e]);
    };

    NodeId GetAndCheckNeib(NodeId p, int e)
    {
        NodeId q = (NodeId)(p + m_nodeShifts[e]);
        return (q<m_nodeNum) ? q : -1;
    };

    ///////////////////////////////////////////////////////////////////////////////////
    //                        Maxflow via dynamic programming                        //
    ///////////////////////////////////////////////////////////////////////////////////

    struct Square;

    void SetSquare(int x, int y, Square& s);
    void MoveSquareRight(Square& s);
    void MoveSquareDown(Square& s);
    void MoveSquareLeft(Square& s);
    void MoveSquareUp(Square& s);

    // Structure VirtualGraph2 defines a graph with two main nodes p, q (and perhaps with other implicit nodes).
    // The graph is defined not by a set of edges, but rather by costs of minimum cuts for every possible
    // configuration of nodes p1, p2:
    //      m_E01 is the cost of the minimum cut when node p is with the source and node q is with the sink
    //      m_E10 is the cost of the minimum cut when node p is with the sink and node q is with the source
    //      m_E11 is the cost of the minimum cut when nodes p, q are with the sink
    // m_E01, m_E10 and m_E11 are normalized (by adding a constant) so that the cost of the minimum cut
    // for configuration '00' is zero.
    struct VirtualGraph2;

    // Functions below take a subgraph in square 's', merge it with given input virtual graphs,
    // and compute the virtual graph on the output side.
    // By convention, virtual graph with main nodes p and q does not include edges
    // from p,q to the terminals and edge p--q.
    // For example, ComputeVirtualGraphFromLeftToRight() computes virtual graph with main nodes
    // s.m_right_top and s.m_right_bottom summarizing graph containing the following:
    //   - graph leftVG with main nodes s.m_left_top and s.m_left_bottom
    //   - all edges withing square s (including edges to terminals), except
    //     for edges s.m_right_top--terminal, s.m_right_bottom--terminal and s.m_right_top--s.m_right_bottom.
    void ComputeVirtualGraphFromRightToLeft(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG);

    void ComputeVirtualGraphFromTopToBottom(Square& s, VirtualGraph2* topVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromRightLeftToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromLeftTopToBottom(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG);
    void ComputeVirtualGraphFromRightLeftTopToBottom(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* bottomVG);

    void ComputeVirtualGraphFromLeftToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* rightVG);
    void ComputeVirtualGraphFromLeftBottomToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* bottomVG, VirtualGraph2* rightVG);
    void ComputeVirtualGraphFromLeftTopToRight(Square& s, VirtualGraph2* leftVG, VirtualGraph2* topVG, VirtualGraph2* rightVG);

    void ComputeVirtualGraphFromBottomToTop(Square& s, VirtualGraph2* bottomVG, VirtualGraph2* topVG);
    void ComputeVirtualGraphFromRightLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* leftVG, VirtualGraph2* topVG);
    void ComputeVirtualGraphFromRightBottomLeftToTop(Square& s, VirtualGraph2* rightVG, VirtualGraph2* bottomVG, VirtualGraph2* leftVG, VirtualGraph2* topVG);

    // Functions below compute maximum flow in a small graph.
    // The last word of function name is square side which defines two main nodes p,q
    // (e.g. for side 'left' p=m_left_top, q=m_left_bottom).
    // The graph is as follows: it is virtual graph sideVG with main nodes p,q
    // plus all edges within square 's' involving at least one of the nodes p,q.
    // Residual capacities of edges within square 's' are modified accordingly.
    // Note that some flow must also be pushed through graph defined by sideVG. Since it is not
    // given by residual capacities, flow that must be pushed through this graph is
    // returned in the last two arguments. For example, for ComputeFlowLeft()
    // left_top_flow is set to the amount of flow that must be pushed from node s.m_left_top,
    // and left_bottom_flow is set to the amount of flow that must be pushed from node s.m_left_bottom.
    // This operation must be done by the caller.
    // 
    // Note that the flow returned in the last two arguments can be either positive or negative.
    // Positive flow means that it must be pushed from the node, negative - that the [negated]
    // flow must be pushed to the node
    void ComputeFlowLeft(Square& s, VirtualGraph2* leftVG, TotalFlowType* left_top_flow, TotalFlowType* left_bottom_flow);
    void ComputeFlowRight(Square& s, VirtualGraph2* rightVG, TotalFlowType* right_top_flow, TotalFlowType* right_bottom_flow);
    void ComputeFlowTop(Square& s, VirtualGraph2* topVG, TotalFlowType* left_top_flow, TotalFlowType* right_top_flow);

    // The function takes two virtual graphs vg0 and vg on nodes p,q and t-links of nodes p,q,
    // merges them together. Flow that must be pushed through vg0 is returned in the last
    // two arguments (similar to ComputeFlow<Side>() functions). T-link capacity of node p
    // is reduced by *p0_flow and t-link capacity of node q is reduced by *q0_flow.
    void ComputeFlowForTwoVirtualGraphs(NodeId p, NodeId q, VirtualGraph2* vg0, VirtualGraph2* vg, TotalFlowType* p0_flow, TotalFlowType* q0_flow);

    // Functions below push given amount of flow from nodes of the given side of
    // square 's' (side=last word in the function name) through edges within 's'
    // going to the other side.
    // For example, PushFlowRight() pushes 'right_top_flow' amount of flow from s.m_right_top
    // through edges s.m_right_top--s.m_left_top and s.m_right_top--s.m_left_bottom,
    // and 'right_bottom_flow' amount of flow from s.m_right_bottom
    // through edges s.m_right_bottom--s.m_left_top and s.m_right_bottom--s.m_left_bottom.
    // Residual capacities of edges must be large enough to hold the amount of flow needed.
    void PushFlowLeft(Square& s, TotalFlowType left_top_flow, TotalFlowType left_bottom_flow);
    void PushFlowBottom(Square& s, TotalFlowType left_bottom_flow, TotalFlowType right_bottom_flow);
    void PushFlowRight(Square& s, TotalFlowType right_top_flow, TotalFlowType right_bottom_flow);
    void PushFlowTop(Square& s, TotalFlowType left_top_flow, TotalFlowType right_top_flow);

    // computes maximum flow on a graph with two nodes - p and q
    void ComputeMaximumFlow2(NodeId p, NodeId q, int e_pq, int e_qp);

    // main DP functions

    // Do DP on subgraphs of size 'm_sizeX' x 2
    HRESULT DynamicProgrammingMaxflowHorz();
    // Do DP on subgraphs of size 2 x 'm_sizeY'
    HRESULT DynamicProgrammingMaxflowVert();

    HRESULT MultipleTreesDynamicProgrammingMaxflowHorz();

    ///////////////////////////////////////////////////////////////////////////////////
    //                                Main algorithm                                 //
    ///////////////////////////////////////////////////////////////////////////////////

    // Returns true if path is found and false if there are no more augmenting paths
    // If path is found, then source_node is set to be a node in the source set,
    // sink_node - node in the sink set, and boundaryEdge - edge from source_node to sink_node
    bool Grow(NodeId& source_node, NodeId& sink_node, int& boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Augment(NodeId source_node, NodeId sink_node, int boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Adopt(NodeId p);

    void ActiveQueueReset(); // empties the queue. In addition, all pointers Node::m_next must be set to NULL
    void ActiveQueueAdd(NodeId p); // add p to the rear of the queue
    NodeId ActiveQueueGetTop();
    void ActiveQueueRemoveTop();

    ProgressCallback m_abortCallback;
    void* m_callbackData;

    // Prevent copy & assign
private:
    GridGraph4( const GridGraph4& );
    GridGraph4& operator=( const GridGraph4& );
}; // end class GridGraph4


/****************************************************************************
*
*   Maxflow algorithm for a collection of regular grid graphs with 4-neighborhood system
*   connected with a small number of extra edges.
*   Algorithm is based on the paper
*     "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Computer Vision",
*     Y. Boykov and V. Kolmogorov, PAMI'04,
*   and my PhD thesis. Additionally, speed-up techniques based on dynamic programming are implemented.
*
*
*
*   Nodes are referenced using struct GridGraph4::NodeId. Nodes are stored in a row-by-row order.
*   Below is an example of fast pass over nodes:
*
*   for (int t=0; t<graph->GetGridNum(); t++)
*   {
*		GridGraph4Plus::NodeId p = graph->GetNodeId(t, 0, 0);
*		for (int y=0; y<graph->GetSizeY(t); y++)
*		for (int x=0; x<graph->GetSizeX(t); x++)
*		{
*			... (your code; p corresponds to node (x, y) in grid t)
*			
*			p ++;
*		}
*   }
*
*   Edges are referenced using integers from 0 to 4
*   (going in clockwise direction: 0 corresponds to shift (1,0), 1 to (0,1), and so on)
*
****************************************************************************/


template <class NLinkFlowType, class TotalFlowType> 
  class CUTOUT_API GridGraph4Plus
{
    struct Node; // Forward decl
public:

    static const int NeighborhoodSize = 4;

    typedef NLinkFlowType T_NLinkFlowType;
    typedef TotalFlowType TLinkFlowType; // capacities of t-links. Should be the same as TotalFlowType

	typedef unsigned short GridId; // an integer in [0,grid_num-1]. Stored at every node, so use 2 bytes to save space
    typedef __int32 NodeId;

    NodeId GetNodeId(GridId t, int x, int y)
    {
        ASSERT(t>=0 && t<m_gridNum && x>=0 && x<m_sizeX[t] && y>=0 && y<m_sizeY[t]);
        return m_nodes[t] + x + y*m_sizeX[t];
    };

    // Default constructor
    GridGraph4Plus(ProgressCallback abortCallback = NULL,
               void* callbackData = NULL
               );
    // Constructor. gridNum must be positive. sizeX and sizeY are arrays of size gridNum. All edge weights are set to 0.
    GridGraph4Plus(int gridNum, int* sizeX, int* sizeY, int extraEdgeNumMax,
                ProgressCallback abortCallback = NULL,
                void* callbackData = NULL
              );
    // Destructor
    ~GridGraph4Plus();
    // Initializes graph. All edge weights are set to zero.
    // If graph is already allocated, it is deallocated first.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise
    HRESULT Allocate(int gridNum, int* sizeX, int* sizeY, int extraEdgeNumMax);
    bool IsAllocated()
    {
        return (m_nodes != NULL);
    }
    void DeAllocate();

    // Get grid dimensions
    int GetSizeX(GridId t)
    {
        ASSERT(IsAllocated());
        return m_sizeX[t];
    }
    int GetSizeY(GridId t)
    {
        ASSERT(IsAllocated());
        return m_sizeY[t];
    }
    // Adds an edge source->p with weight w (if w>0) or p->sink with weight -w (if w<0)
    void AddTEdge(NodeId p, TLinkFlowType w);
    // Adds an edge p->q with weight w
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdge(NodeId p, int edge, NLinkFlowType w);
    // Adds edges p->q and q->p with weights w and w_rev, respectively
    // where q is the neighbor of p determined by e.
    // q must an inner node
    void AddNEdges(NodeId p, int edge, NLinkFlowType w, NLinkFlowType w_rev);
	// Adds edges p->q and q->p (p!=q) between arbitrary nodes p and q. Restrictions:
	// - Can be called at most extraEdgeNumMax times.
	// - The number of extra edges per node cannot exceed 256-7.
	// - Cannot be called after Maxflow().
    void AddExtraEdges(NodeId p, NodeId q, NLinkFlowType w, NLinkFlowType w_rev);
    // returns the weight of TEdge at node p
    TLinkFlowType GetTEdge(NodeId p);
    // returns the weight of NEdge node p and specified edge 
    NLinkFlowType GetNEdge(NodeId p, int edge);
    // Computes maxflow. Can be called several times.
    // Returns S_OK if success and E_OUTOFMEMORY otherwise, S_FALSE if interrupted
    HRESULT Maxflow();
    // After Maxflow() has been called, GetSegmentation(p) can be used to determine the label of p
    // Returns 0 if p is with the source, 1 otherwise
    int GetSegmentation(NodeId p);


/**************************************************************************\
* Implementation
\**************************************************************************/
private:

	typedef __int32 ExtraArcId;

    struct Node
    {
		ExtraArcId		m_extraArcs; // node i has  ((i+1)->m_extraArcs) - i->m_extraArcs) extra arcs
        NodeId          m_next; // queue of active nodes
        NLinkFlowType   m_nCap[NeighborhoodSize];
        TLinkFlowType   m_tCap; // if m_tCap > 0, then residual capacity of edge source->node is m_tCap
                                // if m_tCap < 0, then residual capacity of edge node->sink is -m_tCap
		GridId			m_t;
		unsigned		m_parent : 8;
		unsigned        m_tree : 1; // 0 = source, 1 = sink (if parent!=PARENT_FREE)
		unsigned        m_isActive : 1;
    };

	struct ExtraEdge0 // used only during graph construction. Deallocated by InitExtraArcs()
	{
		NodeId			m_i[2];
		NLinkFlowType	m_nCap[2];
	};

	struct ExtraArc
	{
		NodeId			m_head;
		ExtraArcId		m_sister;
		NLinkFlowType	m_nCap;
	};

    static const unsigned PARENT_TERMINAL;
    static const unsigned PARENT_ORPHAN;
    static const unsigned PARENT_FREE;
    static const unsigned PARENT_EXTRA_ARC;
	// possible values for e=p->m_parent :
	//   [0,NeighborhoodSize-1]  - node neighbors in the grid
	//   [PARENT_EXTRA_ARC, 255] - correspond to extra arcs ( m_extraArcs[p->m_extraArcs + e - PARENT_EXTRA_ARC] )
	//   PARENT_TERMINAL, PARENT_ORPHAN, PARENT_EXTRA_ARC - special values

    /////////////////////////////////////////////////////////////////////////
    //                             Structures                              //
    /////////////////////////////////////////////////////////////////////////

    struct NodeQueue
    {
        NodeQueue();
        ~NodeQueue();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the queue

        HRESULT EnQueue(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId DeQueue();

    private:
        // queue constains at most (m_end - m_start - 1) elements.
        // If m_top == m_rear, then queue is empty
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_top; // pointer to top element in the queue
        NodeId*      m_rear; // pointer to next element after the last element in the queue
    };

    struct NodeStack
    {
        NodeStack();
        ~NodeStack();
        HRESULT Allocate(int size = 1); // returns S_OK if success and E_OUTOFMEMORY otherwise
        void DeAllocate();

        void Reset(); // empties the stack

        HRESULT Push(NodeId p); // returns S_OK if success and E_OUTOFMEMORY otherwise
        NodeId Pop();

    private:
        NodeId*      m_start; // array of size m_end - m_start
        NodeId*      m_end; 
        NodeId*      m_current; // pointer to first element which is NOT in the stack
    };


    /////////////////////////////////////////////////////////////////////////
    //                           Data members                              //
    /////////////////////////////////////////////////////////////////////////

	int			m_gridNum;
    int*        m_sizeX; // grid dimensions
	int*		m_sizeY;

    Node*       m_nodesArray; // array of nodes
	int			m_allocatedNodeNum; // size of m_nodesArray. Some of the nodes are dummy.
    NodeId*     m_nodes; // m_nodes[t] is points to the start of the array's index (within m_nodesArray) for grid t


	int			m_extraEdgeNum, m_extraEdgeNumMax;
	ExtraEdge0*	m_extraEdges0;      // Deallocated by InitExtraEdges()
	ExtraArc*	m_extraArcs;        // Allocated by InitExtraEdges()

    int*        m_nodeShifts; // for function GetNeib(). m_nodeShifts[NeighborhoodSize*t + e].

	// Queue of active nodes. m_ActiveQueueTop is the first node in the list,
    // m_ActiveQueueRear is the last. Node p is in the queue iff p->m_next != NULL.
    // The last node (m_ActiveQueueRear) points to itself.
    // Note: some nodes in the queue may be passive (Node::m_isActive = 0); however,
    // if a node is active, then it is in the queue
    NodeId      m_activeQueueTop;
    NodeId      m_activeQueueRear;
    int         m_activeQueueCount; // Only used in debug mode

    // Queue of orphans.
    NodeQueue   m_orphanQueue;

    // Stack of orphans.
    NodeStack   m_orphanStack;

    /////////////////////////////////////////////////////////////////////////
    //                          Member functions                           //
    /////////////////////////////////////////////////////////////////////////

	HRESULT InitExtraArcs(); // called during the first call to Maxflow().

	bool IS_VALID_ID(NodeId p) { return (p >= 0); }

	unsigned GetReverseEdge(int e)
    {
        return (unsigned)(e ^ 2);
    };

    NodeId GetNeib(NodeId p, int e)
    {
        return p + m_nodeShifts[NeighborhoodSize*m_nodesArray[p].m_t+e];
    };

    ///////////////////////////////////////////////////////////////////////////////////
    //                                Main algorithm                                 //
    ///////////////////////////////////////////////////////////////////////////////////

    // Returns true if path is found and false if there are no more augmenting paths
    // If path is found, then source_node is set to be a node in the source set,
    // sink_node - node in the sink set, and boundaryEdge - edge from source_node to sink_node
    bool Grow(NodeId& source_node, NodeId& sink_node, int& boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Augment(NodeId source_node, NodeId sink_node, int boundaryEdge);
    // Returns S_OK or E_OUTOFMEMORY
    HRESULT Adopt(NodeId p);

    void ActiveQueueReset(); // empties the queue. In addition, all pointers Node::m_next must be set to NULL
    void ActiveQueueAdd(NodeId p); // add p to the rear of the queue
    NodeId ActiveQueueGetTop();
    void ActiveQueueRemoveTop();

    ProgressCallback m_abortCallback;
    void* m_callbackData;

    // Prevent copy & assign
private:
    GridGraph4Plus( const GridGraph4Plus& );
    GridGraph4Plus& operator=( const GridGraph4Plus& );
}; // end class GridGraph4Plus

} // end namespace GraphcutAPI

