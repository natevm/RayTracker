// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
//! \file   cyBVH.h 
//! \author Cem Yuksel
//! 
//! \brief  Bounding Volume Hierarchy class.
//!
//! BVH is a storage class for Bounding Volume Hierarchies.
//!
//-------------------------------------------------------------------------------
// 
// Copyright (c) 2016, Cem Yuksel <cem@cemyuksel.com>
// All rights reserved.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy 
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is 
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
// SOFTWARE.
// 
//-------------------------------------------------------------------------------

#ifndef _CY_BVH_H_INCLUDED_
#define _CY_BVH_H_INCLUDED_

//-------------------------------------------------------------------------------
namespace cy {
//-------------------------------------------------------------------------------

#ifndef CY_BVH_ELEMENT_COUNT_BITS
#define CY_BVH_ELEMENT_COUNT_BITS	3	//!< Determines the number of bits needed to represent the maximum number of elements in a node (8)
#endif

#ifndef CY_BVH_MAX_ELEMENT_COUNT
#define CY_BVH_MAX_ELEMENT_COUNT	(1<<CY_BVH_ELEMENT_COUNT_BITS)	//!< Determines the maximum number of elements in a node (8)
#endif

#define _CY_BVH_NODE_DATA_BITS		(sizeof(unsigned int)*8)
#define _CY_BVH_ELEMENT_COUNT_MASK	((1<<CY_BVH_ELEMENT_COUNT_BITS)-1)
#define _CY_BVH_LEAF_BIT_MASK		((unsigned int)1<<(_CY_BVH_NODE_DATA_BITS-1))
#define _CY_BVH_CHILD_INDEX_BITS	(_CY_BVH_NODE_DATA_BITS-1)
#define _CY_BVH_CHILD_INDEX_MASK	(_CY_BVH_LEAF_BIT_MASK-1)
#define _CY_BVH_ELEMENT_OFFSET_BITS	(_CY_BVH_NODE_DATA_BITS-1-CY_BVH_ELEMENT_COUNT_BITS)
#define _CY_BVH_ELEMENT_OFFSET_MASK	((1<<_CY_BVH_ELEMENT_OFFSET_BITS)-1)

#ifndef NUM_BINS
#define NUM_BINS 16
#endif

//-------------------------------------------------------------------------------

//! Bounding Volume Hierarchy class

class BVH
{
public:

	//!@name Constructor and destructor
	BVH() : nodes(0), elements(0) {}
	virtual ~BVH() { Clear(); }

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Node Access Methods
	//////////////////////////////////////////////////////////////////////////!//!//!

	//! Returns the index of the root node.
	unsigned int GetRootNodeID() const { return 1; }

	//! Returns the bounding box of the node as 6 float values.
	//! The first 3 values are the minimum x, y, and z coordinates and
	//! the last 3 values are the maximum x, y, and z coordinates of the box.
	const float* GetNodeBounds(unsigned int nodeID) const { return nodes[nodeID].GetBounds(); }

	//! Returns true if the node is a leaf node.
	bool IsLeafNode(unsigned int nodeID) const { return nodes[nodeID].IsLeafNode(); }

	//! Returns the index of the first child node (parent must be an internal node).
	unsigned int GetFirstChildNode(unsigned int parentNodeID) const { return nodes[parentNodeID].ChildIndex(); }

	//! Returns the index of the second child node (parent must be an internal node).
	unsigned int GetSecondChildNode(unsigned int parentNodeID) const { return nodes[parentNodeID].ChildIndex()+1; }

	//! Given the first child node index, returns the index of the second child node.
	unsigned int GetSiblingNode(unsigned int firstChildNodeID) const { return firstChildNodeID+1; }

	//! Returns the child nodes of the given node (parent must be an internal node).
	void GetChildNodes(unsigned int parent, unsigned int &child1, unsigned int &child2) const
	{
		child1 = GetFirstChildNode(parent);
		child2 = GetSiblingNode(child1);
	}

	//! Returns the number of elements inside the given node (must be a leaf node).
	unsigned int GetNodeElementCount(unsigned int nodeID) const  { return nodes[nodeID].ElementCount(); }

	//! Returns the list of element inside the given node (must be a leaf node).
	const unsigned int* GetNodeElements(unsigned int nodeID) const { return &elements[nodes[nodeID].ElementOffset()]; }

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Clear and Build Methods
	//////////////////////////////////////////////////////////////////////////!//!//!

	//! Clears the tree structure
	void Clear()
	{
		if (nodes) delete [] nodes;
		nodes = 0;
		if (elements) delete [] elements;
		elements = 0;
	}

	//! Builds the tree structure by recursively splitting the nodes. maxElementsPerNode cannot be larger than 8.
	void Build( unsigned int numElements, unsigned int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT )
	{
		Clear();
		if ( numElements == 0 ) return;
		if ( maxElementsPerNode > CY_BVH_MAX_ELEMENT_COUNT ) maxElementsPerNode = CY_BVH_MAX_ELEMENT_COUNT;
		elements = new unsigned int[numElements];
		for ( unsigned int i=0; i<numElements; i++ ) elements[i] = i;
		Box box;
		box.Init();
		for ( unsigned int i=0; i<numElements; i++ ) {
			Box b;
			GetElementBounds(i,b.b);
			box += b;
		}
		TempNode *tempRoot = new TempNode( numElements, 0, box );
		SplitTempNode(tempRoot,maxElementsPerNode);
		unsigned int numNodes = tempRoot->GetNumNodes();
		nodes = new Node[ numNodes+1 ];
		ConvertTempData( 1, tempRoot, 2 );
		delete tempRoot;
	}

	//////////////////////////////////////////////////////////////////////////!//!//!

protected:

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Methods to be implemented by sub-classes
	//////////////////////////////////////////////////////////////////////////!//!//!

	virtual void  GetElementBounds(unsigned int i, float box[6]) const=0;	//!< Sets box as the i^th element's bounding box.
	virtual float GetElementCenter(unsigned int i, int dimension) const=0;	//!< Returns the center of the i^th element in the given dimension

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Building method that can be overloaded
	//////////////////////////////////////////////////////////////////////////!//!//!

	//! Sorts the given elements of a temporary node while building the BVH hierarchy,
	//! such that first N elements are to be assigned to the first child and the 
	//! remaining elements are to be assigned to the second child node, then returns N.
	//! Returns zero, if the node is not to be split.
	//! The default implementation splits the temporary node down the middle of the
	//! widest axis of its bounding box.
	virtual unsigned int FindSplit(unsigned int elementCount, unsigned int *elements, const float *box, unsigned int maxElementsPerNode )
	{
		return MeanSplit(elementCount,elements,box,maxElementsPerNode);
	}

	//////////////////////////////////////////////////////////////////////////!//!//!

private:

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Internal storage
	//////////////////////////////////////////////////////////////////////////!//!//!

	struct Box
	{
		float b[6];
		Box() { Init(); }
		Box(const Box &box) { for(int i=0; i<6; i++) b[i]=box.b[i]; }
		void Init() { b[0]=b[1]=b[2]=1e30f; b[3]=b[4]=b[5]=-1e30f; }
		void operator += (const Box &box) { for(int i=0; i<3; i++) { if(b[i]>box.b[i])b[i]=box.b[i]; if(b[i+3]<box.b[i+3])b[i+3]=box.b[i+3]; } }
	};

	class Node
	{
	public:
		void SetLeafNode( const Box &bound, unsigned int elemCount, unsigned int elemOffset ) { box=bound; data=(elemOffset&_CY_BVH_ELEMENT_OFFSET_MASK)|((elemCount-1)<<_CY_BVH_ELEMENT_OFFSET_BITS)|_CY_BVH_LEAF_BIT_MASK; }
		void SetInternalNode( const Box &bound, unsigned int chilIndex ) { box=bound; data=(chilIndex&_CY_BVH_CHILD_INDEX_MASK); }
		unsigned int	ChildIndex()	const { return (data&_CY_BVH_CHILD_INDEX_MASK); }									//!< returns the index to the first child (must be internal node)
		unsigned int	ElementOffset()	const { return (data&_CY_BVH_ELEMENT_OFFSET_MASK); }									//!< returns the offset to the first element (must be leaf node)
		unsigned int	ElementCount()	const { return ((data>>_CY_BVH_ELEMENT_OFFSET_BITS)&_CY_BVH_ELEMENT_COUNT_MASK)+1; }	//!< returns the number of elements in this node (must be leaf node)
		bool			IsLeafNode()	const { return (data&_CY_BVH_LEAF_BIT_MASK)>0; }										//!< returns true if this is a leaf node
		const float*	GetBounds()		const { return box.b; }																//!< returns the bounding box of the node
	private:
		Box				box;	//!< bounding box of the node
		unsigned int	data;	//!< node data bits that keep the leaf node flag and the child node index or element count and element offset.
	};

	Node			*nodes;		//!< the tree structure that keeps all the node data (nodeData[0] is not used for cache coherency)
	unsigned int	*elements;	//!< indices of all elements in all nodes

	//////////////////////////////////////////////////////////////////////////!//!//!
	//@ Internal methods for building the BVH tree
	//////////////////////////////////////////////////////////////////////////!//!//!

	//! Temporary node class used for building the hierarchy and then converted to NodeData.
	class TempNode
	{
	public:
		TempNode( unsigned int count, unsigned int offset, const Box &boundBox) : child1(0), child2(0), elementCount(count), elementOffset(offset), box(boundBox) {}
		~TempNode() { if ( child1 ) delete child1; if ( child2 ) delete child2; }

		void Split( unsigned int child1ElementCount, const Box &child1Box, const Box &child2Box )
		{
			child1 = new TempNode(child1ElementCount,elementOffset,child1Box);
			child2 = new TempNode(ElementCount()-child1ElementCount,elementOffset+child1ElementCount,child2Box);
		}
		unsigned int GetNumNodes() const
		{
			unsigned int n = 1;
			if ( child1 ) n += child1->GetNumNodes();
			if ( child2 ) n += child2->GetNumNodes();
			return n;
		}
		bool IsLeafNode() const { return child1==0; }
		unsigned int ElementCount() const { return elementCount; }
		unsigned int ElementOffset() const { return elementOffset; }
		TempNode* GetChild1() { return child1; }
		TempNode* GetChild2() { return child2; }
		const Box& GetBounds() const { return box; }
	private:
		TempNode		*child1, *child2;
		Box				box;
		unsigned int	elementCount;
		unsigned int	elementOffset;
	};

	//! Recursively splits the given temporary node.
	void SplitTempNode(TempNode *tNode, unsigned int maxElementsPerNode)
	{
		const float *box = tNode->GetBounds().b;
		unsigned int *nodeElements = &elements[tNode->ElementOffset()];
		unsigned int child1ElemCount = FindSplit(tNode->ElementCount(),nodeElements,box,maxElementsPerNode);

		// If the FindSplit call does not return a valid split position
		if ( child1ElemCount == 0 || child1ElemCount >= tNode->ElementCount() ) {
			// if we must split anyway
			if ( tNode->ElementCount() > CY_BVH_MAX_ELEMENT_COUNT ) {
				// we split in half arbitrarily.
				child1ElemCount = tNode->ElementCount() / 2;
			} else {
				// otherwise, we reached a leaf node and no more split is necessary.
				return;
			}
		}

		// Compute child bounding boxes
		Box child1Box;
		Box child2Box;
		for ( unsigned int i=0; i<child1ElemCount; i++ ) {
			Box eBox;
			GetElementBounds( nodeElements[i], eBox.b );
			child1Box += eBox;
		}
		for ( unsigned int i=child1ElemCount; i<tNode->ElementCount(); i++ ) {
			Box eBox;
			GetElementBounds( nodeElements[i], eBox.b );
			child2Box += eBox;
		}

		// Split recursively
		tNode->Split( child1ElemCount, child1Box, child2Box );
		SplitTempNode(tNode->GetChild1(),maxElementsPerNode);
		SplitTempNode(tNode->GetChild2(),maxElementsPerNode);
	}

	//! Recursively converts the temporary node data to NodeData.
	unsigned int ConvertTempData( unsigned int nodeID, TempNode *tNode, unsigned int childIndex )
	{
		if ( tNode->IsLeafNode() ) {
			nodes[nodeID].SetLeafNode( tNode->GetBounds(), tNode->ElementCount(), tNode->ElementOffset() );
			return childIndex;
		} else {
			nodes[nodeID].SetInternalNode( tNode->GetBounds(), childIndex );
			unsigned int newChildIndex = ConvertTempData( childIndex, tNode->GetChild1(), childIndex+2 );
			return ConvertTempData( childIndex+1, tNode->GetChild2(), newChildIndex );
		}
	}

	//! Called by the default implementation of FindSplit.
	//! Splits the elements using the widest axis of the given bounding box.
	unsigned int MeanSplit(unsigned int elementCount, unsigned int *nodeElements, const float *box, unsigned int maxElementsPerNode )
	{
		if ( elementCount <= maxElementsPerNode ) return 0;
		float d[3] = { box[3]-box[0], box[4]-box[1], box[5]-box[2] };
		unsigned int sd[3]; // split dimensions
		sd[0] = d[0] >= d[1] ? ( d[0] >= d[2] ? 0 : 2 ) : ( d[1] >= d[2] ? 1 : 2 );	//find the widest axis of the bounding box
		sd[1] = (sd[0]+1) % 3;	//the other two axes
		sd[2] = (sd[0]+2) % 3;
		if ( d[sd[1]] < d[sd[2]] ) { int t=sd[1]; sd[1]=sd[2]; sd[2]=t; }

		unsigned int child1ElemCount = 0;
		for ( int s=0; s<3; s++ ) {	//for each axis, starting with the widest one
			unsigned int splitDim = sd[s];	//get the splitting axis
			float splitPos = 0.5f * ( box[splitDim] + box[splitDim+3] );	//find the midpoint of that axis
			unsigned int i=0, j=elementCount;
			//separate the triangles according to their relative position to the splitting plane
			while ( i<j ) {
				float center = GetElementCenter( nodeElements[i], splitDim );	//get the triangle's centroid
				if ( center <= splitPos ) {	//if it belongs on the left half
					i++;	//add one to the left side counter
				} else {
					j--;	//add one to the right side counter
					unsigned int t = nodeElements[i];	//put triangle i on the right side, grab one from the end of the array to take its place
					nodeElements[i] = nodeElements[j];
					nodeElements[j] = t;
				}
			}
			if ( i < elementCount && i > 0 ) {	//if this splitting axis actually did anything, we'll choose this axis to split on
				child1ElemCount = i;	//the left half will have this many triangles
				break;
			}
		}

		return child1ElemCount;
	}

	//////////////////////////////////////////////////////////////////////////!//!//!
};

//-------------------------------------------------------------------------------

#ifdef _CY_TRIMESH_H_INCLUDED_

//! Bounding Volume Hierarchy for triangular meshes (TriMesh)

class BVHTriMesh : public BVH
{
public:
	//!@name Constructors
	BVHTriMesh() : mesh(0) {}
	BVHTriMesh(const TriMesh *m) { SetMesh(m); }

	//! Sets the mesh pointer and builds the BVH structure.
	void SetMesh(const TriMesh *m, unsigned int maxElementsPerNode=CY_BVH_MAX_ELEMENT_COUNT)
	{
		mesh = m;
		Clear();
		Build(mesh->NF(),maxElementsPerNode);
	}
	
	//Find splits using the Surface Area Heuristic (SAH). Added by Justin Jensen
	virtual unsigned int FindSplit(unsigned int elementCount, unsigned int *elements, const float *box, unsigned int maxElementsPerNode) {
		//if there are fewer elements than the max
		if(elementCount <= maxElementsPerNode){
			//let's make a leaf!
			return 0;
		}

		//pick a number of possible splits, K = somewhere between 8 and 16
		unsigned int K = NUM_BINS;	//NUM_BINS is defined in main.cpp
		//create an array of K unsigned ints, for keeping track of the number of elements in each bin (set each to 0)
		unsigned int numInBin[NUM_BINS];
		//create an array of K bounding boxes, one for each bin
		float binBounds[NUM_BINS][6];
		for (int i = 0; i < NUM_BINS; i++) {
			binBounds[i][0] = BIGFLOAT;
			binBounds[i][1] = BIGFLOAT;
			binBounds[i][2] = BIGFLOAT;
			binBounds[i][3] = -BIGFLOAT;
			binBounds[i][4] = -BIGFLOAT;
			binBounds[i][5] = -BIGFLOAT;

			numInBin[i] = 0;
		}
		
		//calculate the dimensions of the bounding box
		float d[3] = { box[3] - box[0], box[4] - box[1], box[5] - box[2] };
		//find the widest axis of the bounding box, splitAxis
		unsigned int splitAxis = d[0] >= d[1] ? (d[0] >= d[2] ? 0 : 2) : (d[1] >= d[2] ? 1 : 2);

		//setup for assigning triangles to bins
		float epsilon = 0.001f;
		float k0 = box[splitAxis];
		float k1 = ((float)K * (1.0f - epsilon)) / (box[splitAxis + 3] - box[splitAxis]);
		float triBounds[6];
		//for each element
		for (int i = 0; i < elementCount; i++) {
			//Get the triangle's centroid
			float center = GetElementCenter(elements[i], splitAxis);
			//figure out which bin it belongs to. See http://www.sci.utah.edu/~wald/Publications/2007/ParallelBVHBuild/fastbuild.pdf, section 3.3
			unsigned int binID = (int)(k1 * (center - k0));
			//Get the triangle's bounds
			GetElementBounds(elements[i], triBounds);
			//add its bounds to the appropriate bin's bounding box
			binBounds[binID][0] = min(binBounds[binID][0], triBounds[0]);
			binBounds[binID][1] = min(binBounds[binID][1], triBounds[1]);
			binBounds[binID][2] = min(binBounds[binID][2], triBounds[2]);
			binBounds[binID][3] = max(binBounds[binID][3], triBounds[3]);
			binBounds[binID][4] = max(binBounds[binID][4], triBounds[4]);
			binBounds[binID][5] = max(binBounds[binID][5], triBounds[5]);
			//increment the count of the appropriate bin
			numInBin[binID]++;
		}
		
		//calculate the surface area of box, A_LR
		float A_LR = d[0] * d[1] * d[2];
		//calculate the cost of creating a leaf node (not splitting), C_l
		//C_l = A_LR * (float)elementCount
		float C_l = A_LR * (float)elementCount;
		
		//create an array of K-2 floats for keeping track of all possible left-hand box areas, leftBoxAreas
		float leftBoxAreas[NUM_BINS - 1];
		//create a BBox for a running total left-hand box, totalLeftBox
		float totalLeftBox[6] = { BIGFLOAT, BIGFLOAT, BIGFLOAT, -BIGFLOAT, -BIGFLOAT, -BIGFLOAT };
		//for each split from 0 to K-2, K_i
		for (int K_i = 0; K_i < K - 1; K_i++) {
			//get the bounding box for bin K_i
			//add it to totalLeftBox
			totalLeftBox[0] = min(totalLeftBox[0], binBounds[K_i][0]);
			totalLeftBox[1] = min(totalLeftBox[1], binBounds[K_i][1]);
			totalLeftBox[2] = min(totalLeftBox[2], binBounds[K_i][2]);
			totalLeftBox[3] = max(totalLeftBox[3], binBounds[K_i][3]);
			totalLeftBox[4] = max(totalLeftBox[4], binBounds[K_i][4]);
			totalLeftBox[5] = max(totalLeftBox[5], binBounds[K_i][5]);
			//calculate totalLeftBox's surface area, copy that to leftBoxAreas[K_i]
			leftBoxAreas[K_i] = (totalLeftBox[3] - totalLeftBox[0]) * (totalLeftBox[4] - totalLeftBox[1]) * (totalLeftBox[5] - totalLeftBox[2]);
		}

		//keep track of the index of the best split, low_K = -1
		int low_K = -1;
		//keep track of the best split cost so far, low_C = C_l
		float low_C = C_l;
		//keep track of the best number of triangles in the left node, low_N_L
		unsigned int low_N_L = 0;
		//create a BBox for a running total right-hand box, totalRightBox
		float totalRightBox[6] = { BIGFLOAT, BIGFLOAT, BIGFLOAT, -BIGFLOAT, -BIGFLOAT, -BIGFLOAT };
		//set N_R to 0
		unsigned int N_R = 0, N_L = 0;
		float A_R, A_L, C_s;
		//for each split from K-1 to 1, K_i
		for (int K_i = K - 1; K_i > 0; K_i--) {
			//get the bounding box for this bin, add it to totalRightBox
			totalRightBox[0] = min(totalRightBox[0], binBounds[K_i][0]);
			totalRightBox[1] = min(totalRightBox[1], binBounds[K_i][1]);
			totalRightBox[2] = min(totalRightBox[2], binBounds[K_i][2]);
			totalRightBox[3] = max(totalRightBox[3], binBounds[K_i][3]);
			totalRightBox[4] = max(totalRightBox[4], binBounds[K_i][4]);
			totalRightBox[5] = max(totalRightBox[5], binBounds[K_i][5]);
			//get its area, store in A_R
			A_R = (totalRightBox[3] - totalRightBox[0]) * (totalRightBox[4] - totalRightBox[1]) * (totalRightBox[5] - totalRightBox[2]);
			//get the number of triangles in this bin, add to N_R
			N_R += numInBin[K_i];
			//N_L = elementCount - N_R
			N_L = elementCount - N_R;
			//get the bounding box area for the left triangles from leftBoxAreas[K_i-1], store it in A_L
			A_L = leftBoxAreas[K_i - 1];
			//calculate this potential split's cost, C_s
			//C_s = A_L * (float)N_L + A_R * (float)N_R + A_LR
			C_s = A_L * (float)N_L + A_R * (float)N_R + A_LR;
			//if C_s is smaller than C_l
			if (C_s <= low_C) {
				//set low_K to K_i
				low_K = K_i;
				//set low_C to C_s
				low_C = C_s;
				//set low_N_L to N_L
				low_N_L = N_L;
			}
		}

		//if not splitting is better than splitting
		if (low_K == -1) {
			return 0;
		}
		
		//actually separate the triangle indices
		//return the number of triangles in the left half
		unsigned int i = 0, j = elementCount;
		//separate the triangles according to their relative position to the splitting plane
		while (i<j) {
			float center = GetElementCenter(elements[i], splitAxis);	//get the triangle's centroid
			if ((int)(k1 * (center - k0)) < low_K) {	//if it belongs on the left half
				i++;	//add one to the left side counter
			}
			else {
				j--;	//add one to the right side counter
				unsigned int t = elements[i];	//put triangle i on the right side, grab one from the end of the array to take its place
				elements[i] = elements[j];
				elements[j] = t;
			}
		}
		//if (i < elementCount && i > 0) {	//if this splitting axis actually did anything, we'll choose this axis to split on
		//	child1ElemCount = i;	//the left half will have this many triangles
		//	break;
		//}

		return low_N_L;
	}/**/

protected:
	//! Sets box as the i^th element's bounding box.
	virtual void GetElementBounds(unsigned int i, float box[6]) const
	{
		const TriMesh::TriFace &f = mesh->F(i);
		cyPoint3f p = mesh->V( f.v[0] );
		box[0]=box[3]=p.x; box[1]=box[4]=p.y; box[2]=box[5]=p.z;
		for ( int j=1; j<3; j++ ) { // for each triangle
			cyPoint3f p = mesh->V( f.v[j] );
			for ( int k=0; k<3; k++ ) { // for each dimension
				if ( box[k] > p[k] ) box[k] = p[k];
				if ( box[k+3] < p[k] ) box[k+3] = p[k];
			}
		}
	}

	//! Returns the center of the i^th element in the given dimension.
	virtual float GetElementCenter(unsigned int i, int dim) const
	{
		const TriMesh::TriFace &f = mesh->F(i);
		return ( mesh->V(f.v[0])[dim] + mesh->V(f.v[1])[dim] + mesh->V(f.v[2])[dim] ) / 3.0f;
	}

private:
	const TriMesh *mesh;
};

#endif

//-------------------------------------------------------------------------------
} // namespace cy
//-------------------------------------------------------------------------------

typedef cy::BVH cyBVH;	//!< Bounding Volume Hierarchy class

#ifdef _CY_TRIMESH_H_INCLUDED_
typedef cy::BVHTriMesh cyBVHTriMesh;	//!< BVH hierarchy for triangular meshes (TriMesh)
#endif

//-------------------------------------------------------------------------------

#endif

