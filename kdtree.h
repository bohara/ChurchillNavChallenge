#ifndef KDTREE_H
#define KDTREE_H

#include <QDebug>

#include <stddef.h>
#include "point.h"
#include "boundedpqueue.h"

#pragma pack(push, 1)
struct Point
{
    float x;
    float y;
};
#pragma pack(pop)


/***/
template <size_t N>
class KDTree
{
    /* Tree Node structure*/
    struct KDNode
    {
        int sortDim;
        PointV<N> pos;
        KDNode *left;
        KDNode *right;

        KDNode() : sortDim(0),
            left(nullptr),
            right(nullptr) {}
    };

    KDNode *kdRoot;
    size_t dims;
    size_t nodeCount;

    float bestDistance;
    PointV<N> bestGuess;

    /// Not essential!!
    bool branchDirection(const KDNode&_node, const PointV<N>&_p) const;

    KDNode* copyTree(const KDNode *_rootNode)
    {
        if(_rootNode == nullptr)
            return nullptr;

        KDNode* _node = new KDNode();
        _node->pos = _rootNode->pos;
        _node->sortDim = _rootNode->sortDim;

        // Copy left and right subtree
        _node->left = copyTree(_rootNode->left);
        _node->right = copyTree(_rootNode->right);

        return _node;
    }

    void deleteNode(const KDNode *_node)
    {
        if(_node == nullptr) // return when you hit the Child of leaf
            return;

        deleteNode(_node->left);
        deleteNode(_node->right);
        delete _node;
    }

    KDNode* findNode(const PointV<N>&pt) const
    {
        KDNode *head = kdRoot;

        //iterate through left or right child depending upon the value
        while(head != nullptr)
        {
            if(pt == head->pos)
                return head;

            head = (pt[head->sortDim] < head->pos[head->sortDim]) ? head->left :
                                                                    head->right;
        }
        return nullptr;
    }

    KDNode* buildBranch(KDNode **_head, const int&_sortDim,
                        typename std::vector<PointV<N> >::iterator _begin,
                        typename std::vector<PointV<N> >::iterator _end);

    bool checkInRange(const PointV<N>&pt, const PointV<N>&_min, const PointV<N>&_max)
    {
        for(size_t i = 0; i < N; ++i)
        {
            if(pt[i] < _min[i] || pt[i] > _max[i])
                return false;
        }
        return true;
    }

public:
    KDTree();
    ~KDTree();

    // copy constructor
    KDTree(const KDTree& rvTree);

    // construtor given all the points to be inserted in tree
    KDTree(typename std::vector<PointV<N> >::iterator begin,
           typename std::vector<PointV<N> >::iterator end);

    // assignment operator
    KDTree& operator= (const KDTree& rvTree);

    size_t dimension() const;
    bool empty() const;
    size_t size() const;

    void insert(const PointV<N>&pt, const int&value);
    //sbfind(const Point&pt);

    int& at(const PointV<N>&pt);
    const int& at(const PointV<N>&pt) const;
    bool contains(const PointV<N>&pt) const;
    int& operator[] (const PointV<N>&pt);

    //Finds k nearest points and returns the nearest point!!
    PointV<N> kNNValue(const PointV<N>&pt, const size_t& k);
    void findNearestNeighbors(const PointV<N>&pt, KDNode*curNode, BoundedPQueue<PointV<N>> &nnQueue);

    void searchRange(const PointV<N>&min, const PointV<N>&max);
    void findElementsInBox(const KDNode*head, const PointV<N>&min, const PointV<N>&max,
                           std::vector<PointV<N>>&nElements);
};

template <size_t N>
KDTree<N>::KDTree()
{
    dims = N;
    nodeCount = 0;
    kdRoot = nullptr;
}

template <size_t N>
KDTree<N>::KDTree(typename std::vector<PointV<N>>::iterator begin,
                  typename std::vector<PointV<N>>::iterator end)
{
    // Initialize Kd-Tree members
    dims = N;
    nodeCount = 0;
    kdRoot = nullptr;

    // Insert the sorted array of points into the kdtree
    size_t sz = end - begin;
    size_t median =  sz / 2;    //(sz % 2 == 0) ? (sz/2 + 1) : sz/2 + 1;

    // Insert root node
    insert(*(static_cast<typename std::vector<PointV<N>>::iterator>(begin+median)), 0);

    // Divide array into two halves from the median
    if(kdRoot != nullptr)
    {
        // Insert left branch
        buildBranch(&kdRoot->left, kdRoot->sortDim, begin, begin + median);
        // Insert right branch
        buildBranch(&kdRoot->right, kdRoot->sortDim, begin + median + 1, end);
    }
    //qDebug() << sz;
}

template <size_t N>
typename KDTree<N>::KDNode* KDTree<N>::buildBranch(KDNode **_head, const int&_sortDim,
                                                   typename std::vector<PointV<N> >::iterator _begin,
                                                   typename std::vector<PointV<N> >::iterator _end)
{
    // Insert the sorted array of points into the kdtree
    size_t sz = _end - _begin;
    size_t median = sz / 2;   //(sz % 2 == 0) ? (sz/2 + 1) : sz/2 + 1;

    if(sz <= 0 || median > sz)
        return nullptr;

    // Sort the partial array w.r.t corresponding index
    int sortIdx = (_sortDim + 1) % dimension();
    std::sort(_begin, _end, [ sortIdx ](const PointV<N>&a, const PointV<N>&b) -> bool
    { return (a[sortIdx] < b[sortIdx]);   }  );

    // Insert root node
    //insert(*(static_cast<typename std::vector<PointV<N>>::iterator>(begin+median)), 0);
    typename std::vector<PointV<N> >::iterator itr = (static_cast<typename std::vector<PointV<N>>
                                                      ::iterator>(_begin+median));
    if((*_head) == nullptr)
    {
        *_head = new KDNode();
        (*_head)->sortDim = (_sortDim + 1) % dimension();
        (*_head)->left = nullptr;
        (*_head)->right = nullptr;
        (*_head)->pos = (*itr);

        nodeCount++;

        // Insert left branch
        buildBranch(&(*_head)->left, (*_head)->sortDim, _begin, _begin + median);
        // Insert right branch
        buildBranch(&(*_head)->right, (*_head)->sortDim, _begin + median + 1, _end);

        return *_head;
    }
    else
    {
#if 0
        insert(*(static_cast<typename std::vector<PointV<N>>::iterator>(_begin+median)), 0);
#endif
        throw(std::range_error("Should not have reached here!!"));
        return nullptr;
    }
}

template <size_t N>
KDTree<N>::~KDTree()
{
    deleteNode(kdRoot);
}

template <size_t N>
KDTree<N>::KDTree(const KDTree &rvTree)
{
    kdRoot = copyTree(rvTree.kdRoot);    // Copy the tree starting from root

    nodeCount = rvTree.size();
    dims = rvTree.dimension();
}

template <size_t N>
KDTree<N>& KDTree<N>::operator =(const KDTree& rvTree)
{
    if(this != rvTree)
    {
        deleteNode(kdRoot);
        kdRoot = copyTree(rvTree.kdRoot);    // Copy the tree starting from root

        nodeCount = rvTree.size();
        dims = rvTree.dimension();
    }
    return *this;
}

// Assuming there can only be maximum dimentions of 3
template <size_t N>
bool KDTree<N>::branchDirection(const KDNode&_node, const PointV<N>&_p) const
{
    int sortIndex = _node.sortDim;

    if(_p[sortIndex] < _node.pos[sortIndex])
        return true;
    else
        return false;
}

template <size_t N>
size_t KDTree<N>::dimension() const
{
    return N;
}

template <size_t N>
bool KDTree<N>::empty() const
{
    if(kdRoot == nullptr || size() == 0 )
        return true;
    else
        return false;
}

template <size_t N>
void KDTree<N>::insert(const PointV<N> &pt, const int &value)
{
    if(kdRoot == nullptr) // Insert at node, else traverse the child
    {
        kdRoot = new KDNode();
        kdRoot->sortDim = 0;
        kdRoot->left = nullptr;
        kdRoot->right = nullptr;
        kdRoot->pos = pt;

        nodeCount++;
        return;
    }

    KDNode *head = kdRoot;

    //iterate through left or right child depending upon the value
    while(head != nullptr)
    {
        // Check if the incoming point is equal to the head
        // If so replace the associated value
        //if(pt == head->pos)
        //    head->rank = value;

        if(pt[head->sortDim] < head->pos[head->sortDim])      //insert on left child
        {
            if(head->left != nullptr)
            {
                head = head->left;
                continue;
            }
            else
            {
                head->left = new KDNode();
                head->left->pos = pt;
                head->left->sortDim = (head->sortDim + 1) % dims;
                nodeCount++;
                break;  // head = head->left;
            }
        }
        else                              //insert on right child
        {
            if(head->right != nullptr)
            {
                head = head->right;
                continue;
            }
            else
            {
                head->right = new KDNode();
                head->right->pos = pt;
                head->right->sortDim = (head->sortDim + 1) % dims;
                nodeCount++;
                break;  // head = head->right;
            }
        }
    }
}

template <size_t N>
size_t KDTree<N>::size() const
{
    return nodeCount;
}

template <size_t N>
bool KDTree<N>::contains(const PointV<N> &pt) const
{
    return (nullptr != findNode(pt));
}

template <size_t N>
int& KDTree<N>::at(const PointV<N> &pt)
{
    KDNode *target = findNode(pt);
    if(nullptr == target)
        throw std::out_of_range("Point not found in kd-tree.");
    else
        return target->sortDim; // should return associated value
}

template <size_t N>
const int& KDTree<N>::at(const PointV<N> &pt) const
{
    KDNode *target = findNode(pt);
    if(nullptr == target)
        throw std::out_of_range("Point not found in kd-tree.");
    else
        return target->sortDim; // should return associated value
}

template <size_t N>
int& KDTree<N>::operator [](const PointV<N>&pt)
{
    KDNode *target = findNode(pt);
    if(nullptr == target)
        insert(pt, -1);
    else
        return target->sortDim;
}

template <size_t N>
PointV<N> KDTree<N>::kNNValue(const PointV<N> &pt, const size_t &k)
{
    BoundedPQueue<PointV<N> > neighbors(k);      // NN bounded priority queue
    findNearestNeighbors(pt, kdRoot, neighbors);       // test point and root node

    // Print all the kNN points
#ifdef QT_DEBUG
    while(!neighbors.empty())
    {
        PointV<N> _point = neighbors.dequeueMin();
        qDebug() << _point[0] << _point[1];
    }
#endif
}

template <size_t N>
void KDTree<N>::findNearestNeighbors(const PointV<N> &pt, KDNode *curNode,
                                     BoundedPQueue<PointV<N> > &nnQueue)
{
    if(curNode == nullptr)
        return;

    // Insert in the priotity queue
    float distance = Distance(pt, curNode->pos);
    nnQueue.enqueue(curNode->pos, (double)distance);

    if(pt[curNode->sortDim] < curNode->pos[curNode->sortDim])   // left subtree
    {
        findNearestNeighbors(pt, curNode->left, nnQueue);

        // If hypersphere intersect the place, search other side of plane
        if(nnQueue.size() < nnQueue.maxSize() ||
                fabs(double(pt[curNode->sortDim] - curNode->pos[curNode->sortDim])) < nnQueue.worst())
        {   findNearestNeighbors(pt, curNode->right, nnQueue);  }
    }
    else
    {
        findNearestNeighbors(pt, curNode->right, nnQueue);      //righ subtree

        // If hypersphere intersect the place, search other side of plane
        if(nnQueue.size() < nnQueue.maxSize() ||
                fabs(double(pt[curNode->sortDim] - curNode->pos[curNode->sortDim])) < nnQueue.worst())
        {   findNearestNeighbors(pt, curNode->left, nnQueue);   }
    }

}

template <size_t N>
void KDTree<N>::searchRange(const PointV<N> &min, const PointV<N> &max)
{
    std::vector<PointV<N> > pointsInBox;
    findElementsInBox(kdRoot, min, max, pointsInBox);

    // Print all the elements in the list
    // Print all the kNN points
#ifdef QT_DEBUG
    size_t size = pointsInBox.size();
    qDebug() << size;
    for(size_t i = 0; i<size ; ++i)
    {
        PointV<N> _point = pointsInBox[i];
        qDebug() << _point[0] << _point[1];
    }
#endif
}

template <size_t N>
void KDTree<N>::findElementsInBox(const KDNode*head, const PointV<N>&min,
                                  const PointV<N>&max, std::vector<PointV<N>>&nElements)
{
    if(head == nullptr)
        return;

    int sortIndex = head->sortDim;

    // Check the overlapping of Rectangle with current node
    if(head->pos[sortIndex] > min[sortIndex] && head->pos[sortIndex] < max[sortIndex])
    {
        // Node overlaps the search range(rectangle)
        findElementsInBox(head->left, min, max, nElements);
        findElementsInBox(head->right, min, max, nElements);
    }
    else if(head->pos[sortIndex] > min[sortIndex] && head->pos[sortIndex] > max[sortIndex])
    {
        // Range is left of current node
        findElementsInBox(head->left, min, max, nElements);
    }
    else if(head->pos[sortIndex] <= min[sortIndex] && head->pos[sortIndex] <= max[sortIndex])
    {
        // Range is right of the current node
        findElementsInBox(head->right, min, max, nElements);
    }

    // Check if point at current node is in the range, if so add to the list
    if(checkInRange(head->pos, min, max))
        nElements.push_back(head->pos);
}

#endif // KDTREE_H
