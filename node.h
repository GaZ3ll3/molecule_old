/*
 * node.h
 *
 *  Created on: Oct 9, 2016
 *  Updated on: Oct 20, 2016
 *      Author: Yimin Zhong
 */

#ifndef FMM_NODE_H_H
#define FMM_NODE_H_H

#if !defined __extern_always_inline && defined __clang__
# if defined __GNUC_STDC_INLINE__ || defined __GNUC_GNU_INLINE__
#  define __extern_inline extern __inline __attribute__ ((__gnu_inline__))
#  define __extern_always_inline \
  extern __always_inline __attribute__ ((__gnu_inline__))
# else
#  define __extern_inline extern __inline
#  define __extern_always_inline extern __always_inline
# endif
#endif


#include <unordered_set>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

#define __eps 1e-12

using std::unordered_set;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class basePoint {
public:
    double x;
    double y;
    double z;

    basePoint() : x(0.), y(0.), z(0.) {}

    basePoint(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    virtual ~basePoint() {}

    bool operator>=(const basePoint &a) {
        return ( x >= a.x ) && ( y >= a.y ) && (z >= a.z);
    }

    bool operator<=(const basePoint &a) {
        return ( x <= a.x ) && ( y <= a.y ) && (z <= a.z);
    }

    bool operator==(const basePoint &a) {
        return (fabs(x - a.x) < __eps) && (fabs(y - a.y) < __eps) && (fabs(z - a.z) < __eps);
    }
};


class point : public basePoint {
public:
    int triangleId;

    point() : basePoint() { triangleId = -1; }

    point(double _x, double _y, double _z) : basePoint(_x, _y, _z) { triangleId = -1; }

    point(double _x, double _y, double _z, int _id) : basePoint(_x, _y, _z) { triangleId = _id; }

    ~point() {}
};

class baseNode {
public:
    int parent;
    int child[8];

    unordered_set<int> uList;
    unordered_set<int> vList;
    unordered_set<int> wList;
    unordered_set<int> xList;

    int nLevel;

    int nUList;
    int nVList;
    int nWList;
    int nXList;

    int nodeIndex;
    point center;
    point radius;

    int nSource;
    int nTarget;

    vector<int> sourceIndex;
    vector<int> targetIndex;

    bool isLeaf;
    bool isEmpty;
    bool chargeComputed;

    baseNode(int level, int index) {
        parent = -1;
        for (int i = 0; i < 8; ++i) {
            child[i] = -1;
        }
        nLevel = level;
        nodeIndex = index;
        isLeaf = false;
        isEmpty = false;
        chargeComputed = false;
        nSource = 0;
        nTarget = 0;
        nUList = 0;
        nVList = 0;
        nWList = 0;
        nXList = 0;

    }

    virtual ~baseNode() {}

} ;


class node : public baseNode {
public:
    node(int level, int index) : baseNode(level, index) {

    }
    ~node() {}

    /*
     *  some internal members
     */
    vector<point> scaledCnode;
    MatrixXd potential;
    MatrixXd nodeCharge;
    MatrixXd charge;
    MatrixXd nodePotential;
    MatrixXd R;
    MatrixXd L;


};


#endif //FMM_NODE_H_H
