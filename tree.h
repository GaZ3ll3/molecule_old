//
// Created by lurker on 9/29/16.
//

#ifndef FMM_TREE_H
#define FMM_TREE_H

#include "node.h"
#include "measure.h"
#include <chrono>
//#include "omp.h"

class tree : public measure {
public:
    vector<node> dict;
    int maxId;
    int root;
    int nSource;
    int nTarget;

    int rank;
    int maxLevel;

    vector<point> sourceTree;
    vector<point> targetTree;

    point center;
    point radius;

    tree() {
        maxId = -1;
        root = -1;
    }
    ~tree() {

    }

    void populate(vector<point>& _source, vector<point>& _target, int _nSource, int _nTarget, int _rank, int _maxLevel);
    void output(std::string file);

protected:
    void getCenterRadius(vector<point>& _source);
    void assignChildren(int _id, int _maxLevel);
    void buildTree();
    void buildNode(int _id, point& min_p, point& max_p);
    int findNode(int _id, point& p);
    bool isAdjacent(int _aId, int _bId);

};

#endif //FMM_TREE_H
