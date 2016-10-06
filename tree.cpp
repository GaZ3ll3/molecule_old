//
// Created by lurker on 9/29/16.
//

#include "tree.h"

#include <fstream>
#include <cassert>
#include <queue>



void tree::populate(vector<point> &_source, vector<point> &_target, int _nSource, int _nTarget, int _rank,
                    int _maxLevel) {
    this->sourceTree = _source;
    this->targetTree = _target;
    this->nSource = _nSource;
    this->nTarget = _nTarget;
    this->maxLevel = 0;
    this->rank = _rank;
    getCenterRadius(_source);
    this->root = 0;
    this->dict.push_back(node(0, 0));
    this->maxId = root;
    dict[root].nSource = nSource;
    dict[root].nTarget = nTarget;
    dict[root].center = center;
    dict[root].radius = radius;
    dict[root].sourceIndex.resize((unsigned long) nSource);
    dict[root].targetIndex.resize((unsigned long) nTarget);
    for (int i = 0; i < nSource; ++i) {
        dict[root].sourceIndex[i] = i;
    }
    for (int i = 0; i < nTarget; ++i) {
        dict[root].targetIndex[i] = i;
    }


    RUN("initialization", assignChildren(root, _maxLevel));
    RUN("assign lists"  , buildTree());


}

void tree::getCenterRadius(vector<point> &_source) {
    assert(_source.size() > 0);
    double x_max = _source[0].x;
    double x_min = _source[0].x;
    double y_max = _source[0].y;
    double y_min = _source[0].y;
    double z_max = _source[0].z;
    double z_min = _source[0].z;
    for (int i = 0; i < _source.size(); ++i) {
        x_max = std::max(x_max, _source[i].x);
        y_max = std::max(y_max, _source[i].y);
        z_max = std::max(z_max, _source[i].z);
        x_min = std::min(x_min, _source[i].x);
        y_min = std::min(y_min, _source[i].y);
        z_min = std::min(z_min, _source[i].z);
    }
    this->center.x = (x_max + x_min)/2.0;
    this->center.y = (y_max + y_min)/2.0;
    this->center.z = (z_max + z_min)/2.0;
    this->radius.x = (x_max - x_min)/2.0;
    this->radius.y = (y_max - y_min)/2.0;
    this->radius.z = (z_max - z_min)/2.0;
}

void tree::assignChildren(int _id, int _maxLevel) {
    /*
     * when assigning children nodes, the points are not assigned due to storage.
     *
     * Now the limitation of nodes is around 2^24.
     */
    assert(dict.size() > _id); // check if this node exists
    assert(root != -1); // check tree is non-empty

    // check source
    if (dict[_id].nSource == 0) {
        dict[_id].isLeaf = true;
        dict[_id].isEmpty = true;
    }
    else {
        // divide
        if ((dict[_id].nSource <= rank) || (dict[_id].nLevel == _maxLevel)) {
            dict[_id].isLeaf = true;
            if (maxLevel < dict[_id].nLevel) {
                maxLevel = dict[_id].nLevel;
            }
        }
        else {
            // not a leaf
            for (int i = 0; i < 8; ++i) {
                maxId += 1;
                dict[_id].child[i] = maxId;
                dict.push_back(node(dict[_id].nLevel + 1, i));
                dict[maxId].parent = _id;
                dict[maxId].center.x = dict[_id].center.x + ((i & 1) - 0.5) * dict[_id].radius.x;
                dict[maxId].center.y = dict[_id].center.y + (((i >> 1) & 1) - 0.5) * dict[_id].radius.y;
                dict[maxId].center.z = dict[_id].center.z + ((i >> 2) - 0.5) * dict[_id].radius.z;
                dict[maxId].radius.x = dict[_id].radius.x * 0.5;
                dict[maxId].radius.y = dict[_id].radius.y * 0.5;
                dict[maxId].radius.z = dict[_id].radius.z * 0.5;
                dict[maxId].nSource = 0;
                dict[maxId].nTarget = 0;
            }

            /*
             * can be accelerated by **reduce**
             */
            for (int i = 0; i < dict[_id].nSource; ++i) {
                int index = dict[_id].sourceIndex[i];
                int z_bit = sourceTree[index].z < dict[_id].center.z ? 0:1;
                int y_bit = sourceTree[index].y < dict[_id].center.y ? 0:1;
                int x_bit = sourceTree[index].x < dict[_id].center.x ? 0:1;
                int childIndex = 4 * z_bit + 2 * y_bit + x_bit;

                int childId = dict[_id].child[childIndex];
                dict[childId].sourceIndex.push_back(index);
                dict[childId].nSource += 1;
            }

            /*
             * can be accelerated by **reduce**
             */
            for (int i = 0; i < dict[_id].nTarget; ++i) {
                int index = dict[_id].targetIndex[i];
                int z_bit = targetTree[index].z < dict[_id].center.z ? 0:1;
                int y_bit = targetTree[index].y < dict[_id].center.y ? 0:1;
                int x_bit = targetTree[index].x < dict[_id].center.x ? 0:1;
                int childIndex = 4 * z_bit + 2 * y_bit + x_bit;

                int childId = dict[_id].child[childIndex];
                dict[childId].targetIndex.push_back(index);
                dict[childId].nTarget += 1;
            }

            for (int i = 0; i < 8; ++i) {
                assignChildren(dict[_id].child[i], _maxLevel);
            }
        }
    }
}

void tree::buildTree() {
    omp_set_num_threads(4);
    point min_p(dict[root].center.x - dict[root].radius.x,
                dict[root].center.y - dict[root].radius.y,
                dict[root].center.z - dict[root].radius.z);
    point max_p(dict[root].center.x + dict[root].radius.x,
                dict[root].center.y + dict[root].radius.y,
                dict[root].center.z + dict[root].radius.z);
    int i;
#pragma omp parallel for private(i) shared(min_p, max_p) schedule(dynamic)
    for (i = 0; i < dict.size(); ++i) {
        buildNode(i, min_p, max_p);
    }
}

void tree::buildNode(int _id, point &min_p, point &max_p) {
    node& n = dict[_id];
    n.uList.clear();
    n.vList.clear();
    n.wList.clear();
    n.xList.clear();

    // not root
    if (n.parent != -1) {
        node pn = dict[n.parent];
        double dx = n.radius.x;
        double dy = n.radius.y;
        double dz = n.radius.z;
        double xs = pn.center.x - dx;
        double ys = pn.center.y - dy;
        double zs = pn.center.z - dz;

        point cur;

        for (int x_id = -2; x_id < 4; x_id++) {
            for (int y_id = -2; y_id < 4; y_id++) {
                for (int z_id = -2; z_id < 4; z_id++) {
                    cur.x = xs + 2 * x_id * dx;
                    cur.y = ys + 2 * y_id * dy;
                    cur.z = zs + 2 * z_id * dz;

                    // check box and not itself.
                    if (cur <= max_p && cur >= min_p && !(cur == n.center)) {
                        //find node.
                        int curId = findNode(0, cur);
                        bool adj = isAdjacent(_id, curId);
                        node& curNode = dict[curId];

                        if (curNode.nLevel < n.nLevel) {
                            if (adj) {
                                if (curNode.isLeaf) {
                                    n.uList.insert(curId);
                                }
                            }
                            else {
                                n.xList.insert(curId);
                            }
                        }

                        if (curNode.nLevel == n.nLevel) {
                            if (!adj) {
                                n.vList.insert(curId);
                            }
                            else {
                                if (n.isLeaf) {
                                    std::queue<int> rest;
                                    rest.push(curId);
                                    while (!rest.empty()) {
                                        int frontId = rest.front(); rest.pop();
                                        node& frontNode = dict[frontId];
                                        if (!isAdjacent(frontId, _id)) {
                                            n.wList.insert(frontId);
                                        }
                                        else {
                                            if (frontNode.isLeaf) {
                                                n.uList.insert(frontId);
                                            }
                                            else {
                                                for (int i = 0; i < 8; ++i) {
                                                    rest.push(frontNode.child[i]);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (n.isLeaf) {
        n.uList.insert(_id);
    }

    n.nUList = (int) n.uList.size();
    n.nWList = (int) n.wList.size();
    n.nVList = (int) n.vList.size();
    n.nXList = (int) n.xList.size();
}


int tree::findNode(int _id, point &p) {
    node& n = dict[_id];
    if (n.center == p) return _id;
    else {
        if (n.isLeaf) {
            return _id;
        }
        else {
            int x_bit = n.center.x > p.x ? 0 : 1;
            int y_bit = n.center.y > p.y ? 0 : 1;
            int z_bit = n.center.z > p.z ? 0 : 1;
            int id = 4 * z_bit + 2 * y_bit + x_bit;
            return findNode(n.child[id], p);
        }
    }
}

bool tree::isAdjacent(int _aId, int _bId) {
    node& nA = dict[_aId];
    node& nB = dict[_bId];
    double diff_x = fabs(nA.center.x - nB.center.x);
    double diff_y = fabs(nA.center.y - nB.center.y);
    double diff_z = fabs(nA.center.z - nB.center.z);
    double r_x = fabs(nA.radius.x + nB.radius.x);
    double r_y = fabs(nA.radius.y + nB.radius.y);
    double r_z = fabs(nA.radius.z + nB.radius.z);

    bool rdx = r_x >= diff_x - __eps;
    bool rdy = r_y >= diff_y - __eps;
    bool rdz = r_z >= diff_z - __eps;

    bool x_adj = (fabs(diff_x - r_x) < __eps) && (rdy && rdz);
    bool y_adj = (fabs(diff_y - r_y) < __eps) && (rdx && rdz);
    bool z_adj = (fabs(diff_z - r_z) < __eps) && (rdy && rdx);

    return x_adj || y_adj || z_adj;

}


void tree::output(std::string file) {
    std::ofstream file_stream(file);
    if (file_stream.is_open()) {

        for (int i = 0; i < dict.size(); ++i) {
            file_stream << dict[i].center.x << " "
                        << dict[i].center.y << " "
                        << dict[i].center.z << " "
                        << dict[i].radius.x << " "
                        << dict[i].radius.y << " "
                        << dict[i].radius.z << " "
                        << dict[i].nVList << " " << dict[i].nXList << " " << dict[i].nUList <<" "<< dict[i].nWList <<"\n";
        }

        file_stream.close();
    }
    else {
        std::cout << "cannot open file: " << file << std::endl;
    }
}