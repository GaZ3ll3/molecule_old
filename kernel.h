//
// Created by lurker on 9/29/16.
//

#ifndef FMM_KERNEL_H
#define FMM_KERNEL_H


#include "tree.h"
#include <functional>
#include <Eigen/Dense>

using namespace Eigen;

class kernel : public measure {
public:
    tree t;
    MatrixXd chargeTree;
    std::function<double(point&, point&)> eval;
    int rank;

    MatrixXd R[8];

    int nChebyshev;
    VectorXd chebyNode;
    MatrixXd tNode;

    kernel() {}
    ~kernel() {}

    void initialize(int _nChebyshev, vector<point>& _source, vector<point>& _target, double* const _charge, int _nSource,
                    int _nTarget, int _rank, int _maxLevel) {
        // populate the kd-tree.
        t.populate(_source, _target, _nSource, _nTarget, _rank, _maxLevel);

        nChebyshev = _nChebyshev;

        chargeTree = Map<MatrixXd>(_charge, _nSource, 1);

        // nChebyshev^3 nodes are used for interpolation.
        rank = nChebyshev * nChebyshev * nChebyshev;

        chebyNode = VectorXd::Zero(nChebyshev);
        getStandardChebyNodes(nChebyshev, chebyNode);

        tNode = MatrixXd::Zero(nChebyshev, nChebyshev);
        getStandardChebyPoly(nChebyshev, nChebyshev, chebyNode, tNode);

        getTransfer(nChebyshev, chebyNode, tNode, R);
    }


    void getStandardChebyNodes(int _nChebyshev, VectorXd& _chebyNode) {
        _chebyNode = VectorXd::Zero(_nChebyshev);
        for (int i = 0; i < _nChebyshev; ++i) {
            _chebyNode(i) = -cos((i + 0.5) * M_PI/_nChebyshev);
        }
    }

    void getStandardChebyPoly(int _nChebyPoly, int _N, VectorXd& _x, MatrixXd& _T) {
        _T = MatrixXd::Zero(_N, _nChebyPoly);
        _T.col(0) = VectorXd::Ones(_N);
        if (_nChebyPoly > 1) {
            _T.col(1) = _x;
            for (int i = 2; i < _nChebyPoly; ++i) {
                _T.col(i) = 2.0 * _x.cwiseProduct(_T.col(i - 1)) - _T.col(i - 2);
            }
        }
    }

    void getTransferFromParentChebyshevToChildrenChebyshev(int _nChebyshev, VectorXd& _chebyNode, MatrixXd& _tNode, MatrixXd& _transfer) {
        VectorXd childChebyNode (2 * _nChebyshev);
        childChebyNode.segment(0, nChebyshev)          = 0.5 *(_chebyNode - VectorXd::Ones(_nChebyshev));
        childChebyNode.segment(nChebyshev, nChebyshev) = 0.5 *(_chebyNode + VectorXd::Ones(_nChebyshev));
        getStandardChebyPoly(_nChebyshev, 2 * _nChebyshev, childChebyNode, _transfer);
        _transfer = (2.0 * _transfer * _tNode.transpose() - MatrixXd::Ones(2 * _nChebyshev, _nChebyshev))/_nChebyshev;
    }

    void getTransfer(int _nChebyshev, VectorXd& _chebyNode, MatrixXd& _tNode, MatrixXd* R) {
        MatrixXd S;
        S = MatrixXd::Zero(2 * _nChebyshev, _nChebyshev);
        getTransferFromParentChebyshevToChildrenChebyshev(_nChebyshev, _chebyNode, _tNode, S);
        MatrixXd Transfer[2];
        Transfer[0] = S.block(0, 0, _nChebyshev, _nChebyshev);
        Transfer[1] = S.block(_nChebyshev, 0, _nChebyshev, _nChebyshev);
        int _rank = _nChebyshev * _nChebyshev * _nChebyshev;
        for (int i = 0; i < 8; ++i) {
            R[i] = MatrixXd::Zero(_rank, _rank);
        }

        // follow bit representaion.
        for (int i = 0; i < _nChebyshev; ++i) {
            for (int j =0; j < _nChebyshev; ++j) {
                for (int k = 0; k < _nChebyshev; ++k) {
                    for (int l = 0; l < _nChebyshev; ++l) {
                        for (int m = 0; m < _nChebyshev; ++m) {
                            for (int n = 0; n < _nChebyshev; ++n) {
                                for (int id = 0; id < 8; ++id) {
                                    int bit[3];
                                    bit[0] = (id >> 0) & 1;
                                    bit[1] = (id >> 1) & 1;
                                    bit[2] = (id >> 2) & 1;
                                    R[id](i * _nChebyshev * _nChebyshev + j * _nChebyshev + k,
                                         l * _nChebyshev * _nChebyshev + m * _nChebyshev + n) =
                                            Transfer[bit[2]](i, l) * Transfer[bit[1]](j, m) * Transfer[bit[0]](k, n);
                                }
                            }
                        }
                    }
                }
            }
        }
     }

    void getScaledChebyNode(int _nChebyNode, VectorXd& _chebyNode, point& center, point& radius,
                            vector<point>& _scaledCnode) {
        for (int i = 0; i < _nChebyNode; ++i) {
            _scaledCnode.push_back(point(center.x + radius.x * _chebyNode(i),
                                         center.y + radius.y * _chebyNode(i),
                                         center.z + radius.z * _chebyNode(i)));
        }

    }

    void getCharge(int rootId) {
        node& n = t.dict[rootId];
        if(n.chargeComputed){
            return;
        }
        else{
            n.chargeComputed	=	true;
            n.charge		=	MatrixXd::Zero(n.nSource,1);
            for(unsigned long k=0;k<n.nSource; ++k){
                n.charge.row(k)	=	chargeTree.row(n.sourceIndex[k]);
            }
        }
    }

    void getTransferSourceParentToChildren(int _nChebyNode, vector<int>& _sourceIndex, point& _center, point& _radius,
                                     VectorXd& _chebyNode, MatrixXd& _tNode, MatrixXd& R) {

        int N = (int) _sourceIndex.size();
        VectorXd standlocation[3];
        standlocation[0].resize(N);
        standlocation[1].resize(N);
        standlocation[2].resize(N);
        for (int i = 0; i < N; ++i) {
            standlocation[0](i) = (t.sourceTree[_sourceIndex[i]].x - _center.x)/_radius.x;
            standlocation[1](i) = (t.sourceTree[_sourceIndex[i]].y - _center.y)/_radius.y;
            standlocation[2](i) = (t.sourceTree[_sourceIndex[i]].z - _center.z)/_radius.z;
        }

        MatrixXd Transfer[3];
        for (int k = 0; k < 3; ++k) {
            getStandardChebyPoly(_nChebyNode, N, standlocation[k], Transfer[k]);
            Transfer[k] = (2.0 * Transfer[k] * _tNode.transpose() - MatrixXd::Ones(N, _nChebyNode))/_nChebyNode;
        }
        int _rank = _nChebyNode * _nChebyNode * _nChebyNode;
        R = MatrixXd::Zero(N, _rank);
        for (int k = 0; k < N; ++k) {
            for (int i = 0; i < _nChebyNode; ++i) {
                for (int j = 0; j <_nChebyNode; ++j) {
                    for (int l = 0; l< _nChebyNode; ++l) {
                        R(k, l*_nChebyNode * _nChebyNode + j*_nChebyNode + i) =
                        Transfer[0](k, i) * Transfer[1](k, j) * Transfer[2](k, l);
                    }
                }
            }
        }
    }


    void getTransferTargetParentToChildren(int _nChebyNode, vector<int>& _targetIndex, point& _center, point& _radius,
                                           VectorXd& _chebyNode, MatrixXd& _tNode, MatrixXd& R) {

        int N = (int) _targetIndex.size();
        VectorXd standlocation[3];
        standlocation[0].resize(N);
        standlocation[1].resize(N);
        standlocation[2].resize(N);
        for (int i = 0; i < N; ++i) {
            standlocation[0](i) = (t.targetTree[_targetIndex[i]].x - _center.x)/_radius.x;
            standlocation[1](i) = (t.targetTree[_targetIndex[i]].y - _center.y)/_radius.y;
            standlocation[2](i) = (t.targetTree[_targetIndex[i]].z - _center.z)/_radius.z;
        }

        MatrixXd Transfer[3];
        for (int k = 0; k < 3; ++k) {
            getStandardChebyPoly(_nChebyNode, N, standlocation[k], Transfer[k]);
            Transfer[k] = (2.0 * Transfer[k] * _tNode.transpose() - MatrixXd::Ones(N, _nChebyNode))/_nChebyNode;
        }
        int _rank = _nChebyNode * _nChebyNode * _nChebyNode;
        R = MatrixXd::Zero(N, _rank);
        for (int k = 0; k < N; ++k) {
            for (int i = 0; i < _nChebyNode; ++i) {
                for (int j = 0; j <_nChebyNode; ++j) {
                    for (int l = 0; l< _nChebyNode; ++l) {
                        R(k, l*_nChebyNode * _nChebyNode + j*_nChebyNode + i) =
                                Transfer[0](k, i) * Transfer[1](k, j) * Transfer[2](k, l);
                    }
                }
            }
        }
    }


    void kernelEval(vector<point>& _source, vector<point>& _target, MatrixXd& K) {
        K = MatrixXd::Zero(_target.size(), _source.size());
        for (int _s = 0; _s < _source.size(); ++_s) {
            for (int _t = 0; _t < _target.size(); ++_t) {
                K(_t, _s) = this->eval(_source[_s], _target[_t]);
            }
        }
    }



    void kernelEvalIndex(vector<int>& _sourceIndex, vector<int>& _targetIndex, MatrixXd& K) {
        K = MatrixXd::Zero(_targetIndex.size(), _sourceIndex.size());
        for (int _s = 0; _s < _sourceIndex.size(); ++_s) {
            for (int _t = 0; _t < _targetIndex.size(); ++_t) {
                K(_t, _s) = this->eval(
                        this->t.sourceTree[_sourceIndex[_s]],
                        this->t.targetTree[_targetIndex[_t]]
                );
            }
        }
    }

    void kernelEvalChebyshev(int _M, vector<point>& _xv, int _N,  vector<point>& _yv, MatrixXd& K) {
        vector<point> sourceVec;
        vector<point> targetVec;
        K = MatrixXd::Zero(_M * _M * _M, _N* _N * _N);
        for (int k = 0; k < _M; k++) {
            for (int j = 0; j < _M; j++) {
                for (int i = 0; i < _M; i++) {
                    point np(_xv[i].x , _xv[j].y, _xv[k].z);
                    sourceVec.push_back(np);
                }
            }
        }

        for (int k = 0; k < _N; k++) {
            for (int j = 0; j < _N; j++) {
                for (int i = 0; i < _N; i++) {
                    point np(_yv[i].x , _yv[j].y, _yv[k].z);
                    targetVec.push_back(np);
                }
            }
        }

        kernelEval(sourceVec, targetVec, K);
    }


    void run(MatrixXd& potentialMatrix) {
        RUN("up-pass", upPass(0));
        potentialMatrix = MatrixXd::Zero(t.nTarget, 1);
        RUN("down-pass", downPass(0, potentialMatrix));
    }

    void upPass(int rootId) {
        node& n = t.dict[rootId];
        n.scaledCnode.clear();
        n.nodeCharge = MatrixXd::Zero(rank, 1);
        n.nodePotential = MatrixXd::Zero(rank, 1);
        getScaledChebyNode(nChebyshev, chebyNode, n.center, n.radius, n.scaledCnode);

        if (n.isLeaf) {
            // lazy
            getCharge(rootId);
            getTransferSourceParentToChildren(nChebyshev, n.sourceIndex, n.center, n.radius, chebyNode, tNode, n.R);
            getTransferTargetParentToChildren(nChebyshev, n.targetIndex, n.center, n.radius, chebyNode, tNode, n.L);
            n.nodeCharge += n.R.transpose() * n.charge;
        }
        else {
            for (int i = 0; i < 8; ++i) {
                upPass(n.child[i]);
                if (!t.dict[n.child[i]].isEmpty) {
                    n.nodeCharge += R[i].transpose() * t.dict[n.child[i]].nodeCharge;
                }
            }
        }
    }

    void downPass(int rootId, MatrixXd& potential) {
        node& n = t.dict[rootId];
        MatrixXd K;

        VectorXd temp;
        if (n.parent != -1) {
            /*
             * V list
             */
            for (int i : n.vList) {
                if (!t.dict[i].isEmpty) {
                    kernelEvalChebyshev(nChebyshev, t.dict[i].scaledCnode, nChebyshev, n.scaledCnode, K);
                    n.nodePotential += K * t.dict[i].nodeCharge;
                }
            }
            /*
             * X List
             */
            for (int i : n.xList) {
                if (!t.dict[i].isEmpty) {
                    kernelEvalChebyshev(nChebyshev, t.dict[i].scaledCnode, nChebyshev, n.scaledCnode, K);
                    n.nodePotential += K * t.dict[i].nodeCharge;
                }
            }

            /*
             * L2L
             */
            node& p = t.dict[n.parent];
            n.nodePotential += this->R[n.nodeIndex] * p.nodePotential;

        }

        if (n.isLeaf && n.nTarget != 0) {
            n.potential = MatrixXd::Zero(n.nTarget, 1);

            /*
             * U List
             */
            for (int i : n.uList) {
                if (!t.dict[i].isEmpty) {
                    getCharge(i);
                    kernelEvalIndex(t.dict[i].sourceIndex, n.targetIndex, K);
                    n.potential += K * t.dict[i].charge;
                }
            }

            /*
             * W List
             */

            for (int i : n.wList) {
                if (!t.dict[i].isEmpty) {
                    getCharge(i);
                    kernelEvalIndex(t.dict[i].sourceIndex, n.targetIndex, K);
                    n.potential += K * t.dict[i].charge;
                }
            }

            /*
             * L2T
             */
            n.potential += n.L * n.nodePotential;

            /*
             * Finalize
             */
            for (int i = 0; i < n.nTarget; i++) {
                potential.row(n.targetIndex[i]) += n.potential.row(i);
            }
        }

        if (!n.isLeaf) {
            for (int i = 0; i < 8; ++i) {
                downPass(n.child[i], potential);
            }
        }
    }
};


#endif //FMM_KERNEL_H
