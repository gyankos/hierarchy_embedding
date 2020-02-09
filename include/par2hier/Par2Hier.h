//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_PAR2HIER_H
#define HIERARCHY_TESTS_PAR2HIER_H

#include <map>
#include <vector>
#include <naryTree.h>
#include <string_utils.h>
#include <Eigen/Dense>
#include <math_utils.h>
#include <proposal/Proposal.h>

/**
 * C++ Porting of the Java code provided in https://github.com/tteofili/par2hier/ for the paper:
 *
 * Tommaso Teofili: "par2hier: towards vector representations for hierarchical content" ICCS 2017
 *
 *
 */
enum method {
    CLUSTER,SUM
};

class Par2Hier {
    std::map<std::string, std::vector<double>>& originalVectorMap;
    std::map<std::string, std::vector<double>> memoizationMap;
    size_t k;
    method m;
    naryTree& tree;

    static std::vector<double> orginalVector(const Eigen::MatrixXd::RowXpr &l);

    static std::vector<double> originalVector(const Eigen::MatrixXd::ColXpr &l);
public:

    Par2Hier(std::map<std::string, std::vector<double>> &originalVectorMap,  size_t k,
             method m,
             naryTree& tree);

    const std::vector<double>& getPar2HierVector(const std::vector<size_t>& path);

    /*void getTruncatedSVD(const Eigen::MatrixXd& matrix, Eigen::MatrixXd& outMatrix) {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);


        Eigen::MatrixXd matrixU = svd.matrixU();
        Eigen::MatrixXd matrixV = svd.matrixV();
        Eigen::MatrixXd matrixS = svd.singularValues();

        Eigen::MatrixXd truncatedU = matrixU.block(0,0,matrixU.rows()-1,k);
        Eigen::MatrixXd truncatedV = matrixV.block(0,0,k,matrixV.cols()-1);
        Eigen::MatrixXd truncatedS = matrixV.block(0,0,k,k);

        outMatrix = truncatedU * (matrixS * matrixV);
    }*/

    void getTruncatedVT(const Eigen::MatrixXd& matrix, Eigen::MatrixXd& outMatrix, const size_t kDim);

};

#endif //HIERARCHY_TESTS_PAR2HIER_H
