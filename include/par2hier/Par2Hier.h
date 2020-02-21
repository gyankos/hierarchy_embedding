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
    std::mutex g_pages_mutex;
    size_t k;
    method m;
    naryTree& tree;

    static std::vector<double> orginalVector(const Eigen::MatrixXd::RowXpr &l);
    static std::vector<double> originalVector(const Eigen::MatrixXd::ColXpr &l);
public:

    Par2Hier(std::map<std::string, std::vector<double>> &originalVectorMap,  size_t k,
             method m,
             naryTree& tree);

    const std::vector<double>& getPar2HierVector(const std::vector<size_t> &path);

    void getTruncatedVT(const Eigen::MatrixXd& matrix, Eigen::MatrixXd& outMatrix, const size_t kDim);

};

#endif //HIERARCHY_TESTS_PAR2HIER_H
