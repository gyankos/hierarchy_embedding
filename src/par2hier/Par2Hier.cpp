//
// Created by giacomo on 30/12/19.
//

#include "par2hier/Par2Hier.h"


std::vector<double> Par2Hier::orginalVector(const Eigen::MatrixXd::RowXpr &l) {
    std::vector<double> toret;
    for (size_t i = 0, n = l.size(); i<n; i++) {
        toret.emplace_back(l[i]);
    }
    return toret;
}

std::vector<double> Par2Hier::originalVector(const Eigen::MatrixXd::ColXpr &l) {
    std::vector<double> toret;
    for (size_t i = 0, n = l.size(); i<n; i++) {
        toret.emplace_back(l[i]);
    }
    return toret;
}

Par2Hier::Par2Hier(std::map<std::string, std::vector<double>> &originalVectorMap, size_t k, method m, naryTree &tree) : originalVectorMap(originalVectorMap), k(k), m(m), tree{tree} {}

const std::vector<double> &Par2Hier::getPar2HierVector(const std::vector<size_t> &path) {
    std::string currentPath = size_vector_to_string(path);


    g_pages_mutex.lock();
    auto it = memoizationMap.find(currentPath);
    if (it != memoizationMap.end()) {
        std::vector<double>& returno = it->second;
        g_pages_mutex.unlock();
        return returno;
    }
    g_pages_mutex.unlock();

    naryTree* currentNode = &tree;
    size_t parentChildNo = currentNode->getChildNo();

    // Traversing the tree, to get the node of reference
    if (!path.empty()) {
        for (const size_t& currPos : path) {
            currentNode = &currentNode->getIthChild(currPos-1);
            parentChildNo = currentNode->getChildNo();
        }
    }

    // If the selected node is a leaf, then return the same vector from the embedding of choice
    if (parentChildNo == 0) {
        std::lock_guard<std::mutex> lg(g_pages_mutex);
        auto it2 = memoizationMap.find(currentPath);
        if (it2 == memoizationMap.end()) {
            assert(originalVectorMap.find(size_vector_to_string(path)) != originalVectorMap.end());
            return memoizationMap.insert(std::make_pair(currentPath, originalVectorMap[size_vector_to_string(path)])).first->second;
        } else {
            return it2->second;
        }
    }
    else {
        std::vector<double> current = originalVectorMap[currentPath];
        Eigen::MatrixXd matrix(parentChildNo, k);
        for (size_t i = 1, n = parentChildNo; i<=n; i++) {
            std::vector<size_t> childPath{path};
            childPath.emplace_back(i);
            std::vector<double> childVector = getPar2HierVector(childPath);
            for (size_t j = 0; j< std::min(childVector.size(), k); j++) {
                matrix.row(i-1)[j] = childVector[j];
            }
        }

        Eigen::MatrixXd centroids;
        if (parentChildNo > k) {
            getTruncatedVT(matrix, centroids, k);
        } else if (parentChildNo == 1) {
            centroids = matrix.row(0);
        } else {
            getTruncatedVT(matrix, centroids, 1);
        }

        switch (m) {
            case CLUSTER: {
                Eigen::MatrixXd matrix2{centroids.rows()+1, current.size()};
                for (size_t j = 0; j<current.size(); j++) {
                    matrix2.row(0)[j] = current[j];
                }
                for (size_t i = 1, n = centroids.rows(); i<n; i++) {
                    size_t currRowSize = std::min(matrix2.row(0).size(), centroids.row(i-1).size());
                    for (size_t j = 0; j<currRowSize; j++) {
                        matrix2.row(i)[j] = centroids.row(i-1)[j];
                    }
                }
                Eigen::MatrixXd hvMatrix;
                getTruncatedVT(matrix2, hvMatrix, 1);
                current = Par2Hier::originalVector(hvMatrix.col(0));
            }
                break;
            case SUM: {
                for (size_t i = 0; i<centroids.rows(); i++) {
                    current += Par2Hier::orginalVector(centroids.row(i));
                }
            }
                break;
        }

        {
            std::lock_guard<std::mutex> lg(g_pages_mutex);
            auto it23 = memoizationMap.find(currentPath);
            if ((it23) == memoizationMap.end()) {
                return memoizationMap.insert(std::make_pair(currentPath, current)).first->second;
            } else {
                return it23->second;
            }
        }
    }
}

void Par2Hier::getTruncatedVT(const Eigen::MatrixXd &matrix, Eigen::MatrixXd &outMatrix, const size_t kDim) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd matrixV = svd.matrixV();
    Eigen::MatrixXd truncatedV = matrixV.block(0,0,kDim,matrixV.cols()-1);
    outMatrix = truncatedV;
}
