//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_PROPOSAL_UTILS_H
#define HIERARCHY_TESTS_PROPOSAL_UTILS_H

#include <vector>
#include <cmath>
#include <future>

int indexOfSubList(std::vector<size_t>& array, std::vector<size_t>& subarray);
std::vector<std::vector<std::vector<size_t>>> generateCompleteSubgraph(size_t maximumBranchingFactor, size_t maximumHeight);
bool subArrayOf(std::vector<size_t>& array, std::vector<size_t>& subarray);

/**
     * Implements the euclidean distance between two coordinate vectors
     * @param left      Left coordinate vector
     * @param right     Right coordinate vector
     * @return          Distance
     */
double euclideanDistance(const std::vector<double> &left, const std::vector<double> &right);

#endif //HIERARCHY_TESTS_PROPOSAL_UTILS_H
