//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_PROPOSAL_H
#define HIERARCHY_TESTS_PROPOSAL_H

#include <vector>
#include <cmath>
#include <proposal/proposal_utils.h>

class Proposal {

    /**
     * Given a hierarchical tree extracted from a DAG, its maximum assocated branching factor
     */
double maximumBranchingFactor;

    /**
     * Maximum distance represented by the distance between the root (empty vector) and one of the nodes at the fist level
     */
double distanceFactor;

    /**
     * Factor reducing the distance between the one node from the i-th level to the children at the (i+1)-level greater
     * at each level step, such that the nodes from different branches will not be overlapped.
     */
double decayFactor;

public:
    Proposal(double maximumBranchingFactor, double distanceFactor, double decayFactor);


    /**
     * Returns a vector that can be used to compute the euclidean distance from the hierarchical vector
     * @param hierarchyVector   Vector expressed as in the Petermann paper
     * @return
     */
    std::vector<double> fromHierarchyVectorToEuclideanVector(const std::vector<size_t> &hierarchyVector);

    /**
     * Distance threshold values for two given elements
     * @param left
     * @param right
     * @return
     */
 double thresholdValue(const std::vector<size_t> &left, const std::vector<size_t> &right);

    /**
     * @param left      Hierarchy vector
     * @param right     Hierarchy vector
     * @return          Distance
     */
 double distance(const std::vector<size_t> &left, const std::vector<size_t> &right);

    /**
     * This distance is defined as follows: isConsistent(left, right) <=> subArrayOf(left,right) || subArrayOf(right,left)
     *
     * @param left    path index vector
     * @param right   path index vector
     * @return   Whether the two elements are consistent or not.
     */
 bool isConsistent(const std::vector<size_t> &left, const std::vector<size_t> &right);

 double rankingMetric(const std::vector<size_t>& left, const std::vector<size_t >& right) {
     return isConsistent(left, right) ? distance(left, right) : std::numeric_limits<double>::max();
 }

};


#endif //HIERARCHY_TESTS_PROPOSAL_H
