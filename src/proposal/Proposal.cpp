//
// Created by giacomo on 29/12/19.
//

#include "proposal/Proposal.h"

Proposal::Proposal(double maximumBranchingFactor, double distanceFactor, double decayFactor) : maximumBranchingFactor(
        maximumBranchingFactor), distanceFactor(distanceFactor), decayFactor(decayFactor) {}

std::vector<double> Proposal::fromHierarchyVectorToEuclideanVector(const std::vector<size_t> &hierarchyVector) {
    std::vector<double> d;
    d.resize(maximumBranchingFactor, 0.0);
    double currentDecay = 1.0;
    for (size_t i = 0, hierarchyVectorLength = hierarchyVector.size(); i < hierarchyVectorLength; i++) {
        size_t vi = hierarchyVector[i] - 1;
        d[vi] = d[vi] + (distanceFactor / currentDecay);
        currentDecay *= decayFactor;
    }
    return d;
}

double Proposal::thresholdValue(const std::vector<size_t> &left, const std::vector<size_t> &right) {
    int h2 = left.size(), h1 = right.size();
    if (h2 < h1) {
        int tmp = h2;
        h2 = h1;
        h1 = tmp;
    }
    return distanceFactor * (std::pow(decayFactor, h2-h1+1) - 1.0) / (std::pow(decayFactor, h2) * (decayFactor - 1.0));
}

double Proposal::distance(const std::vector<size_t> &left, const std::vector<size_t> &right) {
    auto l = fromHierarchyVectorToEuclideanVector(left);
    auto r = fromHierarchyVectorToEuclideanVector(right);
    return euclideanDistance(l, r);
}

bool Proposal::isConsistent(const std::vector<size_t> &left, const std::vector<size_t> &right) {
    return distance(left, right) <= thresholdValue(left, right);
}
