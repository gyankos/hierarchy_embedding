//
// Created by giacomo on 04/01/20.
//

#include "tests/tree/TestingTreePar.h"

void TestingTreePar::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) {
    for (const std::vector<size_t> &x: subgraph_as_paths) {
        size_t tmp = tree.addChild(x, id);
        originalVectorMap.insert(std::make_pair(size_vector_to_string(x), safekeep.fromHierarchyVectorToEuclideanVector(x)));
        nodes += tmp;
    }
    allPossiblePaths.insert(allPossiblePaths.end(), subgraph_as_paths.begin(), subgraph_as_paths.end());
}

std::vector<size_t> TestingTreePar::getVectorRepresentation(const std::vector<size_t> &current) {
    return current;
}

double TestingTreePar::similarity(const std::vector<size_t> &lhs, const std::vector<size_t> &rhs) {
    std::vector<size_t> L{lhs}, R{rhs};
    auto x = p2h.getPar2HierVector(L);
    auto y = p2h.getPar2HierVector(R);
    double distance =  euclideanDistance(x, y);//safekeep.isConsistent(lhs, rhs) ? euclideanDistance(x, y) : std::numeric_limits<double>::max();
    // From now on, same as mine
    double normalizedDistance = distance / (1+distance);
    double similairity = 1 - normalizedDistance;
    //std::cerr << "sim(" << size_vector_to_string(lhs) << ", " << size_vector_to_string(rhs)<< ") = " << similairity << " via distance " << distanceMetric << " normalized as " << normalizedDistance << std::endl;
    return similairity;
}

void TestingTreePar::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    std::string currentString{size_vector_to_string(current)};
    for (auto & it : allPossiblePaths) {
        //if (it.first != currentString) {
        double  score = similarity(current, it);
        // std::cout << "sim(" << currentString << "," << it.first << ")=" << score << std::endl ;
        map.add(score, size_vector_to_string(it));
        //}
    }
}

void TestingTreePar::finalizeDataIngestion() {
    //p2h.getPar2HierVector({}); // This should initialize the whole metrics, as the root will recursively call all of its childs.
    // As a consequence, the memoization map within it will be stored.
}
