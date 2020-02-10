//
// Created by giacomo on 04/01/20.
//

#ifndef HIERARCHY_TESTS_TESTINGTREEPAR_H
#define HIERARCHY_TESTS_TESTINGTREEPAR_H


#include "tests/TestingTree.h"
#include "par2hier/Par2Hier.h"

class TestingTreePar : public TestingTree<std::vector<size_t>> {
    size_t id = 0; size_t nodes = 0;
    std::map<std::string, std::vector<double>> originalVectorMap;
    naryTree tree{0};
    Par2Hier p2h;
    std::vector<std::vector<size_t>> allPossiblePaths;
    Proposal safekeep;

    std::mutex g_pages_mutex;

protected:
    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) override;

    std::vector<size_t> getVectorRepresentation(const std::vector<size_t> &current) override;

    double similarity(const std::vector<size_t> &lhs, const std::vector<size_t> &rhs) override;

    void generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) override;

    void finalizeDataIngestion() override;

public:
    TestingTreePar(size_t maximumBranchingFactor, size_t maximumHeight, size_t k, method m, double distanceFactor, double decayFactor) : TestingTree(maximumBranchingFactor,
                                                                                                                                                     maximumHeight), p2h{originalVectorMap, k, m, tree}, safekeep{(double )maximumBranchingFactor, distanceFactor, decayFactor} {}

};


#endif //HIERARCHY_TESTS_TESTINGTREEPAR_H
