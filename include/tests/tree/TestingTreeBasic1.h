//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGTREEBASIC1_H
#define HIERARCHY_TESTS_TESTINGTREEBASIC1_H

#include <vector>
#include "tests/TestingTree.h"
#include <naryTree.h>

class TestingTreeBasic1 : public TestingTree<std::vector<double>>  {
    std::map<std::string, std::vector<double>> keyValueMap;
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;

public:
    TestingTreeBasic1(size_t maximumBranchingFactor, size_t maximumHeight);

protected:
    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) override;;
    std::vector<double> getVectorRepresentation(const std::vector<size_t>& current);
    double similarity(const std::vector<double>& lhs, const std::vector<double>& rhs) override;
    void generateTopKCandidates(PollMap<double,std::string>& map, const std::vector<size_t>& current);
    void finalizeDataIngestion() override {}
};

#endif //HIERARCHY_TESTS_TESTINGTREEBASIC1_H
