//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGBASIC1_H
#define HIERARCHY_TESTS_TESTINGBASIC1_H

#include <vector>
#include "Testing.h"
#include <naryTree.h>

class TestingBasic1 : public Testing<std::vector<double>>  {
    std::map<std::string, std::vector<double>> keyValueMap;
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;

public:
    TestingBasic1(size_t maximumBranchingFactor, size_t maximumHeight);

protected:
    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) override;;
    std::vector<double> getVectorRepresentation(const std::vector<size_t>& current);
    double similarity(const std::vector<double>& lhs, const std::vector<double>& rhs) override;
    void generateTopKCandidates(PollMap<double,std::string>& map, const std::vector<size_t>& current);
    void finalizeDataIngestion() override {}
};

/*
 * double test_basic_1(const std::vector<std::vector<size_t>> &ls) {
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
        //std::cout << x <<  "    " << tmp << std::endl ;
        nodes += tmp;
    }
    std::cout << ConceptVector::relevancy_vector(tree, {}) << std::endl;
    for (const std::vector<size_t> &x: ls) {
       std::cout << x << "    " <<  ConceptVector::relevancy_vector(tree, x) << std::endl;
    }
    return 1.0;
}

 */

#endif //HIERARCHY_TESTS_TESTINGBASIC1_H
