//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGTREEBASIC3_H
#define HIERARCHY_TESTS_TESTINGTREEBASIC3_H

#include <vector>
#include "tests/TestingTree.h"
#include <naryTree.h>



class TestingTreeBasic3 : public TestingTree<std::vector<double>>  {

    std::map<std::string, std::vector<double>> keyValueMap;
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;
    double beta, alpha;

public:
    TestingTreeBasic3(size_t maximumBranchingFactor, size_t maximumHeight, double beta, double alpha);

protected:

    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) override;;
    std::vector<double> getVectorRepresentation(const std::vector<size_t>& current);
    double similarity(const std::vector<double>& lhs, const std::vector<double>& rhs) override;
    void generateTopKCandidates(PollMap<double,std::string>& map, const std::vector<size_t>& current);
    void finalizeDataIngestion() override {}
};

/*
 *

double test_basic_3(const std::vector<std::vector<size_t>> &ls) {
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
        //std::cout << x <<  "    " << tmp << std::endl ;
        nodes += tmp;
    }
    std::cout << ConceptVector::multiple_descent_concept_problem_vector(tree, {}, 0.75, 0.5) << std::endl;
    for (const std::vector<size_t> &x: ls) {
        std::cout << x << "    " <<  ConceptVector::multiple_descent_concept_problem_vector(tree, x, 0.75, 0.5) << std::endl;
    }
    return 1.0;
}



double test_basic_4(const std::vector<std::vector<size_t>> &ls) {
    naryTree tree{0};
    size_t nodes = 0;
    size_t id = 0;
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
        //std::cout << x <<  "    " << tmp << std::endl ;
        nodes += tmp;
    }
    std::cout << ConceptVector::multiple_descent_concept_problem_vector(tree, {}, 1, 0.5) << std::endl;
    for (const std::vector<size_t> &x: ls) {
        std::cout << x << "    " <<  ConceptVector::multiple_descent_concept_problem_vector(tree, x, 1, 0.5) << std::endl;
    }
    return 1.0;
}

 */

#endif //HIERARCHY_TESTS_TESTINGTREEBASIC3_H
