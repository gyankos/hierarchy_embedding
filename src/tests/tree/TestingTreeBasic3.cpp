//
// Created by giacomo on 29/12/19.
//

#include "tests/tree/TestingTreeBasic3.h"
#include "tests/tree/TestingTreeBasic1.h"
#include "concept_vector/ConceptVector.h"
#include <math_utils.h>

TestingTreeBasic3::TestingTreeBasic3(size_t maximumBranchingFactor, size_t maximumHeight, double beta, double alpha)
        : TestingTree{maximumBranchingFactor, maximumHeight},
          beta(beta), alpha(alpha) {}

void TestingTreeBasic3::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) {
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
        //std::cout << x <<  "    " << tmp << std::endl ;
        nodes += tmp;
    }
    keyValueMap.insert(std::make_pair("", ConceptVector::multiple_descent_concept_problem_vector(tree, std::vector<size_t>{}, beta, alpha)));
    for (const std::vector<size_t> &x: ls) {
        keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::multiple_descent_concept_problem_vector(tree, x, beta, alpha)));
    }
}

std::vector<double> TestingTreeBasic3::getVectorRepresentation(const std::vector<size_t> &current) {
    return keyValueMap[size_vector_to_string(current)];
}

double TestingTreeBasic3::similarity(const std::vector<double> &lhs, const std::vector<double> &rhs) { return cosine_similarity(lhs, rhs); }

void TestingTreeBasic3::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    std::string currentString{size_vector_to_string(current)};
    const std::vector<double>& currentVector = getVectorRepresentation(current);
    for (auto & it : keyValueMap) {
        //if (it.first != currentString) {
        //double  score = similarity(currentVector, it.second);
        //std::cout << "sim(" << currentString << "," << it.first << ")=" << score << std::endl ;
        map.add(similarity(currentVector, it.second), it.first);
        //}
    }
}
