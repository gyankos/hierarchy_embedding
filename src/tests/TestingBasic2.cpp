//
// Created by giacomo on 29/12/19.
//

#include "tests/TestingBasic2.h"
#include "concept_vector/ConceptVector.h"
#include <math_utils.h>

TestingBasic2::TestingBasic2(size_t maximumBranchingFactor, size_t maximumHeight, double beta)
        : Testing{maximumBranchingFactor, maximumHeight},
          beta(beta) {}

void TestingBasic2::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) {
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
//std::cout << x <<  "    " << tmp << std::endl ;
        nodes += tmp;
    }
    keyValueMap.insert(std::make_pair("", ConceptVector::local_density_problem_vector(tree, std::vector<size_t>{}, beta)));
    for (const std::vector<size_t> &x: ls) {
        keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::local_density_problem_vector(tree, x, beta)));
    }
}

std::vector<double> TestingBasic2::getVectorRepresentation(const std::vector<size_t> &current) {
    return keyValueMap[size_vector_to_string(current)];
}

double TestingBasic2::similarity(const std::vector<double> &lhs, const std::vector<double> &rhs) { return cosine_similarity(lhs, rhs); }

void TestingBasic2::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
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
