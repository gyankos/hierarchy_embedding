//
// Created by giacomo on 29/12/19.
//

#include "tests/tree/TestingTreeBasic1.h"
#include "concept_vector/ConceptVector.h"
#include <math_utils.h>

void TestingTreeBasic1::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &ls) {
    for (const std::vector<size_t> &x: ls) {
        size_t tmp = tree.addChild(x, id);
        nodes += tmp;
    }
    keyValueMap.insert(std::make_pair("", ConceptVector::relevancy_vector(tree, std::vector<size_t>{})));
    for (const std::vector<size_t> &x: ls) {
        keyValueMap.insert(std::make_pair(size_vector_to_string(x), ConceptVector::relevancy_vector(tree, x)));
    }
}

std::vector<double> TestingTreeBasic1::getVectorRepresentation(const std::vector<size_t> &current) {
    return keyValueMap[size_vector_to_string(current)];
}

void TestingTreeBasic1::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    std::string currentString{size_vector_to_string(current)};
    const std::vector<double>& currentVector = getVectorRepresentation(current);
    for (auto & it : keyValueMap) {
        //if (it.first != currentString) {
        //double  score = similarity(currentVector, it.second);
        // std::cout << "sim(" << currentString << "," << it.first << ")=" << score << std::endl ;
        map.add(similarity(currentVector, it.second), it.first);
        //}
    }
}

double TestingTreeBasic1::similarity(const std::vector<double> &lhs, const std::vector<double> &rhs) { return cosine_similarity(lhs, rhs); }

TestingTreeBasic1::TestingTreeBasic1(size_t maximumBranchingFactor, size_t maximumHeight) : TestingTree{maximumBranchingFactor, maximumHeight} {}
