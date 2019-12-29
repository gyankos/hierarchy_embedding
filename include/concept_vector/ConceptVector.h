//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_CONCEPTVECTOR_H
#define HIERARCHY_TESTS_CONCEPTVECTOR_H


#include <vector>
#include <naryTree.h>
#include <iostream>
#include <cmath>

class ConceptVector {
public:
    static std::vector<double> relevancy_vector(naryTree& tree, const std::vector<size_t>& path);
    static std::vector<double> local_density_problem_vector(naryTree &tree, const std::vector<size_t> &path, double beta);
    static std::vector<double>  multiple_descent_concept_problem_vector(naryTree &tree, const std::vector<size_t> &path, double beta, double alpha);
};


#endif //HIERARCHY_TESTS_CONCEPTVECTOR_H
