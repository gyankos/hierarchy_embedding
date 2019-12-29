//
// Created by giacomo on 29/12/19.
//

#include "concept_vector/ConceptVector.h"

std::vector<double> ConceptVector::relevancy_vector(naryTree &tree, const std::vector<size_t> &path) {
    std::vector<double> to_return{};
    to_return.resize((size_t)tree.getDescendantsNo()+1, 0.0);
    naryTree* currentNode = &tree;
    to_return[0] = 1;
    if (path.empty()) {
        for (size_t i = 0; i<tree.getChildNo(); i++) {
            tree.getIthChild(i).setBasicTreeVector1(to_return, tree.getChildNo());
        }
    } else {
        size_t parentChildNo = currentNode->getChildNo();
        for (const size_t& currPos : path) {
            currentNode = &currentNode->getIthChild(currPos-1);
            to_return[currentNode->getId()] = parentChildNo;
            parentChildNo = currentNode->getChildNo();
        }
    }
    return to_return;
}

std::vector<double>
ConceptVector::local_density_problem_vector(naryTree &tree, const std::vector<size_t> &path, double beta) {
    std::vector<double> to_return{};
    to_return.resize((size_t)tree.getDescendantsNo()+1, 0.0);
    std::function<double(double)> fun = [beta](double x) { return (1-std::pow(beta, x))/(1-beta); };
    naryTree* currentNode = &tree;
    to_return[0] = 1;
    if (path.empty()) {
        for (size_t i = 0; i<tree.getChildNo(); i++) {
            tree.getIthChild(i).setBasicTreeVector2(to_return, tree.getChildNo(), fun);
        }
    } else {
        size_t parentChildNo = currentNode->getChildNo();
        for (const size_t& currPos : path) {
            currentNode = &currentNode->getIthChild(currPos-1);
            to_return[currentNode->getId()] = fun(parentChildNo);
            parentChildNo = currentNode->getChildNo();
        }
    }
    return to_return;
}

std::vector<double>
ConceptVector::multiple_descent_concept_problem_vector(naryTree &tree, const std::vector<size_t> &path, double beta,
                                                       double alpha) {
    std::vector<double> to_return{};
    to_return.resize((size_t)tree.getDescendantsNo()+1, 0.0);
    std::function<double(double)> fun = [beta](double x) { return (beta==1) ? beta*x : (1-std::pow(beta, x))/(1-beta); };
    naryTree* currentNode = &tree;
    to_return[0] = 1;
    if (path.empty()) {
        for (size_t i = 0; i<tree.getChildNo(); i++) {
            tree.getIthChild(i).setBasicTreeVector3(to_return, tree.getChildNo(), fun, alpha, 1);
        }
    } else {
        size_t parentChildNo = currentNode->getChildNo();
        double current = alpha;
        for (const size_t& currPos : path) {
            currentNode = &currentNode->getIthChild(currPos-1);
            to_return[currentNode->getId()] = current * fun(parentChildNo);
            parentChildNo = currentNode->getChildNo();
            current *= alpha;
        }
    }
    return to_return;
}
