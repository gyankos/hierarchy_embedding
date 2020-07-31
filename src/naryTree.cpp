//
// Created by giacomo on 29/12/19.
//

#include <iostream>
#include "naryTree.h"

naryTree::naryTree(const size_t &root) : root(root), id{0}, noAllChildren{0}, parent{nullptr} {}

naryTree::naryTree(const naryTree &x) : root{x.root}, children{x.children}, id{x.id}, noAllChildren{x.noAllChildren}, parent{x.parent} {}

naryTree &naryTree::operator=(const naryTree &x) {
    this->root = x.root;
    this->children = x.children;
    this->id = x.id;
    this->noAllChildren = x.noAllChildren;
    this->parent = x.parent;
    return *this;
}

size_t naryTree::addChild(const std::vector<size_t> &path, size_t pathPos, size_t &global_counter) {
    if (path.size() <= pathPos)
        return 0;
    else {
        size_t x = path[pathPos];
        while (children.size()<x) {
            children.emplace_back(0);
        }
        bool  test = (children[x-1].root != x);
        if (test) {
            children[x-1].root = x;
            children[x-1].id = ++global_counter;
        }
        children[x-1].parent = this;
        size_t toRet =  (test ? 1 : 0) + children[x - 1].addChild(path, pathPos+1, global_counter);
        noAllChildren += toRet;
        return toRet;
    }
}

bool naryTree::addChild(const std::vector<size_t> &path, size_t &global_counter) {
    return addChild(path, 0, global_counter);
}

void naryTree::setBasicTreeVector2(std::vector<double> &vector, size_t parentSize, std::function<double(double)> fun) {
    vector[id] = fun(parentSize);
    for (auto& child: children)
        child.setBasicTreeVector2(vector, children.size(), fun);
}

void naryTree::setBasicTreeVector3(std::vector<double> &vector, size_t parentSize, std::function<double(double)> fun,
                                   double alpha, size_t level) {
    //std::cout << levelFactor << "^" << level << "*" << parentSize << " = " << std::pow(levelFactor, level) * parentSize<< std::endl;
    vector[id] = std::pow(alpha, level) * fun(parentSize);
    for (auto& child: children)
        child.setBasicTreeVector3(vector, children.size(), fun, alpha, level + 1);
}
