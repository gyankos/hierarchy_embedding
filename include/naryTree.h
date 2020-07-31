//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_NARYTREE_H
#define HIERARCHY_TESTS_NARYTREE_H


#include <string>
#include <vector>
#include <functional>
#include <tgmath.h>

class naryTree {
    size_t root, id;
    std::vector<naryTree> children;
    size_t noAllChildren;
    naryTree* parent;

public:
    naryTree(const size_t &root);
    naryTree(const naryTree& x);
    naryTree& operator=(const naryTree& x);
    size_t getParentChildNo() const {
        return parent ? parent->getChildNo() : 1;
    }
    size_t getChildNo() const  { return children.size(); }
    naryTree& getIthChild(size_t u) { return children[u]; }
    bool addChild(const std::vector<size_t> &path, size_t &global_counter);
    size_t getId() const { return id; }
    size_t getDescendantsNo() const { return  noAllChildren; }

    void setBasicTreeVector1(std::vector<double> &vector, size_t parentSize) {
        vector[id] = parentSize;
        for (auto& child: children)
            child.setBasicTreeVector1(vector, children.size());
    }


    void setBasicTreeVector2(std::vector<double> &vector, size_t parentSize, std::function<double(double)> fun);

    void setBasicTreeVector3(std::vector<double> &vector, size_t parentSize, std::function<double(double)> fun, double alpha, size_t level);

private:
    size_t addChild(const std::vector<size_t> &path, size_t pathPos, size_t &global_counter);
};


#endif //HIERARCHY_TESTS_NARYTREE_H
