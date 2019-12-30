//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_NODE_H
#define HIERARCHY_TESTS_NODE_H

#include <vector>

class Node {
public:
    Node(const int id, const int level);;
    ~Node() {};

    void set_level(const int level);;
    void AddParent(int p_idx);

    void AddChild(int c_idx);

    inline const std::vector<int>& parent_idx() const { return parent_idx_; }
    inline const std::vector<int>& child_idx() const { return child_idx_; }
    inline const int id() const { return id_; }
    inline const int level() const { return level_; }

    bool operator<(const Node& node);

private:
    // entity or category id
    int id_;
    int level_;

    // parents' index in hierarchy
    std::vector<int> parent_idx_;
    // children's index in hierarchy
    std::vector<int> child_idx_;
};


#endif //HIERARCHY_TESTS_NODE_H
