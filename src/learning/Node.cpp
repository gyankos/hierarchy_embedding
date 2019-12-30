//
// Created by giacomo on 29/12/19.
//

#include "learning/Node.h"

Node::Node(const int id, const int level) : id_(id), level_(level) {}

void Node::set_level(const int level) { level_ = level; }

void Node::AddParent(int p_idx) {
    for (int i = 0; i < parent_idx_.size(); ++i) {
        if (parent_idx_[i] ==  p_idx) { return; }
    }
    parent_idx_.push_back(p_idx);
}

void Node::AddChild(int c_idx) {
    for (int i = 0; i < child_idx_.size(); ++i) {
        if (child_idx_[i] ==  c_idx) { return; }
    }
    child_idx_.push_back(c_idx);
}

bool Node::operator<(const Node &node) {
    return level_ < node.level();
}
