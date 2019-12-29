//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_COUT_UTILS_H
#define HIERARCHY_TESTS_COUT_UTILS_H

#include <iostream>
#include <vector>
#include <set>

std::ostream& operator<<(std::ostream& os, const std::vector<size_t>& x);

std::ostream& operator<<(std::ostream& os, const std::vector<double>& x);

std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<size_t>>& x);

std::ostream &operator<<(std::ostream &os, const std::set<std::string> &x);

#endif //HIERARCHY_TESTS_COUT_UTILS_H
