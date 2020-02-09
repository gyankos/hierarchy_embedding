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

//std::ostream &operator<<(std::ostream &os, const std::set<std::string> &x);

template <typename t> std::ostream &operator<<(std::ostream &os, const std::set<t> &x) {
    os << "{ ";
    size_t i = 0, n = x.size();
    for (auto it = x.begin(); it != x.end(); ++it) {
        os << *it << (i == n-1 ? " }" : ", "); i++;
    }
    if (x.empty()) os << "}";
    return os;
}

#endif //HIERARCHY_TESTS_COUT_UTILS_H
