//
// Created by giacomo on 29/12/19.
//

#include "cout_utils.h"

std::ostream &operator<<(std::ostream &os, const std::vector<size_t> &x) {
    os << "{ ";
    for (size_t i = 0, n = x.size(); i<n; i++) {
        os << x[i] << (i == n-1 ? " }" : ", ");
    }
    if (x.empty()) os << "}";
    return os;
}

std::ostream &operator<<(std::ostream &os, const std::vector<double> &x) {
    os << "{ ";
    for (size_t i = 0, n = x.size(); i<n; i++) {
        os << x[i] << (i == n-1 ? " }" : ", ");
    }
    if (x.empty()) os << "}";
    return os;
}

std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<size_t>> &x) {
    os << "{ ";
    for (size_t i = 0, n = x.size(); i<n; i++) {
        os << x[i] << (i == n-1 ? " }" : ", ");
    }
    if (x.empty()) os << "}";
    return os;
}

std::ostream &operator<<(std::ostream &os, const std::set<std::string> &x) {
    os << "{ ";
    size_t i = 0, n = x.size();
    for (auto it = x.begin(); it != x.end(); ++it) {
        os << *it << (i == n-1 ? " }" : ", "); i++;
    }
    if (x.empty()) os << "}";
    return os;
}
