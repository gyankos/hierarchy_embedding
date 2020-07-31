//
// Created by giacomo on 21/02/20.
//

#ifndef HIERARCHY_TESTS_RESULTMAP_H
#define HIERARCHY_TESTS_RESULTMAP_H


#include <ostream>

struct result_map {
    double path_length_size = 0;
    double spearman = 0;
    double precision_leqK = 0;
    double recall_gtK = 0;
    double milliseconds = 0;

    result_map() = default;
    result_map(const struct result_map& x) = default;
    result_map& operator=(const struct result_map& x) = default;
    result_map(result_map&& x) = default;

    result_map& operator+=(const result_map& rhs) {
        path_length_size += rhs.path_length_size;
        spearman += rhs.spearman;
        precision_leqK += rhs.precision_leqK;
        recall_gtK += rhs.recall_gtK;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const result_map &map) {
        os << map.spearman << " (spe)," << map.precision_leqK << " (prec)," << map.recall_gtK << " (rec)," << map.milliseconds;
        return os;
    }
};

#endif //HIERARCHY_TESTS_RESULTMAP_H
