//
// Created by giacomo on 21/02/20.
//

#ifndef HIERARCHY_TESTS_RESULTMAP_H
#define HIERARCHY_TESTS_RESULTMAP_H


struct result_map {
    double path_length_size = 0;
    double spearman = 0;
    double precision_leqK = 0;
    double precision_narrow = 0;
    double ncdg = 0;
    double recall_gtK = 0;
    double smallerNotCandidate = 0;

    result_map& operator+=(const result_map& rhs) {
        path_length_size += rhs.path_length_size;
        spearman += rhs.spearman;
        precision_leqK += rhs.precision_leqK;
        precision_narrow += rhs.precision_narrow;
        ncdg += rhs.ncdg;
        recall_gtK += rhs.recall_gtK;
        smallerNotCandidate += rhs.smallerNotCandidate;
        return *this;
    }
};

#endif //HIERARCHY_TESTS_RESULTMAP_H
