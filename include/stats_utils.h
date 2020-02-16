//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_STATS_UTILS_H
#define HIERARCHY_TESTS_STATS_UTILS_H

#include <map>

template <typename K> double Spearmam_rank_correlation_coefficient(const std::map<K,size_t>& a,                                const std::map<K,size_t>& b ) {
    double dSum = 0.0;
    for(auto it_m1 = a.cbegin(), end_m1 = a.cend(),
                it_m2 = b.cbegin(), end_m2 = b.cend();
        it_m1 != end_m1 &&  it_m2 != end_m2; it_m1++, it_m2++)
    {
        dSum += (double)std::abs((double)it_m1->second - (double)it_m2->second);
//        std::cout << it_m1->first << " = " << it_m2->second / it_m1->second << std::endl;
    }
    return 1.0- (6*dSum/(a.size()*(a.size()*a.size()-1)));
}

template <typename K> double normalized_discounted_cumulative_gain(const std::map<double, std::set<K>> &poll, std::map<K,size_t>& b ) {
    double NDCG_pos = 1.0;
    double NDCG_dcg_p = 0.0;
    double NDCG_idcg_p = 0.0;
    for (auto it = poll.rbegin(); it != poll.rend(); it++) {
        if (it->second.empty()) continue;
        double div = log2(++NDCG_pos);
        // average rel
        double toSum = 0.0;
        for (auto k : it->second)
            toSum += (double)b[k];
        double reli = toSum / (double)it->second.size();     // Average ranking
        toSum = (1.0 + reli) / div;                           // Actual formula
        NDCG_dcg_p += toSum;
        NDCG_idcg_p += 1.0 + (std::pow(2, reli) - 1.0) / div;
    }
    return NDCG_dcg_p / NDCG_idcg_p;
}

#endif //HIERARCHY_TESTS_STATS_UTILS_H
