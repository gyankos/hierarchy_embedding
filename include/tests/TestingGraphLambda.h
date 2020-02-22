//
// Created by giacomo on 21/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H
#define HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H

#include <functional>
#include <map>
#include "ResultMap.h"

/**
 * Class simplifying the task of ranking the goodness of a solution
 * @tparam ForComparison
 */
template <typename ForComparison>
struct TestingGraphLambda {
    const size_t& current;
    std::function<ForComparison(size_t)>& getVectorRepresentation;
    std::function<double(const ForComparison&,const ForComparison&)>& similarity;
    std::map<size_t, size_t>& currentRankMap;

    PollMap<double, size_t> pollMap;
    ForComparison vcurrent;

    TestingGraphLambda(const size_t& current, std::map<size_t, size_t>& currentRankMap, std::function<ForComparison(size_t)>& getVectorRepresentation,
                       std::function<double(const ForComparison&,const ForComparison&)>& similarity) :
            current{current}, currentRankMap{currentRankMap}, getVectorRepresentation{getVectorRepresentation}, similarity{similarity} {
        vcurrent = getVectorRepresentation(current);
    }

    size_t map_size = 0;
    size_t              noCandidates;



    void operator()(size_t first, size_t second) {
        ///std::cout << "(1) " << first << "," << second << std::endl;
        auto vfirst = getVectorRepresentation(first);
        double score = similarity(vcurrent,vfirst);
        //if (score == 0) return; // Zero score means that the current element is not a candidate

        // Spearman
        //expectedMeanPollRank += second;

        pollMap.add(score, first);
        map_size++;
    }

    void finalize(struct result_map& maps, int maxLength) {
        maxLength = std::max(maxLength, 0);

        // Determining the k for the precision at <=k and recall >k
        size_t kSplit = 0;
        for (const auto& it : currentRankMap) {
            kSplit = std::max(kSplit, it.second);
        }

        // Globals: generating the candidates that are just the ones returned by the algorithm itself
        maps.path_length_size += 1.0;
        const std::map<double, std::set<size_t>> &castor_et_pollux = pollMap.getPoll();

        double sumX = 0, sumY = 0;

        size_t inferredPollCurrentPos = 1;
        size_t overallX = 0;
        std::vector<double> calculateX;
        std::vector<double>    calculateY;
        auto it = castor_et_pollux.rbegin();

        std::unordered_set<size_t>    candidates;
        std::unordered_set<size_t>    candidates_recall;

        double retrieved_for_precision = 0, recallToBeZero = 0.0;
        double precisionNumber = 0, recallNumber = 0;

        double overK = 0;
        while (it != castor_et_pollux.rend()) {
            size_t N = it->second.size();

            overallX += N;
            for (const auto& x : it->second) {
                auto it2 = currentRankMap.find(x);
                if (inferredPollCurrentPos > kSplit) {
                    recallNumber++;
                    if (it2 != currentRankMap.end()) {
                        if(candidates_recall.insert(x).second)
                            recallToBeZero+=1.0;
                    }
                } else {
                    precisionNumber++;
                    if (it2 != currentRankMap.end()) {
                        if(candidates.insert(x).second)
                            retrieved_for_precision+=1.0;
                    }
                }

                double expectedPollCurrentPos = (it2 != currentRankMap.end()) ? ((double)it2->second)/*-expectedMeanPollRank*/ : (((double)maxLength) /*- expectedMeanPollRank*/);
                calculateY.emplace_back(expectedPollCurrentPos);
                calculateX.emplace_back(inferredPollCurrentPos);

                sumY += expectedPollCurrentPos;
                sumX += inferredPollCurrentPos;
            }
            inferredPollCurrentPos++;
            it++;
        }

        maps.precision_leqK += ((double)retrieved_for_precision) / (precisionNumber);
        maps.recall_gtK    += (recallNumber == 0.0) ? 0.0 : ((double)recallToBeZero) / (recallNumber);

        // Spearman
        double numeratore = 0.0;
        double                 calculateXSquaredSummed = 0.0;
        double                 calculateYSquaredSummed = 0.0;
        sumY /= ((double)overallX);
        sumX /= ((double)overallX);
        for (size_t i = 0; i<overallX; i++) {
            double tmpX = (((double) calculateX[i]) - sumX);
            double tmpY = (((double) calculateY[i]) - sumY);
            calculateXSquaredSummed += std::pow(tmpX, 2.0);
            calculateYSquaredSummed += std::pow(tmpX, 2.0);
            numeratore += (tmpX) * (tmpY);
        }
        maps.spearman += (numeratore > 0) ? numeratore / std::sqrt(calculateXSquaredSummed * calculateYSquaredSummed) : 0;
    }

};



#endif //HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H
