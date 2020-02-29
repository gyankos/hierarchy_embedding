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

    /**
     *
     * @param current                           Leaf element from which we're going to assess the distance
     * @param currentRankMap                    Reference to the expected rank map
     * @param getVectorRepresentation           Lambda function returning the memoized representation of the data structure
     * @param similarity                        Lambda function for the similarity scorer
     */
    TestingGraphLambda(const size_t& current, std::map<size_t, size_t>& currentRankMap, std::function<ForComparison(size_t)>& getVectorRepresentation,
                       std::function<double(const ForComparison&,const ForComparison&)>& similarity) :
            current{current}, currentRankMap{currentRankMap}, getVectorRepresentation{getVectorRepresentation}, similarity{similarity} {
        vcurrent = getVectorRepresentation(current);
    }

    size_t map_size = 0;

    /**
     * Callback function that is called to evaluate the distance among all the elements within the data structure.
     * pollMap will contain all the score sorted, and then the
     *
     * @param first
     * @param unused_to_remove
     */
    void operator()(size_t first, size_t unused_to_remove) {
        auto vfirst = getVectorRepresentation(first);
        double score = similarity(vcurrent, vfirst);
        pollMap.add(score, first);
        map_size++;
    }

    /**
     *
     * @param maps          Map where to write the final result for the current thread
     * @param maxLength     Maximum length of the elements out of the boundaries of the K-th path
     */
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
        while (it != castor_et_pollux.rend()) { // Iterating over the current implementation's ranked elements from the smaller to the largest
            size_t N = it->second.size();

            overallX += N;
            for (const auto& x : it->second) {
                auto it2 = currentRankMap.find(x);
                if (inferredPollCurrentPos > kSplit) { // If I am already past the K elements, then I can only intercept the real positive values with a recall measure
                    recallNumber++;
                    if (it2 != currentRankMap.end()) {
                        if(candidates_recall.insert(x).second)
                            recallToBeZero+=1.0;
                    }
                } else {
                    precisionNumber++;
                    if (it2 != currentRankMap.end()) { // If I am below the K elements, I suppose that I can use this information to reconstruct the path
                        if(candidates.insert(x).second)
                            retrieved_for_precision+=1.0;
                    }
                }

                // Determining the rank that is provided in the data structure: if it was not there, then it means that it is not a valid candidate, and hence it will be mapped in the last ranked position
                double expectedPollCurrentPos = (it2 != currentRankMap.end()) ? ((double)it2->second) : (((double)maxLength));
                calculateY.emplace_back(expectedPollCurrentPos);
                calculateX.emplace_back(inferredPollCurrentPos);
                sumY += expectedPollCurrentPos;
                sumX += inferredPollCurrentPos;
            }
            inferredPollCurrentPos++;
            it++;
        }

        // Calculating precision and recall
        maps.precision_leqK += ((double)retrieved_for_precision) / (precisionNumber);
        maps.recall_gtK    += (recallNumber == 0.0) ? 0.0 : ((double)recallToBeZero) / (recallNumber);

        // Last, calculating correlation coefficient
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
