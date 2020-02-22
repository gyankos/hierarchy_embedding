//
// Created by giacomo on 21/02/20.
//

#ifndef HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H
#define HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H

#include <functional>
#include <map>
#include "ResultMap.h"

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

    /*double min_minmax = std::numeric_limits<double>::max();
    std::vector<size_t> minCandidates_minmax;
    double max_minmax = std::numeric_limits<double>::min();
    std::vector<size_t> maxCandidates_minmax;*/

    //std::set<size_t>    precision;
    size_t              retrieved_for_precision = 0;
    size_t              noCandidates;


    double              expectedMeanPollRank = 0.0;

    void operator()(size_t first, size_t second) {
        //std::cout << first << "," << second << std::endl;
        auto vfirst = getVectorRepresentation(first);
        double score = similarity(vcurrent,vfirst);
        //if (score == 0) return; // Zero score means that the current element is not a candidate

        // Spearman
        expectedMeanPollRank += second;

        // MinMax computation
        //if (first != current) {
        /// ???????
#if 0
        if (score < min_minmax) {
            minCandidates_minmax.clear();
            minCandidates_minmax.emplace_back(first);
            min_minmax = score;
        } else if (score == min_minmax) {
            minCandidates_minmax.emplace_back(first);
        }
        if (score > max_minmax) {
            maxCandidates_minmax.clear();
            maxCandidates_minmax.emplace_back(first);
        } else if (score == max_minmax) {
            maxCandidates_minmax.emplace_back(first);
        }
        max_minmax = std::max(max_minmax, score);
#endif
        //}
        // MiniMax computation

        // Globals
        pollMap.add(score, first);
        map_size++;
    }

    void finalize(struct result_map& maps, int maxLength) {
        maxLength = std::max(maxLength, 0);
        size_t kSplit = currentRankMap.size();
        double recallToBeZero = 0.0;
        //pollMap.resize(kSplit);

#if 0
        // MiniMax computation

        maps.recall += recallToBeZero / minCandidates_minmax.size();
#endif
        //maps.smallerNotCandidate +=  maxCandidates_minmax.size() / (1.0 + maxCandidates_minmax.size());
        // MiniMax computation


        // Globals: generating the candidates that are just the ones returned by the algorithm itself
        maps.path_length_size += 1.0;
        //PollMap<double,size_t> pollMap{map_size};
        //generateTopKCandidates(pollMap, current);
        const std::map<double, std::set<size_t>> &castor_et_pollux = pollMap.getPoll();

        // Spearman + NCDG
        double NDCG_pos = 1.0;
        double NDCG_dcg_p = 0.0;
        double NDCG_idcg_p = 0.0;

        expectedMeanPollRank /= ((double)map_size);

        size_t inferredPollCurrentPos = 1;
        double inferredMeanPollRank = 0.0;
        size_t overallX = 0;
        std::vector<double> precalculateX;
        std::vector<double>    calculateY;
        double                 calculateYSquaredSummed = 0.0;
        auto it = castor_et_pollux.rbegin();

        std::unordered_set<size_t>    candidates;
        std::unordered_set<size_t>    candidates_recall;
        double overK = 0;
        while (it != castor_et_pollux.rend()) {
            size_t N = it->second.size();
            //bool doNCDG = false/*(N > 0)*/;
            //double div;
            //double toSum = 0.0;
            //double reli;

            overallX += N;
            for (const auto& x : it->second) {
                auto it2 = currentRankMap.find(x);
                if (inferredPollCurrentPos > kSplit) {
                    overK += N;
                    if (it2 != currentRankMap.end()) {
                        if(candidates_recall.insert(x).second)
                            recallToBeZero+=1.0;
                    }
                } else {
                    if (it2 != currentRankMap.end()) {
                        if(candidates.insert(x).second)
                            retrieved_for_precision+=1.0;
                    }
                }

                double tmpY = (it2 != currentRankMap.end()) ? ((double)it2->second)-expectedMeanPollRank : (((double)maxLength) - expectedMeanPollRank);
                calculateYSquaredSummed += std::pow(tmpY, 2.0);
                calculateY.emplace_back(tmpY);
                precalculateX.emplace_back(inferredPollCurrentPos);
                //toSum += (double)currentRankMap[x];
            }
            for (size_t i = 0; i<N; i++)
                inferredMeanPollRank += (inferredPollCurrentPos * N);
            inferredPollCurrentPos++;
            it++;
        }

        /*for (const auto& it : castor_et_pollux) {
            size_t N = it.second.size();
            //bool doNCDG = (N > 0);
            //double div;
            //double toSum = 0.0;
            //double reli;

            overallX += N;
            for (const auto& x : it.second) {
                if (currentRankMap.find(x) != currentRankMap.end()) {
                    if(candidates.insert(x).second)
                        retrieved_for_precision+=1.0;
                }

                double tmpY = ((double)currentRankMap[x])-expectedMeanPollRank;
                calculateYSquaredSummed += std::pow(tmpY, 2.0);
                calculateY.emplace_back(tmpY);
                precalculateX.emplace_back(inferredPollCurrentPos);
                //toSum += (double)currentRankMap[x];
            }
            for (size_t i = 0; i<N; i++)
                inferredMeanPollRank += (inferredPollCurrentPos * N);
            inferredPollCurrentPos++;

        }*/
        inferredMeanPollRank /= ((double)overallX);

        /// Precision
        maps.precision += ((double)retrieved_for_precision)/((double)currentRankMap.size());
        //std::cout << maps.precision << std::endl;
        maps.recall    += (overK == 0.0) ? 0.0 : ((double)recallToBeZero)/(overK);

        double numeratore = 0.0;
        double                 calculateXSquaredSummed = 0.0;
        for (size_t i = 0; i<overallX; i++) {
            double tmpX = (precalculateX[i] - inferredMeanPollRank);
            calculateXSquaredSummed += std::pow(tmpX, 2.0);
            numeratore += tmpX * calculateY[i];
        }

        maps.spearman += (overallX > 1) ? numeratore / std::sqrt(calculateXSquaredSummed * calculateYSquaredSummed) : 1.0-abs(inferredMeanPollRank-expectedMeanPollRank)/(1+abs(inferredMeanPollRank-expectedMeanPollRank));
    }

};



#endif //HIERARCHY_TESTS_TESTINGGRAPHLAMBDA_H
