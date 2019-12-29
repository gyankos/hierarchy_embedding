//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_TESTING_H
#define HIERARCHY_TESTS_TESTING_H

#include <vector>
#include <string>
#include <proposal/proposal_utils.h>
#include <string_utils.h>
#include <map>
#include <set>
#include <iostream>
#include "PollMap.h"
#include "stats_utils.h"
#include "thread_pool.h"

#define         DEBUG           (false)

struct result_map {
    double path_length_size;
    double spearman;
    double precision;
    double precision_narrow;
    double ncdg;
    double recall;
    double smallerNotCandidate;
};

template <typename ForComparison>
class Testing {
    size_t maximumBranchingFactor, maximumHeight;

    double Spearman(const std::vector<size_t>& current, const std::vector<std::vector<size_t>> &candidates) {
        std::set<std::string> candidatesSet;
        double noCandidates = candidates.size();
        double retrieved = 0.0;
        std::string currentString = size_vector_to_string(current);
        std::map<std::string, size_t> toretMap, currentRankMap;

        for (auto& y : candidates) {
            std::string tmp = size_vector_to_string(y);
            candidatesSet.emplace(size_vector_to_string(y));
            currentRankMap[tmp] = candidates.size() - y.size() + 1;
        }

        //if (DEBUG) std::cout << "Candidates: " << candidatesSet << " for element " << current << std::endl;
        PollMap<double,std::string> pollMap{candidatesSet.size()};
        generateTopKCandidates(pollMap, current);
        pollMap.getRankedPoll(toretMap);

        for (std::map<std::string, size_t>::iterator it = toretMap.begin(); it != toretMap.end(); it++) {
            if (candidatesSet.find(it->first) == candidatesSet.end()) {
                currentRankMap[it->first] = candidates.size() + 2;
            }
        }

        return Spearmam_rank_correlation_coefficient<std::string>(toretMap, currentRankMap);
    }

    double NDCG(const std::vector<size_t>& current, const std::vector<std::vector<size_t>> &candidates) {
        std::set<std::string> candidatesSet;
        double noCandidates = candidates.size();
        double retrieved = 0.0;
        std::string currentString = size_vector_to_string(current);
        std::map<std::string, size_t> toretMap, currentRankMap;

        for (auto& y : candidates) {
            std::string tmp = size_vector_to_string(y);
            candidatesSet.emplace(size_vector_to_string(y));
            currentRankMap[tmp] =  y.size();
        }

        //if (DEBUG) std::cout << "Candidates: " << candidatesSet << " for element " << current << std::endl;
        PollMap<double,std::string> pollMap{candidatesSet.size()};
        generateTopKCandidates(pollMap, current);
        const std::map<double, std::set<std::string>> &poll = pollMap.getPoll();

        return normalized_discounted_cumulative_gain<std::string>(poll, currentRankMap);
    }

    double Precision(const std::vector<size_t>& current, const std::vector<std::vector<size_t>> &candidates) {
        std::set<std::string> candidatesSet;
        double retrieved = 0.0;
        std::string currentString = size_vector_to_string(current);
        std::map<std::string, size_t> toretMap, currentRankMap;

        for (auto& y : candidates) {
            std::string tmp = size_vector_to_string(y);
            candidatesSet.emplace(size_vector_to_string(y));
            currentRankMap[tmp] = candidates.size() - y.size() + 1;
        }

        //if (DEBUG) std::cout << "Candidates: " << candidatesSet << " for element " << current << std::endl;
        PollMap<double,std::string> pollMap{candidatesSet.size()};
        generateTopKCandidates(pollMap, current);
        std::multimap<double, std::string> mmap;
        pollMap.getPoll(mmap);
        double noCandidates = mmap.size();
        for (auto it = mmap.rbegin(); it != mmap.rend(); it++) {
            if (candidatesSet.find(it->second) != candidatesSet.end()) {
                retrieved+=1.0;
            }
        }
        return retrieved / noCandidates;
    }

    double Precision_narrow(const std::vector<size_t>& current, const std::vector<std::vector<size_t>> &candidates) {
        std::set<std::string> candidatesSet;
        double retrieved = 0.0;
        std::string currentString = size_vector_to_string(current);
        std::map<std::string, size_t> toretMap, currentRankMap;

        for (auto& y : candidates) {
            std::string tmp = size_vector_to_string(y);
            candidatesSet.emplace(size_vector_to_string(y));
            currentRankMap[tmp] = candidates.size() - y.size() + 1;
        }

        //if (DEBUG) std::cout << "Candidates: " << candidatesSet << " for element " << current << std::endl;
        double currentlyGot = 0.0;
        PollMap<double,std::string> pollMap{candidatesSet.size()};
        generateTopKCandidates(pollMap, current);
        pollMap.getNarrowRankedPoll(toretMap);

        for (auto& cp : toretMap) {
            if (candidatesSet.find(cp.first) != candidatesSet.end()) {
                //std::cout << "<" << cp.first << ", " << cp.second << ">" << std::endl;
                retrieved += 1.0;
            }
            currentlyGot += 1.0;
            //if (DEBUG) std::cout << "<" << cp.first << ", " << cp.second << ">" << std::endl;
        }
        //std::cout << std::endl;
        return retrieved / currentlyGot;
    }

public:
    Testing(size_t maximumBranchingFactor, size_t maximumHeight) : maximumBranchingFactor(maximumBranchingFactor),
                                                                   maximumHeight(maximumHeight) {}


    std::pair<double,double> minmaxCandidates(const std::vector<size_t>& current, const std::vector<std::vector<std::vector<size_t>>>& z, const std::vector<std::vector<size_t>> &candidates) {
        double min = std::numeric_limits<double>::max();
        std::vector<std::vector<size_t >> minCandidates;
        double max = std::numeric_limits<double>::min();
        std::vector<std::vector<size_t >> maxCandidates;

        // Getting all the elements required for Precision
        std::set<std::string> candidatesSet;
        std::string currentString = size_vector_to_string(current);
        for (auto& y : candidates) {
            std::string tmp = size_vector_to_string(y);
            candidatesSet.emplace(size_vector_to_string(y));
        }

        std::string currentStirng = size_vector_to_string(current);
        for (auto& candidates : z) {
            for (auto &x : candidates) {
                std::string x_tmp = size_vector_to_string(x);
                if (x_tmp != currentStirng) {
                    double z = similarity(getVectorRepresentation(current), getVectorRepresentation(x));
                    if (z < min) {
                        minCandidates.clear();
                        minCandidates.emplace_back(x);
                        min = z;
                    } else if (z == min) {
                        minCandidates.emplace_back(x);
                    }
                    if (z > max) {
                        maxCandidates.clear();
                        maxCandidates.emplace_back(x);
                    } else if (z == max) {
                        maxCandidates.emplace_back(x);
                    }
                    max = std::max(max, z);
                }
            }
        }

        double recallToBeZero = 0.0;
        for (const std::vector<size_t>& x : minCandidates) {
            std::string current = size_vector_to_string(x);
            if (candidatesSet.find(current) != candidatesSet.end())
                recallToBeZero += 1.0;
        }

        return std::make_pair(recallToBeZero/candidatesSet.size(), 1.0/(1.0+maxCandidates.size()));
    }



    void run(std::vector<std::vector<std::vector<size_t>>>& ls) {
        // Generating all the possible nodes for the hierarchy of choice
        //const std::vector<std::vector<size_t>> &ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
        std::cout << "Tree Generation" << std::endl;
        for (auto& element : ls)
            initialize_hierarchy_with_all_paths(element);

        double path_length_size = 0, spearman = 0,  precision = 0,  precision_narrow = 0, ncdg = 0, smallerNotCandidate = 0, recall =0;

        std::vector<std::future<struct result_map >> futures;
        thread_pool pool(ls.size());

        std::cout << "Now Computing..." << std::endl;
        // Performing the test for each node in the hierarchy, root excluded
        for (auto& y: ls) {
            futures.push_back(pool.execute([this, ls](const std::vector<std::vector<size_t>>& y) {
                struct result_map maps;

                for (auto& x : y) {
                    size_t current_path_length = x.size();

                    if (current_path_length != maximumHeight) continue;

                    // Increment the path length size by one and, if not present, insert 1.0
                    maps.path_length_size += 1.0;

                    // Generating all the candidates, represented as indices, as relevant is-a nodes for the current candidate.
                    const std::vector<std::vector<size_t>> &candidates = generateAllPossibleSubpaths(x);
                    std::pair<double,double> minmax = this->minmaxCandidates(x, ls, candidates);

                    // Testing the precision at k representation
                    //double currentPrecision = this->Spearman(x, candidates);
                    //if (DEBUG)  std::cout << "Current precision of " << currentPrecision << " for length " << current_path_length << std::endl;
                    maps.spearman += this->Spearman(x, candidates); // Getting the classification difference between the two
                    maps.precision += this->Precision(x, candidates);
                    maps.precision_narrow += this->Precision_narrow(x, candidates);
                    maps.ncdg += this->NDCG(x, candidates);
                    maps.recall += minmax.first;
                    maps.smallerNotCandidate += minmax.second;
                }

                return maps;
            }, y));
        }



        std::cout << "Summing up the maps..." << std::endl;
        for (auto& y : futures) {
            auto x = y.get();
#define update_local(field, currFuture)           (field) += currFuture .  field
            update_local(path_length_size, x);
            update_local(spearman, x);
            update_local(precision, x);
            update_local(precision_narrow, x);
            update_local(ncdg, x);
            update_local(recall, x);
            update_local(smallerNotCandidate, x);
        }
        futures.clear();



        std::cout << "Spearman, for top-k                                              " << spearman / path_length_size << std::endl;
        //print_result_maps(path_length_size, spearman);
        std::cout << "Precision, for top-k scores                                      " << precision / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision);
        std::cout << "Precision, for top-k elements                                    " << precision_narrow / path_length_size << std::endl;
        std::cout << "Recall, for non top-k elements (wrongly matching the candidates) " << recall / path_length_size << std::endl;
        std::cout << "More similar than worst top-1 not candidate, not current element " << smallerNotCandidate / path_length_size << std::endl;
        //print_result_maps(path_length_size, precision_narrow);
    }

    void
    print_result_maps(const std::map<size_t, double> &keys_with_sizes,
                      const std::map<size_t, double> &keys_with_values) const {
        for(auto it_m1 = keys_with_sizes.cbegin(), end_m1 = keys_with_sizes.cend(),
                    it_m2 = keys_with_values.cbegin(), end_m2 = keys_with_values.cend();
            it_m1 != end_m1 &&  it_m2 != end_m2; it_m1++, it_m2++)
        {
            std::cout << it_m1->first << " = " << it_m2->second / it_m1->second << std::endl;
        }
    }

protected:
    virtual void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) = 0;
    virtual ForComparison getVectorRepresentation(const std::vector<size_t>& current) = 0;
    virtual double similarity(const ForComparison& lhs, const ForComparison& rhs) = 0;
    virtual void  generateTopKCandidates(PollMap<double,std::string>& map, const std::vector<size_t>& current) = 0;
};


#endif //HIERARCHY_TESTS_TESTING_H
