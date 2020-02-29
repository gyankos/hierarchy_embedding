//
// Created by giacomo on 29/12/19.
//

#include <algorithm>
#include <multithreaded/thread_pool.h>
#include <multithreaded/MultithreadWrap.h>
#include "proposal/proposal_utils.h"

int indexOfSubList(std::vector<size_t> &array, std::vector<size_t> &subarray) {
    std::vector<size_t>::iterator pos = std::search(
            array.begin(), array.end(),
            subarray.begin(), subarray.end());

    if(pos == array.end())
        return -1;
            else
    // found at pos
        return pos - array.begin();
    /*size_t sourceSize = array.size();
    size_t targetSize = subarray.size();
    if (targetSize > sourceSize)
        return -1;
    size_t maxCandidate = sourceSize - targetSize;

    size_t candidate = 0;
    while ( candidate <= maxCandidate) {
        bool doContinue = false;
        for (size_t i=0, j=candidate; i<targetSize; i++, j++)
            if (subarray[i] != array[j]) {
                candidate++;
                doContinue = true;
                continue;// nextCand;  // Element mismatch, try next cand
            }
        if (!doContinue) return candidate;  // All elements of candidate matched target
    }
    return -1;  // No candidate matched the target*/
}

void
generateCompleteSubgraph(std::vector<std::vector<size_t>> &finalList, size_t currentHeight, size_t maximumBranching,
                         std::vector<size_t> &currentArray) {
    if (currentHeight > 0) {
        for (int i = 1; i<=maximumBranching; i++) {
            std::vector<size_t> current = currentArray;
            current.emplace_back(i);
            finalList.emplace_back(current);
            generateCompleteSubgraph(finalList, currentHeight-1, maximumBranching, current);
        }
    }
}

bool subArrayOf(std::vector<size_t> &array, std::vector<size_t> &subarray) {
    if (array.size() < subarray.size()) {
        return subArrayOf(subarray, array);
    }
    return indexOfSubList(array, subarray) == 0;
}



std::vector<std::vector<std::vector<size_t>> > generateCompleteSubgraph(size_t maximumBranchingFactor, size_t maximumHeight) {
    std::vector<std::future<std::vector<std::vector<size_t>> >> futures;
    unsigned int nthreads = std::thread::hardware_concurrency() / 2;
    MultithreadWrap<std::vector<std::vector<size_t>>> wrapper{nthreads, IS_MULTITHREADED};
    //thread_pool pool(nthreads);
    if (maximumBranchingFactor > 0) {
        for (size_t i = 1; i<=maximumBranchingFactor; i++) {
            /*futures.push_back(pool.execute([maximumBranchingFactor, maximumHeight](size_t i) {
                std::vector<std::vector<size_t>>  result;
                std::vector<size_t> current;
                current.emplace_back(i);
                result.emplace_back(current);
                generateCompleteSubgraph(result, maximumHeight-1, maximumBranchingFactor, current);
                return result;
                }, i));*/
            wrapper.poolExecute([maximumBranchingFactor, maximumHeight](size_t i) {
                std::vector<std::vector<size_t>>  result;
                std::vector<size_t> current;
                current.emplace_back(i);
                result.emplace_back(current);
                generateCompleteSubgraph(result, maximumHeight-1, maximumBranchingFactor, current);
                return result;
            }, i);
        }
    }
    return wrapper.foreach();
}

double euclideanDistance(const std::vector<double> &left, const std::vector<double> &right) {
    double sum = 0.0;
    int M = left.size();
    int N = right.size();
    for (int i = 0, n = std::max(M, N); i<n; i++) {
        sum += std::pow((i<M ? left[i] : 0.0)- (i<N ? right[i] : 0.0), 2);
    }
    return std::sqrt(sum);
}

double poincarreDistance(const std::vector<double> &left, const std::vector<double> &right) {
    double num = std::pow(euclideanDistance(left, right), 2);
    double den = (1-vectorNorm2(left.begin(), left.end())) * (1-vectorNorm2(right.begin(), right.end()));
    return std::acosh(1+2 * (num/den));
}

