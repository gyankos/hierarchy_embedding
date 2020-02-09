//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_TESTINGTREEPROPOSAL_H
#define HIERARCHY_TESTS_TESTINGTREEPROPOSAL_H

#include <tests/TestingTree.h>
#include <proposal/Proposal.h>

class TestingTreeProposal : public TestingTree<std::vector<size_t >> {
    Proposal prop;
    std::vector<std::vector<size_t>> allPossiblePaths;
public:
    TestingTreeProposal(size_t maximumBranchingFactor, double distanceFactor, double decayFactor, size_t  maxHeight);
protected:
    void initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subgraph_as_paths) override;
    std::vector<size_t> getVectorRepresentation(const std::vector<size_t> &current) override;
    double similarity(const std::vector<size_t> &lhs, const std::vector<size_t> &rhs) override;
    void generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) override;
    void finalizeDataIngestion() override {}
};

/**
 * Old testing, for checking the correctness of the implementation, independently form the metric of choice, that
 * is the thing used in the upper representation
 *
 * @param maximumBranchingFactor
 * @param distanceFactor
 * @param decayFactor
 * @param ls
 * @return
 */
double test_my_implementation(size_t maximumBranchingFactor, size_t distanceFactor, size_t decayFactor, const std::vector<std::vector<size_t>> &ls);


void test_my_proposal();;

#endif //HIERARCHY_TESTS_TESTINGTREEPROPOSAL_H
