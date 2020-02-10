#include <iostream>

#include <tests/tree/TestingTreeLearning.h>
#include <tests/tree/TestingTreeProposal.h>
#include "tests/tree/TestingTreeExternal.h"


#include <tests/TestingGraph.h>
#include <tests/tree/TestingTreePar.h>

/**
 * This generic function computes all the operations in function for all the possibile parameters
 * @tparam T
 * @param maxBranch
 * @param maxHeight
 * @param function
 */
template <class T> void perform_test(size_t maxBranch, size_t maxHeight, T function) {
    for (size_t i =/* 2*/maxBranch; i<=maxBranch; i++) {
        for (size_t j = /*2*/maxHeight; j<=maxHeight; j++) {
            function(i, j);
        }
    }
}

void write_poincarre_files(size_t maximumBranchingFactor, size_t maximumHeight) {
    std::string filename{"usecase_" + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv"};
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
        std::map<std::string, std::set<std::string>> calculateClosure;

        for (const auto& x : ls) {
            for (std::vector<size_t> y : x) {
                std::string current = size_vector_to_string(y);
                y.pop_back();
                while (!y.empty()) {
                    std::string upper = size_vector_to_string(y);
                    calculateClosure[current].insert(upper);
                    y.pop_back();
                }
                calculateClosure[current].insert("e");
            }
        }

        myfile << "id1,id2,weight" << std::endl;
        for (const auto& x : calculateClosure) {
            for (const auto& y : x.second) {
                myfile << x.first << ',' << y << ",1" << std::endl;
            }
        }
        myfile.close();
    }
}

void test_proposal(size_t maximumBranchingFactor, size_t maximumHeight) {
    std::cerr << maximumBranchingFactor << ", " << maximumHeight << std::endl;
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    //size_t maximumBranchingFactor, double distanceFactor, double decayFactor, size_t  maxHeight
    TestingTreeProposal tepee{maximumBranchingFactor, 3, (double)2.0, maximumHeight};
    tepee.run(ls);
}

void external_testing(size_t maximumBranchingFactor, size_t maximumHeight) {
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    std::cout << maximumBranchingFactor << ", " << maximumHeight << " Poincarré@7" << std::endl;
    {
        TestingTreeExternal p{maximumBranchingFactor, maximumHeight, "/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/results/usecase_" + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv.7pth.txt"};
        p.run(ls);
    }
    std::cout << maximumBranchingFactor << ", " << maximumHeight << " Poincarré@50" << std::endl;
    {
        TestingTreeExternal p{maximumBranchingFactor, maximumHeight, "/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/results/usecase_" + std::to_string(maximumBranchingFactor) + "_" + std::to_string(maximumHeight) + ".csv.50pth.txt"};
        p.run(ls);
    }
}

void testing_learning(size_t maximumBranchingFactor, size_t maximumHeight) {
    std::cout << maximumBranchingFactor << ", " << maximumHeight << " EEEL@" << 7 << std::endl;
    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    TestingTreeLearning testing_dimension_as_branching{maximumBranchingFactor, maximumHeight, 7};
    testing_dimension_as_branching.run(ls);

    std::cout << maximumBranchingFactor << ", " << maximumHeight << " EEEL@" << 50 << std::endl;
    TestingTreeLearning testing_dimension_mid_branching_100{maximumBranchingFactor, maximumHeight,
                                                            50};
    testing_dimension_mid_branching_100.run(ls);
}

void parhier_testing(size_t maximumBranchingFactor, size_t maximumHeight) {
    double distanceFactor = 3;
    double decayFactor = 2;

    auto ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    {
        std::cout << maximumBranchingFactor << ", " << maximumHeight << " ParHierCluster k=3 @" << std::endl;
        TestingTreePar testing{maximumBranchingFactor, maximumHeight, 3, CLUSTER, 3.0, 2.0};
        testing.run(ls);
    }
    {
        std::cout << maximumBranchingFactor << ", " << maximumHeight << " ParHierSum k=3 @" << std::endl;
        TestingTreePar testing{maximumBranchingFactor, maximumHeight, 3, SUM, 3.0, 2.0};
        testing.run(ls);
    }
}

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <JLLemma.h>
#include <tests/tree/TestingTreeBasic.h>
#include "Graph.h"
#include "tests/graph/TestingGraphProposal.h"

int main() {


#if 0
    // Correspondence between the node in adj and the node in the graph, so that we can reconstruct all the vectors to be associated to it
    std::unordered_map<size_t, size_t> treeToGraphMorphism;
    // Tree representation of the data structure, with new element ids
    std::unordered_map<size_t, std::unordered_set<size_t>> tree;
    // path string generated from the adj representation. This can be used to then generate multiple possible elements
    std::map<size_t, std::string> treeIdToPathString;

    // Filling up all the three data structures provided up above.
    g.generateNaryTree(treeToGraphMorphism, tree, treeIdToPathString);
#endif

#if 1
    std::ifstream file{"noun.txt"};
    Graph g{file};
    file.close();
    double distanceFactor = 3;
    double decayFactor = 2;
    TestingGraphProposal proposal(distanceFactor, decayFactor, g);
    proposal.run(false);
#endif

    ///perform_test(7, 5, parhier_testing);

    /*std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(7, 5);
    TestingProposal testing_dimension_as_branching{7, 3, 2, 5};
    testing_dimension_as_branching.run(ls);*/

    //perform_test(7, 5, external_testing);
    //size_t maximumBranchingFactor, size_t maximumHeight,size_t k, method m, double distanceFactor, double decayFactor
   /* ParTesting testing{6, 2, 3, CLUSTER, 3.0, 2.0};
    auto ls = generateCompleteSubgraph(6, 2);
    testing.run(ls);
*/
   // perform_test(5, 5, external_testing);
    /*ExternalTesting p{5, 4, "/media/giacomo/Data/hierarchy_paper/projects/poincare-embeddings/results/usecase_5_4.csv.7pth.txt"};
    std::vector<std::vector<std::vector<size_t>> > ls = generateCompleteSubgraph(5, 4);
    p.run(ls);*/

    /*
    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;
    std::cout << "Generating the complete tree" << std::endl;

    perform_test(7, 5, write_poincarre_files);*/

    //std::cout << generateAllPossibleSubpaths({1,2,3}) << std::endl;
    //std::cout << test_basic_4(ls) << std::endl;
    //std::cout << test_my_implementation(maximumBranchingFactor, distanceFactor, decayFactor, ls) << std::endl;
    return 0;
}
