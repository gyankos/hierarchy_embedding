#include <iostream>
#include <vector>
#include <proposal/Proposal.h>
#include "cout_utils.h"
#include "learning/EEE_Engine.h"
#include <concept_vector/TestingBasic1.h>
#include <concept_vector/TestingBasic2.h>
#include <concept_vector/TestingBasic3.h>
#include <fixed_bimap.h>
#include <learning/Blob.h>


double test_my_implementation(size_t maximumBranchingFactor, size_t distanceFactor, size_t decayFactor, const std::vector<std::vector<size_t>> &ls) {
    Proposal hd(maximumBranchingFactor, distanceFactor, decayFactor);
    std::cout << "Starting the testing..." << std::endl;
    double tests = 0.0;
    double correct = 0.0;
    for (int i = 0, n = ls.size(); i<n; i++) {
        std::cout << i << std::endl;
        auto veci = ls[i];
        for (int j = 0; j<i; j++) {
            auto vecj = ls[j];
            bool consistenceInference = hd.isConsistent(veci, vecj);
            tests++;
            if ((subArrayOf(veci, vecj)) || (subArrayOf(vecj, veci))) {
                if (!consistenceInference)
                    std::cerr << ("ERROR! the algorithm says that is not consistent, while it should be") << std::endl;
                else
                    correct++;
            } else {
                if (consistenceInference) {
                    std::cerr << veci << " vs. " << vecj << std::endl;
                    bool consistenceInference = hd.isConsistent(veci, vecj);
                    std::cerr << ("ERROR! the algorithm says that is consistent, while it is not") << std::endl;
                } else
                    correct++;
            }
        }
    }
    return (correct/tests);
}




void testing_basic_implementation() {
    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;

    size_t maximumHeight = 4;
    std::cout << "Generating the complete tree" << std::endl;
    std::vector<std::future<std::vector<std::vector<size_t>>>> tmp = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    std::vector<std::vector<std::vector<size_t>>> ls;
    for (auto& element : tmp)
        ls.emplace_back(element.get());

    TestingBasic1 testingBasic1{maximumBranchingFactor, maximumHeight};
    testingBasic1.run(ls);

    TestingBasic2 testingBasic2{maximumBranchingFactor, maximumHeight, 0.75};
    testingBasic2.run(ls);

    TestingBasic3 testingBasic3{maximumBranchingFactor, maximumHeight, 0.75, 0.5};
    testingBasic3.run(ls);

    TestingBasic3 testingBasic4{maximumBranchingFactor, maximumHeight, 1, 0.5};
    testingBasic4.run(ls);
};


int main() {


    size_t maximumBranchingFactor = 5;
    double distanceFactor = 3;
    int decayFactor = 2;
    size_t maximumHeight = 3;

    std::cout << "Generating the complete tree" << std::endl;
    naryTree tree{0};
    size_t id = 0, num_entities_and_classes = 1; //num_entities_and_classes is initialized with 1, because the root will be added automatically

    std::cout << "Generating the positive examples" << std::endl;
    std::map<std::string, std::set<std::string>> positiveExamples;
    fixed_bimap<std::string, size_t>             bimap;             // Bijection between path as a string and the id
    std::set<std::string> allNodes;                                 // Containing all the nodes within the hierarchy

    // Adding the root node
    bimap.put("", 0);
    for (auto& element : generateCompleteSubgraph(maximumBranchingFactor, maximumHeight)) {
        for (auto& x : element.get()) {
            std::string x_val = size_vector_to_string(x);
            allNodes.emplace(x_val);
            for (auto& y : generateAllPossibleSubpaths(x)) {
                std::string y_val = size_vector_to_string(y);
                positiveExamples[x_val].emplace(y_val);
                positiveExamples[y_val].emplace(x_val);
            }
            // Getting the id of the new element that is now added: id
            num_entities_and_classes += tree.addChild(x, id);
            // Adding the id of the associated element, and associating that to the string. That is going ot be used to determine the best common ancestor.
            bimap.put(x_val, id);
        }
    }

    // Init category hierarchy
    EEEngine engine{EDIAG, maximumBranchingFactor, tree, bimap};
    engine.entity_category_hierarchy_.InitHierarchy(num_entities_and_classes);
    engine.InitEntityCategories(num_entities_and_classes);
    engine.ReadHierarchyIdWithLevel(tree);

    std::cout << "Generating the negative examples" << std::endl;
    std::map<std::string, std::set<std::string>> negativeExamples;
    for (const std::string& x : allNodes) {
        const std::set<std::string>& set = positiveExamples[x];
        std::set<std::string>& difference = negativeExamples[x];
        std::set_difference(set.begin(), set.end(), allNodes.begin(), allNodes.end(), std::inserter(difference, difference.begin()));
    }

    // Finished to generate the

    //std::cout << generateAllPossibleSubpaths({1,2,3}) << std::endl;
    //std::cout << test_basic_4(ls) << std::endl;
    //std::cout << test_my_implementation(maximumBranchingFactor, distanceFactor, decayFactor, ls) << std::endl;
    return 0;
}