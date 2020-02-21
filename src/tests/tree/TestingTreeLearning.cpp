//
// Created by giacomo on 30/12/19.
//

#include "tests/tree/TestingTreeLearning.h"

TestingTreeLearning::TestingTreeLearning(size_t maximumBranchingFactor, size_t maximumHeight, size_t vectorDimension)
        : TestingTree(maximumBranchingFactor,
                      maximumHeight), dimEmbedding(vectorDimension <= 0 ? std::max((size_t)2, maximumBranchingFactor) : vectorDimension) {
    ee_engine_ = nullptr;
    trainer = nullptr;
}

void TestingTreeLearning::initialize_hierarchy_with_all_paths(const std::vector<std::vector<size_t>> &subpaths) {
    for (auto& path : subpaths/*.get()*/) {
        std::string x_val = size_vector_to_string(path);
        allNodes.emplace(x_val);
        for (auto& subpath : generateAllPossibleSubpaths(path)) {
            std::string y_val = size_vector_to_string(subpath);
            positiveExamples[x_val].emplace(y_val);
            positiveExamples[y_val].emplace(x_val);
        }
        // Getting the id of the new element that is now added: id
        num_entities_and_classes += tree.addChild(path, id);
        // Adding the id of the associated element, and associating that to the string. That is going ot be used to determine the best common ancestor.
        bimap.put(x_val, id);
    }
}

void TestingTreeLearning::finalizeDataIngestion() {
    // Init category hierarchy
    std::cout << "Initializing the EEEL Engine with some tweaks" << std::endl;
    if (ee_engine_) {
        delete ee_engine_;
        ee_engine_ = nullptr;
    }
    if (trainer) {
        delete trainer;
        trainer = nullptr;
    }
    ee_engine_ = new EEEngine(FULL, dimEmbedding, tree, bimap, num_entities_and_classes, nullptr);
    //EEEngine engine{FULL, dimEmbedding, tree, bimap, num_entities_and_classes};

    std::cout << "Generating the negative examples" << std::endl;
    std::map<std::string, std::set<std::string>> negativeExamples;
    for (const std::string& x : allNodes) {
        const std::set<std::string>& set = positiveExamples[x];
        std::set<std::string>& difference = negativeExamples[x];
        std::set_difference(allNodes.begin(), allNodes.end(), set.begin(), set.end(), std::inserter(difference, difference.begin()));
    }

    std::vector<Datum> batch = ee_engine_->generateData(positiveExamples, negativeExamples);

    // Finished to generate the
    std::cout << "Initializing the EEEL Trainer" << std::endl;
    trainer = new HierarchyLearning (FULL, dimEmbedding, num_entities_and_classes, num_entities_and_classes);

    //batch_learning(batch, trainer, 10);
    //for (int i = 0; i<10; i++)
    trainer->Solve_single(batch);
    //std::cout << trainer.ComputeObjective_single(batch) << std::endl;
}

size_t TestingTreeLearning::getVectorRepresentation(const std::vector<size_t> &current) {
    return bimap.getValue(size_vector_to_string(current));
}

double TestingTreeLearning::similarity(const size_t &lhs, const size_t &rhs) {
    Path path = ee_engine_->entity_category_hierarchy().FindPathBetweenEntities(lhs, rhs);
    path.RefreshAggrDistMetric(trainer->copyCategories());
    return trainer->ComputeDist(lhs, rhs, path);
}

void TestingTreeLearning::generateTopKCandidates(PollMap<double, std::string> &map, const std::vector<size_t> &current) {
    size_t  entity = getVectorRepresentation(current);
    for (size_t e_id = 0; e_id < num_entities_and_classes; ++e_id) {
        double dist = similarity(entity, e_id);
        map.add(dist, bimap.getKey(e_id));
    }
}

void testing_batch_learning_method() {

    size_t maximumBranchingFactor = 4;
    size_t maximumHeight = 3;


    std::vector<std::vector<std::vector<size_t>>> ls = generateCompleteSubgraph(maximumBranchingFactor, maximumHeight);
    TestingTreeLearning testing_dimension_as_branching{maximumBranchingFactor, maximumHeight, 7};
    testing_dimension_as_branching.run(ls);

    TestingTreeLearning testing_dimension_mid_branching_100{maximumBranchingFactor, maximumHeight,
                                                            50};
    testing_dimension_mid_branching_100.run(ls);
}

void batch_learning(const std::vector<Datum> &batch, HierarchyLearning &trainer, int iterations) {
    for (int i = iterations; i > 0; i--) {
        // split vector into sub-vectors each of size n
        size_t n = std::thread::hardware_concurrency() * i;

        // determine number of sub-vectors of size n
        int size = (batch.size() - 1) / n + 1;

        // create array of vectors to store the sub-vectors
        std::vector<Datum> vec[size];

        // each iteration of this loop process next set of n elements
        // and store it in a vector at k'th index in vec
        for (int k = 0; k < size; ++k)
        {
            // get range for next set of n elements
            auto start_itr = std::next(batch.cbegin(), k*n);
            auto end_itr = std::next(batch.cbegin(), k*n + n);

            // allocate memory for the sub-vector
            vec[k].resize(n);

            // code to handle the last sub-vector as it might
            // contain less elements
            if (k*n + n > batch.size()) {
                end_itr = batch.cend();
                vec[k].resize(batch.size() - k*n);
            }

            // copy elements from the input range to the sub-vector
            std::copy(start_itr, end_itr, vec[k].begin());
        }

        for (size_t i = 0; i<size; i++) {
            trainer.Solve_single(vec[i]);
        }

        double  objFinal = 0.0;
        for (size_t i = 0; i<size; i++) {
            objFinal += trainer.ComputeObjective_single(vec[i]);
        }
        std::cout << objFinal << std::endl;

    }
}
