//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_EEE_ENGINE_H
#define HIERARCHY_TESTS_EEE_ENGINE_H

// Author: Zhi-Ting Hu, Po-Yao Huang
// Date: 2014.10.26
#pragma once

#include <vector>
#include <stdint.h>
#include <utility>
#include <iostream>
#include <fstream>
#include <learning/Dataset.h>
#include <naryTree.h>
#include "Hierarchy.h"


#define         GET_ENTITY_ID(i)            ((i)*2)
#define         GET_CATEGORY_ID(i)          ((i)*2+1)

class EEEngine {
public:
    EEEngine(DistMetricMode metric, const size_t dim_embedding);

    ~EEEngine();

    void Start();

    // for analysis
    Hierarchy &entity_category_hierarchy() {
        return entity_category_hierarchy_;
    }

    inline int num_entity() { return num_entity_category; }

    //Dataset &train_data() { return train_data_; }

//private:    // private functions

    // TODO: Given a the current batch, all the negative samples are the pairs of the objects that are not within the hierarchy
    //void SampleNegEntities(Datum &datum);

    // TODO: the minibatch will be the set of all the elements that are within the hierarchy
    /*void ThreadCreateMinibatch(const std::vector<int> &next_minibatch_data_idx,
                               std::vector<Datum> &next_minibatch);*/

    /*inline void CopyMinibatch(const std::vector<Datum*>& source,
                              std::vector<Datum*>& target);
    inline void ClearMinibatch(std::vector<Datum*>& minibatch);*/

    void ReadEntityCategoryFile(size_t num_entities_and_categories);

    //void ReadEntityAncestorFile_txt(const std::string &filename);

    // for neg sampling
    void BuildNoiseDistribution();

    //int RandSampleNegEntity();

    Hierarchy entity_category_hierarchy_;

    /**
     * Storing the hierarchy accordingly to the tree information.
     *
     * @param tree              Tree from which we want to store the hierarchies
     * @param levelId           Level id from which the visit is started: initially zero if the tree is the tree node.
     */
    void ReadHierarchyIdWithLevel(naryTree &tree, size_t levelId = 0);

private:

    //Dataset train_data_;
    //Dataset test_data_;
    int num_train_data_;
    //int num_test_data_;


    int32_t num_neg_sample_;
    size_t num_entity_category;

    // for neg sampling
    //std::vector<double> entity_freq_;
    double freq_sum_;
};


#endif //HIERARCHY_TESTS_EEE_ENGINE_H
