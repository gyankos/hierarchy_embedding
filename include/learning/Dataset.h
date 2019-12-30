//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_DATASET_H
#define HIERARCHY_TESTS_DATASET_H
#if 0

#include <set>
#include <map>
#include <learning/Datum.h>

class Dataset {
    DistMetricMode metric;
    size_t dim_embedding;
public:
    Dataset(DistMetricMode metric, const size_t dim_embedding);;

    ~Dataset() {};

    // Must call AddPair() to finish adding new datum
    void AddDatum(int entity_i, int entity_o, int count);

    //Datum* datum(int idx) { return &(data_[idx]); }
    Datum datum(int idx);

    const std::set<int>& positive_entities(const int entity_i);


    inline void AddPair(const int entity_i, const int entity_o);

private:

private:
    int num_neg_sample_;

    //vector<Datum> data_;
    // < (entity_i, entity_o) >
    std::vector<std::pair<int, int> > data_;
    // < count >
    //std::vector<int> count_;

    // entity_i => (entitiy_o => path)
    //map<int, map<int, Path*> > entity_pair_path_;
    // entity_i => { entity_o }
    std::map<int, std::set<int> > entity_pairs_;
    std::set<int> empty_set_;
};
#endif
#endif //HIERARCHY_TESTS_DATASET_H
