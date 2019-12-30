//
// Created by giacomo on 29/12/19.
//

#if 0
#include "learning/Dataset.h"

Dataset::Dataset(DistMetricMode metric, const size_t dim_embedding) : metric{metric} , dim_embedding{dim_embedding} {
    num_neg_sample_ = 50;
}

void Dataset::AddDatum(int entity_i, int entity_o, int count) {
    data_.push_back(std::make_pair(entity_i, entity_o));
    count_.push_back(count);
    AddPair(entity_i, entity_o);
}

Datum Dataset::datum(int idx) {
    return Datum(metric, dim_embedding, data_[idx].first, data_[idx].second,count_[idx], num_neg_sample_);
}

const std::set<int> &Dataset::positive_entities(const int entity_i) {
    std::map<int, std::set<int> >::const_iterator it = entity_pairs_.find(entity_i);
    if (it != entity_pairs_.end()) {
        return it->second;
    } else {
        return empty_set_;
    }
}

void Dataset::AddPair(const int entity_i, const int entity_o) {
    entity_pairs_[entity_i].insert(entity_o);
}
#endif