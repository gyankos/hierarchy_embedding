//
// Created by giacomo on 28/12/19.
//

#ifndef HIERARCHY_TESTS_DATUM_H
#define HIERARCHY_TESTS_DATUM_H

#include <map>
#include "Path.h"

class Datum {
public:
    // must call AddPath() to finish initialization
    Datum(DistMetricMode metric, const size_t dim_embedding, const size_t entity_i, const size_t entity_o,
          const size_t num_neg_sample = 0, const size_t count = 1);


    Datum() : Datum{FULL, 100, 0, 0} {}

    void AddPath(Path& path);

    const size_t entity_i() { return entity_i_; }
    const size_t entity_o() { return entity_o_; }
    const size_t count() { return count_; }

    Path& category_path();
    /// used in optimization
    void AddNegSample(size_t neg_idx, size_t neg_entity_id, Path& path);

    void ClearEntityGrads();
    void ClearCategoryGrads();

    Blob& entity_i_grad() { return entity_i_grad_; }
    Blob& entity_o_grad() { return entity_o_grad_; }
    const int neg_entity(const int neg_idx) { return neg_entity_id_[neg_idx]; }
    std::vector<Path>& neg_category_paths() { return neg_category_paths_; }
    Path& neg_category_path(const int neg_idx) { return neg_category_paths_[neg_idx]; }
    Blob& neg_entity_grad(const int neg_idx) { return neg_entity_grads_[neg_idx]; }
    Blob& category_grad(const int category_id) { return category_grads_[category_index_[category_id]]; }
    const std::vector<Blob>& category_grads() { return category_grads_; }
    const std::map<size_t, size_t>& category_index() { return category_index_; }

private:

    DistMetricMode m;
    size_t dim_embedding;
    /// Observed data

    size_t entity_i_;
    size_t entity_o_;
    size_t count_;
    Path category_path_;

    /// Used in optimization

    Blob entity_i_grad_;
    Blob entity_o_grad_;

    // negative samples
    std::vector<int> neg_entity_id_;
    // paths between entity_i_ and neg_entities
    std::vector<Path> neg_category_paths_;
    std::vector<Blob> neg_entity_grads_;

    // category_id => index in category_grads_
    // Note: includes categories in BOTH category_path_ and neg_category_paths_
    //   Clear it whenever negative sampling, and re-insert categories in
    //   category_path_
    // TODO: or we can define seperate categroy_grads for category_path_
    std::map<size_t, size_t> category_index_;
    std::vector<Blob> category_grads_;
};


#endif //HIERARCHY_TESTS_DATUM_H
