//
// Created by giacomo on 29/12/19.
//

#ifndef HIERARCHY_TESTS_PATH_H
#define HIERARCHY_TESTS_PATH_H

#include <map>
#include "Blob.h"

class Path {
public:
    Path(DistMetricMode m, size_t dim_embedding) : aggr_dist_metric_(m, dim_embedding, dim_embedding) {};
    Path(const Path& x) = default;

    Path& operator=(const Path& x) = default;

    void RefreshAggrDistMetric(const std::vector<Blob>& categories);

    void AddCategoryNode(const int category_id, const float weight);

    double scale_2() const { return scale_2_; }

    std::vector<size_t>& category_nodes() { return category_nodes_; }
    double category_node_weight(size_t category_id) const;
    const Blob& aggr_dist_metric() const { return aggr_dist_metric_; }


private:
    std::vector<size_t> category_nodes_;
    // category_id => weight in the path
    std::map<size_t, double> category_node_weights_;
    // aggregrated distance metrix
    Blob aggr_dist_metric_;
    double scale_2_;
};



#endif //HIERARCHY_TESTS_PATH_H
