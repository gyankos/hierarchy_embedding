//
// Created by giacomo on 28/12/19.
//

#ifndef HIERARCHY_TESTS_BLOB_H
#define HIERARCHY_TESTS_BLOB_H

#include <stdexcept>
#include <vector>
#include <cmath>

enum DistMetricMode {
    DIAG, // diagonal matrix
    EDIAG, // diagonal matirx with tied diagonal elements
    FULL // full PSD matrix
};

class Blob {
    DistMetricMode m;
public:
    Blob(DistMetricMode metric, const size_t num_row, const size_t num_col = 1);
    Blob(const Blob& x);
    Blob& operator=(const Blob& x);

    //~Blob() { delete data_; };
    // read only
    const std::vector<double> data() const { return data_; }
    // mutable
    //float* mutable_data() { return data_; }

    inline size_t count() const { return count_; }
    inline size_t num_row() const { return num_row_; }
    inline size_t num_col() const { return num_col_; }

    inline double data_at(const size_t r, const size_t c = 0) const {
        return (data_[offset(r, c)]);
    }
    inline double get(size_t i) const {
        return data_[i];
    }

    void incr_data_at(size_t i, const double t) {
        data_[i] += t;
    }
    void init_data_at(const double init_val, const size_t r, const size_t c = 0);

    void ClearData();

    void CopyFrom(const Blob& source, const double coeff = 1.0);

    // data_ = source * coeff + data_ * d_coeff
    void Accumulate(const Blob& source, const double coeff = 1.0,
                    const double d_coeff = 1.0);

    void Normalize();

    void Rectify();

private:
    size_t offset(const size_t r, const size_t c) const;


private:
    std::vector<double> data_;

    // first dimention
    size_t num_row_;
    // second dimention
    // = 1 if data_ is a vector
    size_t num_col_;
    // number of elements
    size_t count_;
};


#endif //HIERARCHY_TESTS_BLOB_H
