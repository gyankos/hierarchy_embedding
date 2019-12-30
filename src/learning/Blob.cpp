//
// Created by giacomo on 28/12/19.
//

#include "learning/Blob.h"

Blob::Blob(DistMetricMode metric, const size_t num_row, const size_t num_col)
        : num_row_(num_row), num_col_(num_col), m{metric} {
    if (num_col == 1)
        count_ = num_row_;
    else {
        switch (metric) {
            case DIAG:
                count_ = num_row_;
                break;
            case EDIAG:
                count_ = 1;
                break;
            case FULL:
                count_ = num_row_ * num_col_;
                break;
            default:
                throw std::runtime_error("Unknown metric mode");
        }
    }
    data_.resize(count_);
    ClearData();
}

size_t Blob::offset(const size_t r, const size_t c) const {
    if (num_col_ == 1) {
        return r;
    } else {
        switch (m) {
            case DIAG:
                return r;
            case EDIAG:
                return 0;
            case FULL:
                return r * num_col_ + c;
            default:
                throw std::runtime_error("Unknown metric mode");
        }
    }
}


void Blob::init_data_at(const double init_val, const size_t r, const size_t c) {
    data_[offset(r, c)] = init_val;
}

void Blob::ClearData() {
    std::fill(data_.begin(), data_.end(), 0);
}

void Blob::CopyFrom(const Blob &source, const double coeff) {
    for (int idx = 0; idx < count_; ++idx) {
        data_[idx] = source.data_[idx] * coeff;
    }
}

void Blob::Accumulate(const Blob &source, const double coeff, const double d_coeff) {
    for (int idx = 0; idx < count_; ++idx) {
        data_[idx] = source.data_[idx] * coeff + data_[idx] * d_coeff;
    }
}

void Blob::Normalize() {
    float sum = 0;
    for (int idx = 0; idx < count_; ++idx) {
        sum += data_[idx] * data_[idx];
    }
    sum = sqrt(sum);
    for (int idx = 0; idx < count_; ++idx) {
        data_[idx] /= sum;
    }
}

void Blob::Rectify() {
    for (int idx = 0; idx < count_; ++idx) {
        data_[idx] = (data_[idx] > 0 ? data_[idx] : 0);
    }
}

Blob::Blob(const Blob &x) {
    m = x.m;
    data_ = x.data_;
    num_row_ = x.num_row_;
    num_col_ = x.num_col_;
    count_ = x.count_;
}

Blob &Blob::operator=(const Blob &x) {
    m = x.m;
    data_ = x.data_;
    num_row_ = x.num_row_;
    num_col_ = x.num_col_;
    count_ = x.count_;
    return *this;
}
