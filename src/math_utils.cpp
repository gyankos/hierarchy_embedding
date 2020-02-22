//
// Created by giacomo on 28/12/19.
//

#include "math_utils.h"

std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    std::vector<double> toret;
    const size_t lS = lhs.size(), rS = rhs.size();
    const size_t M = std::max(lS, rS);
    for (size_t i = 0; i<M; i++) {
        if (i<lS && i<rS) {
            toret.emplace_back(lhs[i]+rhs[i]);
        } else if (i<lS) {
            toret.emplace_back(lhs[i]);
        } else {
            toret.emplace_back(rhs[i]);
        }
    }
    return toret;
}

std::vector<double> VectorPow(const std::vector<double> &lhs, double exp) {
    std::vector<double> bar;
    bar.resize(lhs.size());
    std::transform (lhs.begin(), lhs.end(), bar.begin(), [exp](double x) { return std::pow(x, exp); });
    return bar;
}

std::vector<double> operator*(const std::vector<std::vector<double>> &matrix, const std::vector<double> array) {
    std::vector<double> toret;
    toret.resize(matrix.size());
    for (size_t i = 0, m = matrix.size(); i<m; i++) {
        double sum = 0.0;
        for (size_t j = 0, n = std::max(matrix[i].size(), array.size()); j<n; j++)
            sum += matrix[i][j] * array[j];
        toret[i] = sum;
    }
    return toret;
}

double fastsigmoid(const double x) { return (double)0.5 * (Tanh::getInstance().fasttanh(0.5 * x) + 1.); }

double cosine_similarity(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    size_t n = std::min(lhs.size(), rhs.size());
    double score = 0.0, normL = 0.0, normR = 0.0;
    for (size_t i = 0; i<n; i++) {
        score += (lhs[i]*rhs[i]);
        normL += (lhs[i]*lhs[i]);
        normR += (rhs[i]*rhs[i]);
    }
    return score / (std::sqrt(normL)*std::sqrt(normR));

}


std::vector<double> operator+=(std::vector<double> &lhs, const std::vector<double> &rhs) {
    lhs = lhs+rhs;
    return lhs;
}

double max_cosine_similarity(const std::vector<std::vector<double>> &lhs, const std::vector<std::vector<double>> &rhs) {
    double sim = 0;
    for (auto x : lhs) {
        for (auto y: lhs) {
            sim = std::max(sim, cosine_similarity(x,y));
        }
    }
    return sim;
}

Tanh &Tanh::getInstance() {
    static Tanh    instance; // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
}

double Tanh::fasttanh(const double &x) {
    if (std::isnan(x)) return x; // tanh(nan)=nan
    int is_inf=std::isinf(x);
    if (is_inf>0) return 1; // tanh(inf)=1
    if (is_inf<0) return -1; // tanh(-inf)=-1
    if(x>0)
    {
        if(x>MAXTANHX)
            return double(tanhtable[TANHTABLESIZE-1]);
        else
        {
            int i;
            DOUBLE_TO_INT( double(x*((TANHTABLESIZE-1)/MAXTANHX)), i);
            return double(tanhtable[i]);
        }
    }
    else
    {
        double nx = -x;
        if(nx>MAXTANHX)
            return double(-tanhtable[TANHTABLESIZE-1]);
        else
        {
            int i;
            DOUBLE_TO_INT( double(nx*((TANHTABLESIZE-1)/MAXTANHX)), i);
            return double(-tanhtable[i]);
        }
    }
}

Tanh::Tanh() {
    double scaling = MAXTANHX/(TANHTABLESIZE-1);
    for(int i=0; i<TANHTABLESIZE; i++)
        tanhtable.emplace_back(tanh(i*scaling));
}
