//
// Created by giacomo on 28/12/19.
//

#ifndef HIERARCHY_TESTS_MATH_UTILS_H
#define HIERARCHY_TESTS_MATH_UTILS_H

#include <vector>
#include <algorithm>
#include <cmath>

std::vector<double> operator+(const std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> operator+=(std::vector<double>& lhs, const std::vector<double>& rhs);
std::vector<double> VectorPow(const std::vector<double>& lhs, double exp);

// MatrixVectorMultiply
std::vector<double> operator*(const std::vector<std::vector<double>>& matrix, const std::vector<double> array);

#define TANHTABLESIZE 5000
#define MAXTANHX 10.
#if defined(LINUX) && !defined(__INTEL_COMPILER)  // note: intel compiler on SGI does not like that
#define DOUBLE_TO_INT(in,out) __asm__ __volatile__ ("fistpl %0" : "=m" (out) : "t" (in) : "st")
#else
#define DOUBLE_TO_INT(in, out)  out = int(round(in))
#endif

class Tanh
{
public:
    static Tanh& getInstance();

    double fasttanh(const double&x);

private:
    std::vector<double> tanhtable;
    Tanh();                    // Constructor? (the {} brackets) are needed here.

    // C++ 11
    // =======
    // We can use the better technique of deleting the methods
    // we don't want.
public:
    Tanh(Tanh const&)               = delete;
    void operator=(Tanh const&)  = delete;


};


inline double fastsigmoid(const double x);

double cosine_similarity(const std::vector<double>& lhs, const std::vector<double>& rhs);
#endif //HIERARCHY_TESTS_MATH_UTILS_H
