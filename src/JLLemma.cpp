//
// Created by giacomo on 30/12/19.
//

#include "JLLemma.h"

double dimension_extimate(size_t dataSize, double epsilon) {
    double denominator = ( std::pow(epsilon,2) / 2 - std::pow(epsilon, 3) / 3 );
    return ( floor((4 * std::log(dataSize) / denominator)));
}
