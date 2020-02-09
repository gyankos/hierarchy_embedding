//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_JLLEMMA_H
#define HIERARCHY_TESTS_JLLEMMA_H

 // C++ porting of the RandPro Library: https://cran.r-project.org/web/packages/RandPro/RandPro.pdf
#include <cstring>
#include <cmath>

/**
 * Function to determine the required number of dimension for generating the projection matrix
#'
#' Johnson-Lindenstrauss (JL) lemma is the heart of random projection.
#' The lemma states that a small set of points in a high-dimensional space
#' can be embedded into low dimensional space in such a way that distances between the
#' points are nearly preserved.
#' The lemma has been used in dimensionality reduction, compressed sensing, manifold learning and graph embedding.
#' The epsilon is the error tolerant parameter and it is inversely
#' proportional to the accuracy of the result. The higher error tolerant level decreases the number of
#' dimension and also the computation complexity with the marginal loss of accuracy.
 *
 * @param dataSize
 * @param epsilon
 * @return
 */

double dimension_extimate(size_t dataSize, double epsilon = 0.1);

#if 0
enum JLProjections {
    /**
     * The default projection function is "gaussian". In probability theory, Gaussian distribution is also
#'  called as normal distribution. It is a continuous probability distribution used to represent real-valued random
#'  variables. The elements in the random matrix are drawn from N(0,1/k), N is a Natural number and
#'  k value calculated based on JL - Lemma using dimension() function.
     */
    Gaussian,

    /**
     * Achlioptas matrix is easy to generate and also the 2/3rd of the matrix was filled
#' with zero which makes it as more sparse and cut-off the 2/3rd computation.
     */
    Achlioptas,

    /**
     * This method generalizes the achlioptas method and generate very sparse random matrix
#'  to improve the computational speed up of random projection.
     */
    Li,

    /**
     * In this method, the matrix was generated using the equal probability
#' distribution with the elements [-1, 1].
     */
    probability
};

/**
 * Forms the Projection Matrix
 *
 * @param rows
 * @param cols
 * @param JLT
 * @param eps
 * @param projection
 */
void formProjectionMatrix(size_t rows,size_t cols,bool JLT,double eps=0.1,JLProjections projection= Gaussian) {

}
#endif

#endif //HIERARCHY_TESTS_JLLEMMA_H
