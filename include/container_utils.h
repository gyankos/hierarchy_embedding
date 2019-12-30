//
// Created by giacomo on 30/12/19.
//

#ifndef HIERARCHY_TESTS_CONTAINER_UTILS_H
#define HIERARCHY_TESTS_CONTAINER_UTILS_H

/**
 * Provides the intersection between collection 1 and collection 2 stored in d_first and transformed using t.
 *
 * @tparam InputIt1             Type of the first iterator
 * @tparam InputIt2             Type of the last iterator
 * @tparam OutputIt             Inserter type for the output collection
 * @tparam Transformation       Transformation lambda providing a transformation of the element to be inserted
 * @param first1                Iterator to the start of the first collection
 * @param last1                 Iterator to the end of the first collection
 * @param first2                Iterator to the start of the second collection
 * @param last2                 Iterator to the end of the second collection
 * @param d_first               See {@OutputIt}
 * @param t                     See {@Transformation}
 * @return
 */
template<class InputIt1, class InputIt2, class OutputIt, class Transformation>
OutputIt set_intersection_with_transformation(InputIt1 first1, InputIt1 last1,
                                              InputIt2 first2, InputIt2 last2,
                                              OutputIt d_first, Transformation t)
{
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) {
            ++first1;
        } else  {
            if (!(*first2 < *first1)) {
                *d_first++ = t(*first1++);
            }
            ++first2;
        }
    }
    return d_first;
}

#endif //HIERARCHY_TESTS_CONTAINER_UTILS_H
