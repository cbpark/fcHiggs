/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_UTILS_H_
#define FCHIGGS_SRC_UTILS_H_

#include <cmath>
#include <utility>
#include "constants.h"

namespace fchiggs {
double getRandom();

inline double lambda12(const double x, const double y, const double z) {
    double lambda = x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * z * x;
    return std::sqrt(lambda);
}

inline double costh(const double delta) { return -1.0 + getRandom() * delta; }

inline std::pair<double, double> sigma(const double sum_w,
                                       const double sum_w_sq,
                                       const unsigned int n) {
    const double xsec = sum_w / n;  // cross section
    const double variance = sum_w_sq / n - xsec * xsec;
    const double err = std::sqrt(variance / n);  // error
    return std::make_pair(xsec * PBCONV, err * PBCONV);
}
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_UTILS_H_
