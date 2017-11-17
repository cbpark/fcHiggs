/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "angles.h"
#include <cmath>
#include "constants.h"

namespace fchiggs {
double mpi_half_pi_half(double alpha) {
    while (alpha > PIHALF) { alpha -= PI; }
    while (alpha < -PIHALF) { alpha += PI; }
    return alpha;
}

void Angles::initBetas() {
    double beta = std::atan(tan_beta_);
    cos_beta_ = std::cos(beta);
    sin_beta_ = cos_beta_ * tan_beta_;

    double alpha = beta + acos(cos_alpha_beta_);
    alpha = mpi_half_pi_half(alpha);  // -pi/2 <= alpha <= pi/2
    sin_alpha_ = std::sin(alpha);
    cos_alpha_ = std::cos(alpha);

    sin_alpha_beta_ = std::sin(alpha - beta);
}
}  // namespace fchiggs
