/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "couplings.h"
#include <cmath>
#include "angles.h"
#include "constants.h"

namespace fchiggs {
double Hup::calc_h33u(const Angles &ang, const double y33u) {
    double h33u =
        1.0 -
        VEW2 * ang.cos_beta() * ang.cos_beta() * y33u * y33u / (2.0 * MT2);
    h33u *= SQRT2 * MT / (VEW * ang.sin_beta());
    return h33u;
}

void Hdown::init(const Angles &ang) {
    const double cfac = SQRT2 * MB / (VEW * ang.sin_beta());
    h13d_ = cfac * (VUD * VUB + VCD * VCB);
    h23d_ = cfac * (VUS * VUB + VCS * VCB);
    h33d_ = cfac * (VUB * VUB + VCB * VCB);
}

void HQuartic::init_lambda(const double mh1, const double mh2, const double mu,
                           const double vs) {
    double sin_alpha_sq = std::pow(ang_.sin_alpha(), 2);
    double cos_alpha_sq = std::pow(ang_.cos_alpha(), 2);
    double sin_beta_sq = std::pow(ang_.sin_beta(), 2);
    double cos_beta_sq = std::pow(ang_.cos_beta(), 2);

    double mh1_sq = mh1 * mh1, mh2_sq = mh2 * mh2;
    double mu_vs = SQRT2 * mu * vs;

    lambda1_ = (2 * mh1_sq * sin_alpha_sq + 2 * mh2_sq * cos_alpha_sq -
                mu_vs * ang_.tan_beta()) /
               (4 * VEW2 * cos_beta_sq);
    lambda2_ = (2 * mh1_sq * cos_alpha_sq + 2 * mh2_sq * sin_alpha_sq -
                mu_vs / ang_.tan_beta()) /
               (4 * VEW2 * sin_beta_sq);

    lambda34_ =
        (2 * (-mh1_sq + mh2_sq) * ang_.sin_alpha() * ang_.cos_alpha() + mu_vs) /
        (8 * VEW2 * ang_.sin_beta() * ang_.cos_beta());
}

double HQuartic::trilinear() const {
    double term1 = lambda1_ * ang_.sin_alpha() * ang_.cos_beta() +
                   lambda2_ * ang_.cos_alpha() * ang_.sin_beta();
    term1 *= 6 * ang_.sin_alpha() * ang_.cos_alpha();

    double cos_2alpha =
        std::pow(ang_.cos_alpha(), 2) - std::pow(ang_.sin_alpha(), 2);
    double cos_alpha_beta =
        ang_.cos_alpha() * ang_.cos_beta() - ang_.sin_alpha() * ang_.sin_beta();

    double term2 = 3 * cos_alpha_beta * cos_2alpha - ang_.cos_alpha_beta();
    term2 *= lambda34_;

    return term1 + term2;
}
}  // namespace fchiggs
