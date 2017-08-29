/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "couplings.h"
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
}  // namespace fchiggs
