/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_pphq.h"
#include <cmath>
#include "kinematics.h"
#include "utils.h"

namespace fchiggs {
double dsigma_dt(const double shat, const double mh, const double mqin,
                 const double mqout, const double alpha_s, const double g,
                 const double gtilde) {
    if (shat < std::pow(mh + mqout, 2)) { return 0.0; }

    CM22 pmom{shat, mh, mqin, mqout};
    double s = shat, t = pmom.t_hat();
    double mh2 = mh * mh, mqin2 = mqin * mqin, mqout2 = mqout * mqout;

    double F1 = s * t - mqin2 * mqout2, F2 = s + t - mqin2 - mqout2;
    double G1 = mh2 - mqout2 - s, G2 = mh2 - mqin2 - t;
    double SS = s - mqin2, TT = t - mqout2;
    double g2 = g * g, gt2 = gtilde * gtilde;

    double sigma =
        (g2 + gt2) * ((2 * F1 - F2 * F2 - 2 * G1 * G2) / (SS * TT) +
                      2 * mqin2 * G1 / (SS * SS) + 2 * mqout2 * G2 / (TT * TT));
    sigma += (g2 - gt2) * 4.0 * mqin * mqout * mh2 / (SS * TT) *
             (1.0 - F1 * F2 / (mh2 * SS * TT));

    sigma *= alpha_s / (8.0 * NC * SS * SS);
    return sigma;
}

double dsigma_dcos(const double shat, const double mh, const double mqin,
                   const double mqout, const double alpha_s, const double g,
                   const double gtilde) {
    double dsigma = dsigma_dt(shat, mh, mqin, mqout, alpha_s, g, gtilde);
    double jacobian = 0.5 * (shat - mqin * mqin) *
                      lambda12(shat, mh * mh, mqout * mqout) / shat;
    return dsigma * jacobian;
}
}  // namespace fchiggs
