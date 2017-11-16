/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_ppht_charged.h"
#include <memory>
#include "LHAPDF/LHAPDF.h"
#include "angles.h"
#include "couplings.h"
#include "initial_states.h"
#include "sigma_pphq.h"

namespace fchiggs {
double dsigma_dcos_bg_ht(const double shat, const double mh,
                         const double alpha_s, const Hup &hu, const VHd &v,
                         const Angles &ang) {
    double lamL =
        SQRT2 * MB * ang.tan_beta() * VTB / VEW - v.VHd33() / ang.cos_beta();
    double lamR =
        -VTB * (SQRT2 * MT * ang.tan_beta() / VEW - hu.c33() / ang.cos_beta());
    double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
    double sigma = dsigma_dcos(shat, mh, MB, MT, alpha_s, g, gtilde);
    return sigma;
}

double dsigma_dcos_ht(std::shared_ptr<LHAPDF::PDF> pdf, const InitPartons &p,
                      const double mu, const double mh, const double alpha_s,
                      const Hup &hu, const Hdown &hd, const Angles &ang) {
    const double x1 = p.x1(), x2 = p.x2();
    const double shat = p.shat();
    const double pdf_g = pdf->xfxQ(21, x2, mu);

    VHd v{hd};

    double sigma = (pdf->xfxQ(5, x1, mu) + pdf->xfxQ(-5, x1, mu)) * pdf_g *
                   dsigma_dcos_bg_ht(shat, mh, alpha_s, hu, v, ang);
    return sigma / (x1 * x2);
}
}  // namespace fchiggs
