/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_pphq_charged.h"
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

double dsigma_dcos_ug_hb(const double shat, const double mh,
                         const double alpha_s, const VHd &v, const Angles &ang,
                         const UQuark &type) {
    double lamL = 0, lamR = 0;
    if (type == UQuark::Up) {
        lamL = SQRT2 * MB * ang.tan_beta() * VUB / VEW -
               v.VHd13() / ang.cos_beta();
        double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
        return dsigma_dcos(shat, mh, 0, MB, alpha_s, g, gtilde);
    } else if (type == UQuark::Charm) {
        lamL = SQRT2 * MB * ang.tan_beta() * VCB / VEW -
               v.VHd23() / ang.cos_beta();
        lamR = - SQRT2 * MC * ang.tan_beta() * VCB / VEW;
        double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
        return dsigma_dcos(shat, mh, MC, MB, alpha_s, g, gtilde);
    }
}

double dsigma_dcos_hb(std::shared_ptr<LHAPDF::PDF> pdf, const InitPartons &p,
                      const double mu, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang) {
    const double x1 = p.x1(), x2 = p.x2();
    const double shat = p.shat();
    const double pdf_g = pdf->xfxQ(21, x2, mu);

    VHd v{hd};

    // u g --> H b
    auto q_typ = UQuark::Up;
    double sigma = (pdf->xfxQ(2, x1, mu) + pdf->xfxQ(-2, x1, mu)) * pdf_g *
                   dsigma_dcos_ug_hb(shat, mh, alpha_s, v, ang, q_typ);

    // c g --> H b
    q_typ = UQuark::Charm;
    sigma += (pdf->xfxQ(4, x1, mu) + pdf->xfxQ(-4, x1, mu)) * pdf_g *
             dsigma_dcos_ug_hb(shat, mh, alpha_s, v, ang, q_typ);

    return sigma / (x1 * x2);
}
}  // namespace fchiggs
