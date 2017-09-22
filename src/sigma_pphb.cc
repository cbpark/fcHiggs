/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_pphb.h"
#include "LHAPDF/LHAPDF.h"
#include "angles.h"
#include "constants.h"
#include "couplings.h"
#include "initial_states.h"
#include "kinematics.h"
#include "pdf.h"
#include "utils.h"

namespace fchiggs {
/**
 * differential cross secion for qin(p1) g(k1) --> qout(p2) H(k2) process.
 */
double dsigma_dt(const double shat, const double mh, const double mqin,
                 const double mqout, const double alpha_s, const double g,
                 const double gtilde) {
    if (shat < std::pow(mh + mqout, 2)) { return 0.0; }

    CM22 pmom{shat, mh, mqin, mqout};
    const double s = shat, t = pmom.t_hat();
    const double mh2 = mh * mh, mqin2 = mqin * mqin, mqout2 = mqout * mqout;

    const double F1 = s * t - mqin2 * mqout2, F2 = s + t - mqin2 - mqout2;
    const double G1 = mh2 - mqout2 - s, G2 = mh2 - mqin2 - t;
    const double SS = s - mqin2, TT = t - mqout2;
    const double g2 = g * g, gt2 = gtilde * gtilde;

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
    const double dsigma = dsigma_dt(shat, mh, mqin, mqout, alpha_s, g, gtilde);
    const double jacobian = 0.5 * (shat - mqin * mqin) *
                            lambda12(shat, mh * mh, mqout * mqout) / shat;
    return dsigma * jacobian;
}

double dsigma_dcos_dg(const double shat, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang, const DQuark &type) {
    double g = ang.sin_alpha_beta() / (2 * SQRT2 * ang.cos_beta());
    if (type == DQuark::Down) {
        g *= hd.c13();
    } else if (type == DQuark::Strange) {
        g *= hd.c23();
    } else {
        g *= 0;
    }
    double gtilde = g;
    return dsigma_dcos(shat, mh, 0.0, MB, alpha_s, g, gtilde);
}

double dsigma_dcos_bg(const double shat, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang) {
    double lambda_b = SQRT2 * MB * ang.cos_alpha() / (VEW * ang.cos_beta()) +
                      hd.c33() * ang.sin_alpha_beta() / ang.cos_beta();
    double g = lambda_b / (2 * SQRT2);
    double gtilde = g;
    return dsigma_dcos(shat, mh, MB, MB, alpha_s, g, gtilde);
}

double dsigma_dcos_hb(std::shared_ptr<LHAPDF::PDF> pdf, const InitPartons &p,
                      const double mu, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang) {
    const double x1 = p.x1(), x2 = p.x2();
    const double shat = p.shat();
    const double pdf_g = pdf->xfxQ(21, x2, mu);

    // d g --> H b
    auto q_typ = DQuark::Down;
    double sigma = (pdf->xfxQ(1, x1, mu) + pdf->xfxQ(-1, x1, mu)) * pdf_g *
                   dsigma_dcos_dg(shat, mh, alpha_s, hd, ang, q_typ);

    // s g --> H b
    q_typ = DQuark::Strange;
    sigma += (pdf->xfxQ(3, x1, mu) + pdf->xfxQ(-3, x1, mu)) * pdf_g *
             dsigma_dcos_dg(shat, mh, alpha_s, hd, ang, q_typ);

    // b g --> H b
    sigma += (pdf->xfxQ(5, x1, mu) + pdf->xfxQ(-5, x1, mu)) * pdf_g *
             dsigma_dcos_bg(shat, mh, alpha_s, hd, ang);

    return sigma / (x1 * x2);
}
}  // namespace fchiggs2
