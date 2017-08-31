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
double dsigma_dt_dg(const double shat, const double mh, const double alpha_s,
                    const Hdown &hd, const Angles &ang, const DQuark &type) {
    if (shat < std::pow(mh + MB, 2)) { return 0.0; }

    const double coeff = alpha_s / (32.0 * NC * shat * shat);
    double coup = ang.sin_alpha_beta() / ang.cos_beta();
    if (type == DQuark::Down) {
        coup *= hd.c13();
    } else if (type == DQuark::Strange) {
        coup *= hd.c23();
    } else {
        coup *= 0.0;
    }

    CM22 pmom{shat, mh, 0.0};
    const double s = shat, t = pmom.t_hat(), mh2 = mh * mh;
    const double amp2 = -t / s - s / t + MB2 * (1.0 / s + 1.0 / t) -
                        MB2 * mh2 / (t * t) +
                        2.0 * (s - mh2 + MB2) * (mh2 - t) / (s * t);

    return coeff * coup * coup * amp2;
}

double dsigma_dcos_dg(const double shat, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang, const DQuark &type) {
    const double jacobian = 0.5 * lambda12(shat, mh * mh, MB2);
    const double dsigma_dt = dsigma_dt_dg(shat, mh, alpha_s, hd, ang, type);
    return dsigma_dt * jacobian;
}

double dsigma_dt_bg(const double shat, const double mh, const double alpha_s,
                    const Hdown &hd, const Angles &ang) {
    if (shat < std::pow(mh + MB, 2)) { return 0.0; }

    const double coeff = alpha_s / (32.0 * NC * std::pow(shat - MB2, 2));
    const double lambda_b =
        (SQRT2 * MB * ang.cos_alpha() / VEW + hd.c33() * ang.sin_alpha_beta()) /
        ang.cos_beta();

    CM22 pmom{shat, mh, MB};
    const double s = shat, t = pmom.t_hat(), mh2 = mh * mh;
    const double amp2 = -t / s - s / t - 2.0 * MB2 * (1.0 / s + 1.0 / t) -
                        MB2 * (mh2 - MB2) * (1.0 / (s * s) + 1.0 / (t * t)) +
                        2.0 * (s - mh2 + MB2) * (mh2 - MB2 - t) / (s * t) +
                        8.0 * MB2 * (mh2 - 2.0 * MB2) / (s * t);
    return coeff * lambda_b * lambda_b * amp2;
}

double dsigma_dcos_bg(const double shat, const double mh, const double alpha_s,
                      const Hdown &hd, const Angles &ang) {
    const double jacobian =
        0.5 * (1 - MB2 / shat) * lambda12(shat, mh * mh, MB2);
    const double dsigma_dt = dsigma_dt_bg(shat, mh, alpha_s, hd, ang);
    return dsigma_dt * jacobian;
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
