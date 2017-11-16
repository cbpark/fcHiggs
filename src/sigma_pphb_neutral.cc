/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_pphb_neutral.h"
#include "LHAPDF/LHAPDF.h"
#include "angles.h"
#include "constants.h"
#include "couplings.h"
#include "initial_states.h"
#include "pdf.h"
#include "sigma_pphq.h"

namespace fchiggs {
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
}  // namespace fchiggs
