/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "sigma_pph.h"
#include <cmath>
#include <complex>
#include <memory>
#include "LHAPDF/LHAPDF.h"
#include "angles.h"
#include "constants.h"
#include "couplings.h"
#include "initial_states.h"

using std::complex;

namespace fchiggs {
complex<double> ftau(const double tau) {
    if (tau < 1) {
        const double beta = std::sqrt(1 - tau);
        const complex<double> arg{std::log((1 + beta) / (1 - beta)), -PI};
        return -0.25 * std::pow(arg, 2);
    }
    const double arg = std::asin(1.0 / std::sqrt(tau));
    return {arg * arg, 0};
}

complex<double> fTriangle(const double tau) {
    return tau * (complex<double>(1.0, 0) + (1 - tau) * ftau(tau));
}

/*
 * Based on Eq.(3.57) in arXiv:hep-ph/0503172.
 * tau is 1/tau_Q. The 3/4 factor is rescaled to 3/2.
 */
double sigma0(const double mh, const double alpha_s, const Hup &hu,
              const Hdown &hd, const Angles &ang) {
    // SM couplings:
    // const double coup_u = 1.0;
    // const double coup_d = 1.0;
    const double c1 = ang.cos_alpha() / ang.cos_beta();
    const double c2 = VEW * ang.sin_alpha_beta() / (SQRT2 * ang.cos_beta());
    const double coup_u = c1 + c2 * hu.c33() / MT;
    const double coup_d = c1 + c2 * hd.c33() / MB;

    const double mh2 = mh * mh;

    complex<double> a12tau = coup_u * fTriangle(4 * MT2 / mh2);
    a12tau += coup_d * fTriangle(4 * MB2 / mh2);
    const double a12tau_sq = std::norm(a12tau);

    const double coeff = alpha_s * alpha_s * mh2 / (576 * PI * VEW2);

    return coeff * (9.0 / 4.0) * a12tau_sq;
}

double delta(const double shat, const double mh, const double gammah) {
    const double sgammah = shat * gammah / mh;
    return (1.0 / PI) * sgammah /
           (std::pow(shat - mh * mh, 2) + sgammah * sgammah);
}

double sigma_gg(const double mh, const double alpha_s, const Hup &hu,
                const Hdown &hd, const Angles &ang) {
    return sigma0(mh, alpha_s, hu, hd, ang);
}

double sigma_bb(const double mh, const Hdown &hd, const Angles &ang) {
    const double coeff = PI * MB2 / (2.0 * NC * NC * VEW2);
    const double coup =
        ang.cos_alpha() / ang.cos_beta() +
        VEW * ang.sin_alpha_beta() * hd.c33() / (SQRT2 * MB * ang.cos_beta());
    const double beta2 = 1 - 4 * MB2 / (mh * mh);
    return coeff * coup * coup * std::sqrt(beta2);
}

double sigma_db(const Hdown &hd, const Angles &ang, const DQuark &type) {
    const double coeff = PI / (8.0 * NC * NC);
    double coup = ang.sin_alpha_beta() / ang.cos_beta();
    if (type == DQuark::Down) {
        coup *= hd.c13();
    } else if (type == DQuark::Strange) {
        coup *= hd.c23();
    } else {
        coup *= 0.0;
    }
    return coeff * coup * coup;
}

double dsigma_h(std::shared_ptr<LHAPDF::PDF> pdf, const InitPartons &p,
                const double mu, const double mh, const double gammah,
                const double alpha_s, const Hup &hu, const Hdown &hd,
                const Angles &ang) {
    const double x1 = p.x1(), x2 = p.x2();
    const double shat = p.shat();

    // g g --> H
    double sigma = pdf->xfxQ(21, x1, mu) * pdf->xfxQ(21, x2, mu) *
                   sigma_gg(mh, alpha_s, hu, hd, ang);

    // b b --> H
    sigma +=
        pdf->xfxQ(5, x1, mu) * pdf->xfxQ(-5, x2, mu) * sigma_bb(mh, hd, ang);

    const double pdf_b = pdf->xfxQ(5, x2, mu), pdf_bbar = pdf->xfxQ(-5, x2, mu);

    // d b --> H
    auto q_typ = DQuark::Down;
    sigma += (pdf->xfxQ(1, x1, mu) * pdf_bbar + pdf->xfxQ(-1, x1, mu) * pdf_b) *
             sigma_db(hd, ang, q_typ);

    // s b --> H
    q_typ = DQuark::Strange;
    sigma += (pdf->xfxQ(3, x1, mu) * pdf_bbar + pdf->xfxQ(-3, x1, mu) * pdf_b) *
             sigma_db(hd, ang, q_typ);

    return sigma * delta(shat, mh, gammah) / (x1 * x2);
}
}  // namespace fchiggs
