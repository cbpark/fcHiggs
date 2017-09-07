/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "gamma_h.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include "angles.h"
#include "constants.h"
#include "couplings.h"

using std::complex;
using std::setw;

namespace fchiggs {
double gamma_bd(const double mh, const Hdown &c, const Angles &ang,
                const DQuark &typ) {
    if (mh <= MB) { return 0; }

    double coeff = NC * std::pow(ang.sin_alpha_beta(), 2) /
                   (32.0 * PI * std::pow(ang.cos_beta(), 2));
    if (typ == DQuark::Down) {
        coeff *= c.c13();
    } else if (typ == DQuark::Strange) {
        coeff *= c.c23();
    } else {
        coeff *= 0;
    }

    double fac = mh * std::pow(1.0 - MB2 / (mh * mh), 2);
    return coeff * fac;
}

double gamma_qq(const double mh, const double mq, const double lambda_q) {
    if (mh <= 2.0 * mq) { return 0; }

    double coeff = NC * lambda_q * lambda_q / (16.0 * PI);
    double fac = mh * std::pow(1.0 - 4.0 * mq * mq / (mh * mh), 1.5);
    return coeff * fac;
}

double gamma_cc(const double mh, const Angles &ang) {
    double lambda_c = SQRT2 * MC * ang.cos_alpha() / (VEW * ang.cos_beta());
    return gamma_qq(mh, MC, lambda_c);
}

double gamma_bb(const double mh, const Hdown &c, const Angles &ang) {
    double lambda_b = SQRT2 * MB * ang.cos_alpha() / (VEW * ang.cos_beta()) +
                      c.c33() * ang.sin_alpha_beta() / ang.cos_beta();
    return gamma_qq(mh, MB, lambda_b);
}

double gamma_tt(const double mh, const Hup &c, const Angles &ang) {
    double lambda_t = SQRT2 * MT * ang.cos_alpha() / (VEW * ang.cos_beta()) +
                      c.c33() * ang.sin_alpha_beta() / ang.cos_beta();
    return gamma_qq(mh, MT, lambda_t);
}

double gamma_ll(const double mh, const double ml, const Angles &ang) {
    if (mh <= 2.0 * ml) { return 0; }

    double ml2 = ml * ml;
    double coeff = ml2 * std::pow(ang.cos_alpha(), 2) /
                   (8.0 * PI * VEW2 * std::pow(ang.cos_beta(), 2));
    double fac = mh * std::pow(1.0 - 4 * ml2 / (mh * mh), 1.5);
    return coeff * fac;
}

double gamma_mumu(const double mh, const Angles &ang) {
    return gamma_ll(mh, MMU, ang);
}

double gamma_tautau(const double mh, const Angles &ang) {
    return gamma_ll(mh, MTAU, ang);
}

double gamma_vv(const double mh, const double mv, const double coeff) {
    if (mh <= 2.0 * mv) { return 0; }

    double mh2 = mh * mh, mv2 = mv * mv;
    double fac = std::sqrt(1.0 - 4.0 * mv2 / mh2) *
                 (1.0 - 4.0 * mv2 / mh2 + 12.0 * mv2 * mv2 / (mh2 * mh2));
    return coeff * fac;
}

double gamma_ww(const double mh, const Angles &ang) {
    double coeff = std::pow(mh, 3) * std::pow(ang.cos_alpha_beta(), 2) /
                   (16.0 * PI * VEW2);
    return gamma_vv(mh, MW, coeff);
}

double gamma_zz(const double mh, const Angles &ang) {
    double coeff = std::pow(mh, 3) * std::pow(ang.cos_alpha_beta(), 2) /
                   (32.0 * PI * VEW2);
    return gamma_vv(mh, MZ, coeff);
}

double gamma_zpzp(const double mh, const double mzp, const double gx,
                  const Angles &ang) {
    double coeff = std::pow(gx, 4) * std::pow(mh, 3) * VEW2 *
                   std::pow(ang.sin_beta() * ang.sin_alpha(), 2) /
                   (2592.0 * PI * MZ2 * MZ2);
    return gamma_vv(mh, mzp, coeff);
}

complex<double> f_tau(const double tau) {
    if (tau > 1) {
        double beta = std::sqrt(1 - 1.0 / tau);
        complex<double> arg{std::log((1 + beta) / (1 - beta)), -PI};
        return -0.25 * std::pow(arg, 2);
    }
    double arg = std::asin(std::sqrt(tau));
    return {arg * arg, 0};
}

complex<double> loop_12(const double tau) {
    complex<double> num = complex<double>{tau, 0} + (tau - 1) * f_tau(tau);
    return 2.0 / (tau * tau) * num;
}

complex<double> loop_1(const double tau) {
    double tau2 = tau * tau;
    complex<double> num1{2 * tau2 + 3 * tau, 0};
    complex<double> num2 = 3 * (2 * tau - 1) * f_tau(tau);
    return -1.0 / tau2 * (num1 + num2);
}

double gamma_aa(const double mh, const Hup &cup, const Hdown &cdown,
                const Angles &ang) {
    double coeff = ALPHA * ALPHA * std::pow(mh, 3) / (256.0 * PI3 * VEW2);

    double mh2 = mh * mh;
    complex<double> arg_t = NC * (4.0 / 9) * cup.c33() * VEW / (SQRT2 * MT) *
                            loop_12(mh2 / (4.0 * MT2));
    complex<double> arg_b = NC * (1.0 / 9) * cdown.c33() * VEW / (SQRT2 * MB) *
                            loop_12(mh2 / (4.0 * MB2));
    complex<double> arg_tau =
        ang.cos_alpha() / ang.cos_beta() * loop_12(mh2 / (4.0 * MTAU2));
    complex<double> arg_w = ang.cos_alpha_beta() * loop_1(mh2 / (4.0 * MW2));
    double loop_fac = std::norm(arg_t + arg_b + arg_tau + arg_w);

    return coeff * loop_fac;
}

double gamma_gg(const double mh, const Hup &cup, const Hdown &cdown,
                const double alpha_s) {
    double coeff = alpha_s * alpha_s * std::pow(mh, 3) / (72.0 * PI3 * VEW2);

    double mh2 = mh * mh;
    complex<double> arg_t =
        cup.c33() * VEW / (SQRT2 * MT) * loop_12(mh2 / (4 * MT2));
    complex<double> arg_b =
        cdown.c33() * VEW / (SQRT2 * MB) * loop_12(mh2 / (4 * MB2));
    double loop_fac = (9.0 / 16) * std::norm(arg_t + arg_b);

    return coeff * loop_fac;
}

double gamma_hh(const double mh, const double mh_sm, const double ghhh) {
    if (mh <= 2.0 * mh_sm) { return 0; }

    double coeff = ghhh * ghhh * VEW2 / (32.0 * PI * mh);
    double fac = std::sqrt(1.0 - 4 * mh_sm * mh_sm / (mh * mh));
    return coeff * fac;
}

void HiggsDecayWidth::init_gamma(const double mh_sm, const double mzp,
                                 const double alpha_s, const double gx,
                                 const double ghhh, const Hup &cup,
                                 const Hdown &cdown, const Angles &ang) {
    // the factor 2 is to take into account the charge conjugation.
    gamma_bd_ = 2 * gamma_bd(mh_, cdown, ang, DQuark::Down);
    gamma_bs_ = 2 * gamma_bd(mh_, cdown, ang, DQuark::Strange);
    gamma_cc_ = gamma_cc(mh_, ang);
    gamma_bb_ = gamma_bb(mh_, cdown, ang);
    gamma_tt_ = gamma_tt(mh_, cup, ang);
    gamma_mumu_ = gamma_mumu(mh_, ang);
    gamma_tautau_ = gamma_tautau(mh_, ang);
    gamma_ww_ = gamma_ww(mh_, ang);
    gamma_zz_ = gamma_zz(mh_, ang);
    gamma_zpzp_ = gamma_zpzp(mh_, mzp, gx, ang);
    gamma_aa_ = gamma_aa(mh_, cup, cdown, ang);
    gamma_gg_ = gamma_gg(mh_, cup, cdown, alpha_s);
    gamma_hh_ = gamma_hh(mh_, mh_sm, ghhh);

    gamma_total_ = gamma_bd_ + gamma_bs_ + gamma_cc_ + gamma_bb_ + gamma_tt_;
    gamma_total_ += gamma_mumu_ + gamma_tautau_;
    gamma_total_ += gamma_ww_ + gamma_zz_;
    gamma_total_ += gamma_zpzp_;
    gamma_total_ += gamma_aa_ + gamma_gg_;
    gamma_total_ += gamma_hh_;
}

void printOutput(const std::string &mode, const double br) {
    std::cout << "H --> " + mode << ":\t" << br << '\n';
}

void HiggsDecayWidth::printBR() const {
    printOutput("bq (q = d, s)", br_bq());
    printOutput("cc", br_cc());
    printOutput("bb", br_bb());
    printOutput("tt", br_tt());
    printOutput("mu+mu-", br_mumu());
    printOutput("tau+tau-", br_tautau());
    printOutput("ww", br_ww());
    printOutput("zz", br_zz());
    printOutput("z'z'", br_zpzp());
    printOutput("aa", br_aa());
    printOutput("gg", br_gg());
    printOutput("hh", br_hh());
}

std::ostream &operator<<(std::ostream &os, const HiggsDecayWidth &hdec) {
    int width = 12;
    int pre = 8;

    os << std::right << std::fixed << std::setprecision(2) << setw(7)
       << hdec.mh_;
    os << std::setprecision(pre) << setw(width) << hdec.br_bq() << setw(width)
       << hdec.br_cc() << setw(width) << hdec.br_bb() << setw(width)
       << hdec.br_tt();
    os << setw(width) << hdec.br_mumu() << setw(width) << hdec.br_tautau();
    os << setw(width) << hdec.br_ww() << setw(width) << hdec.br_zz();
    os << setw(width) << hdec.br_zpzp();
    os << setw(width) << hdec.br_aa() << setw(width) << hdec.br_gg();
    os << setw(width) << hdec.br_hh();

    return os;
}
}  // namespace fchiggs
