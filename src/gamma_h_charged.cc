/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "gamma_h_charged.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include "angles.h"
#include "constants.h"
#include "couplings.h"
#include "utils.h"

using std::setw;

namespace fchiggs {
double gamma_qq(const double mh, const double mq1, const double mq2,
                const double g, const double gtilde) {
    double coeff = NC * mh / (8 * PI);

    double g2 = g * g, gtilde2 = gtilde * gtilde;
    double mh_sq = mh * mh, mq1_sq = mq1 * mq1, mq2_sq = mq2 * mq2;
    double coup = (g2 + gtilde2) * (1.0 - (mq1_sq + mq2_sq) / mh_sq) -
                  2.0 * (g2 - gtilde2) * mq1 * mq2 / mh_sq;
    double beta = lambda12(1.0, mq1_sq / mh_sq, mq2_sq / mh_sq);

    return coeff * coup * beta;
}

double gamma_tb(const double mh, const Hup &hu, const VHd &v,
                const Angles &ang) {
    double lamL =
        SQRT2 * MB * ang.tan_beta() * VTB / VEW - v.VHd33() / ang.cos_beta();
    double lamR =
        -VTB * (SQRT2 * MT * ang.tan_beta() / VEW - hu.c33() / ang.cos_beta());
    double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
    double gamma = gamma_qq(mh, MT, MB, g, gtilde);
    return gamma;
}

double gamma_cb(const double mh, const VHd &v, const Angles &ang) {
    double lamL =
        SQRT2 * MB * ang.tan_beta() * VCB / VEW - v.VHd23() / ang.cos_beta();
    double lamR = -VCB * SQRT2 * MC * ang.tan_beta() / VEW;
    double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
    double gamma = gamma_qq(mh, MC, MB, g, gtilde);
    return gamma;
}

double gamma_ub(const double mh, const VHd &v, const Angles &ang) {
    double lamL =
        SQRT2 * MB * ang.tan_beta() * VUB / VEW - v.VHd13() / ang.cos_beta();
    double lamR = 0.0;
    double g = (lamL + lamR) / 2.0, gtilde = (lamL - lamR) / 2.0;
    double gamma = gamma_qq(mh, 0.0, MB, g, gtilde);
    return gamma;
}

double gamma_lnu(const double mh, const double ml, const Angles &ang) {
    double ml2 = ml * ml;
    double coeff = ml2 * std::pow(ang.tan_beta(), 2) * mh / (8 * PI * VEW2);
    double beta = std::pow(1.0 - ml2 / (mh * mh), 2);
    double gamma = coeff * beta;
    return gamma;
}

double gamma_taunu(const double mh, const Angles &ang) {
    return gamma_lnu(mh, MTAU, ang);
}

double gamma_munu(const double mh, const Angles &ang) {
    return gamma_lnu(mh, MMU, ang);
}

double gamma_wh(const double mh, const double mh_sm, const Angles &ang) {
    double coeff =
        G2 * std::pow(ang.cos_alpha_beta(), 2) * mh * mh * mh / (64 * PI * MW2);

    double mh2 = mh * mh, mh_sm2 = mh_sm * mh_sm;
    double beta = std::pow(lambda12(1.0, MW2 / mh2, mh_sm2 / mh2), 3);

    double gamma = coeff * beta;
    return gamma;
}

void ChargedHiggsDecayWidth::init_gamma(const double mh_sm, const Hup &hu,
                                        const VHd &v, const Angles &ang) {
    gamma_tb_ = gamma_tb(mh_, hu, v, ang);
    gamma_cb_ = gamma_cb(mh_, v, ang);
    gamma_ub_ = gamma_ub(mh_, v, ang);

    gamma_taunu_ = gamma_taunu(mh_, ang);
    gamma_munu_ = gamma_munu(mh_, ang);

    gamma_wh_ = gamma_wh(mh_, mh_sm, ang);

    gamma_total_ = gamma_tb_ + gamma_cb_ + gamma_ub_;
    gamma_total_ += gamma_taunu_ + gamma_munu_;
    gamma_total_ += gamma_wh_;
}

void printOutput(const std::string &mode, const double br) {
    std::cout << "H^\\pm --> " + mode << ":\t" << br << '\n';
}

void ChargedHiggsDecayWidth::printBR() const {
    printOutput("t b", br_tb());        // (2)
    printOutput("c b", br_cb());        // (3)
    printOutput("u b", br_ub());        // (4)
    printOutput("tau nu", br_taunu());  // (5)
    printOutput("mu nu", br_munu());    // (6)
    printOutput("W h", br_wh());        // (7)
}

std::ostream &operator<<(std::ostream &os, const ChargedHiggsDecayWidth &hdec) {
    int width = 12;
    int pre = 8;

    os << std::right << std::fixed << std::setprecision(2) << setw(7)
       << hdec.mh_;
    os << std::setprecision(pre) << setw(width) << hdec.br_tb() << setw(width)
       << hdec.br_cb() << setw(width) << hdec.br_ub();
    os << setw(width) << hdec.br_taunu() << setw(width) << hdec.br_munu();
    os << setw(width) << hdec.br_wh();

    return os;
}
}  // namespace fchiggs
