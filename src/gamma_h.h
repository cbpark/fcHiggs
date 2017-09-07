/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_GAMMA_H_H_
#define FCHIGGS_SRC_GAMMA_H_H_

#include <ostream>
#include "angles.h"
#include "couplings.h"
#include "utils.h"

namespace fchiggs {
using GH3 = ValueType<double>;
using GZPX = ValueType<double>;

class HiggsDecayWidth {
private:
    double gamma_bd_, gamma_bs_;
    double gamma_bb_, gamma_tt_;
    double gamma_tautau_;
    double gamma_ww_, gamma_zz_;
    double gamma_zpzp_;
    double gamma_aa_, gamma_gg_;
    double gamma_hh_;

    double gamma_total_;

public:
    HiggsDecayWidth() = delete;
    HiggsDecayWidth(const double mh, const double mh_sm, const double mzp,
                    const double alpha_s, const GZPX &gx, const GH3 &ghhh,
                    const Hup &cup, const Hdown &cdown, const Angles &ang) {
        init_gamma(mh, mh_sm, mzp, alpha_s, gx.value, ghhh.value, cup, cdown,
                   ang);
    }

    double br_bd() const { return gamma_bd_ / gamma_total_; }
    double br_bs() const { return gamma_bs_ / gamma_total_; }
    double br_bb() const { return gamma_bb_ / gamma_total_; }
    double br_tt() const { return gamma_tt_ / gamma_total_; }
    double br_tautau() const { return gamma_tautau_ / gamma_total_; }
    double br_ww() const { return gamma_ww_ / gamma_total_; }
    double br_zz() const { return gamma_zz_ / gamma_total_; }
    double br_zpzp() const { return gamma_zpzp_ / gamma_total_; }
    double br_aa() const { return gamma_aa_ / gamma_total_; }
    double br_gg() const { return gamma_gg_ / gamma_total_; }
    double br_hh() const { return gamma_hh_ / gamma_total_; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const HiggsDecayWidth &hdec);

private:
    void init_gamma(const double mh, const double mh_sm, const double mzp,
                    const double alpha_s, const double gx, const double ghhh,
                    const Hup &cup, const Hdown &cdown, const Angles &ang);
};
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_GAMMA_H_H_
