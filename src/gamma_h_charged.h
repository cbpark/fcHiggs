/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_GAMMA_H_CHARGED_H_
#define FCHIGGS_SRC_GAMMA_H_CHARGED_H_

#include <ostream>
#include "couplings.h"

namespace fchiggs {
class ChargedHiggsDecayWidth {
private:
    double mh_;
    double gamma_tb_, gamma_cb_, gamma_ub_;
    double gamma_taunu_, gamma_munu_;
    double gamma_wh_;

    double gamma_total_;

public:
    ChargedHiggsDecayWidth() = delete;
    ChargedHiggsDecayWidth(const double mh, const double mh_sm, const Hup &hu,
                           const VHd &v, const Angles &ang)
        : mh_(mh) {
        init_gamma(mh_sm, hu, v, ang);
    }

    double br_tb() const { return gamma_tb_ / gamma_total_; }
    double br_cb() const { return gamma_cb_ / gamma_total_; }
    double br_ub() const { return gamma_ub_ / gamma_total_; }
    double br_taunu() const { return gamma_taunu_ / gamma_total_; }
    double br_munu() const { return gamma_munu_ / gamma_total_; }
    double br_wh() const { return gamma_wh_ / gamma_total_; }

    void printBR() const;

    friend std::ostream &operator<<(std::ostream &os,
                                    const ChargedHiggsDecayWidth &hdec);

private:
    void init_gamma(const double mh_sm, const Hup &hu, const VHd &v,
                    const Angles &ang);
};
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_GAMMA_H_CHARGED_H_
