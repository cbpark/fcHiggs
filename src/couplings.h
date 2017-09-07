/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_COUPLINGS_H_
#define FCHIGGS_SRC_COUPLINGS_H_

#include "angles.h"
#include "utils.h"

namespace fchiggs {
enum class DQuark { Down, Strange, Bottom };

class Hup {
private:
    double h31u_, h32u_, h33u_;

public:
    Hup() = delete;
    Hup(const Angles &ang, const double y33u)
        : h31u_{0.0}, h32u_{0.0}, h33u_{calc_h33u(ang, y33u)} {}

    double c31() const { return h31u_; }
    double c32() const { return h32u_; }
    double c33() const { return h33u_; }

private:
    double calc_h33u(const Angles &ang, const double y33u);
};

class Hdown {
private:
    double h13d_, h23d_, h33d_;

public:
    Hdown() = delete;
    explicit Hdown(const Angles &ang) { init(ang); }

    double c13() const { return h13d_; }
    double c23() const { return h23d_; }
    double c33() const { return h33d_; }

private:
    void init(const Angles &ang);
};

using Mu = ValueType<double>;
using Vs = ValueType<double>;

class HQuartic {
private:
    Angles ang_;
    double lambda1_, lambda2_, lambda34_;

public:
    HQuartic() = delete;
    HQuartic(const double mh1, const double mh2, const Mu &mu, const Vs &vs,
             const Angles &ang)
        : ang_(ang) {
        init_lambda(mh1, mh2, mu.value, vs.value);
    }

    double lambda1() const { return lambda1_; }
    double lambda2() const { return lambda2_; }
    double lambda34() const { return lambda34_; }

    double trilinear() const;

private:
    void init_lambda(const double mh1, const double mh2, const double mu,
                     const double vs);
};
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_COUPLINGS_H_
