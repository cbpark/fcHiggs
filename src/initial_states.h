/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_INITIAL_STATES_H_
#define FCHIGGS_SRC_INITIAL_STATES_H_

#include <cmath>

namespace fchiggs {
class InitPartons {
private:
    double shat_;
    double ymax_;
    double x1_, x2_;

public:
    InitPartons() = delete;
    InitPartons(const double s, const double shat)
        : shat_{shat}, ymax_{-0.5 * std::log(shat_ / s)} {
        init(s);
    }

    double x1() const { return x1_; }
    double x2() const { return x2_; }
    double shat() const { return shat_; }
    double delta_y() const { return 2 * ymax_; }

private:
    void init(const double s);
};

class Rho {
private:
    double m_, gamma_, s_;
    double rho1_, rho2_;

public:
    Rho() = delete;
    Rho(const double qmin, const double qmax, const double m,
        const double gamma, const double s)
        : m_{m}, gamma_{gamma}, s_{s} {
        init(qmin * qmin, qmax * qmax);
    }

    double delta() const { return rho2_ - rho1_; }

    double jacobian(const double val) const;

    /** \hat{s} = M * Gamma * tan(rho) + M^2 */
    double shat(const double val) const;

    friend double rhoValue(const Rho &rho);

private:
    void init(const double qmin2, const double qmax2);
};

double rhoValue(const Rho &rho);
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_INITIAL_STATES_H_
