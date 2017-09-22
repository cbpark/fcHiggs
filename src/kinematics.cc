/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "kinematics.h"
#include <cmath>
#include <ostream>
#include "constants.h"
#include "utils.h"

namespace fchiggs {
std::ostream &operator<<(std::ostream &os, const FourMomentum &p) {
    os << "e = " << p.e_ << ", px = " << p.px_ << ", py = " << p.py_
       << ", pz = " << p.pz_;
    return os;
}

FourMomentum boostZ(const FourMomentum &p, const double beta) {
    const double gamma = 1.0 / std::sqrt(1 - beta * beta);
    const double gb = gamma * beta;
    return {gamma * p.e_ - gb * p.pz_, p.px_, p.py_,
            -gb * p.e_ + gamma * p.pz_};
}

double phi_mpi_pi(double phi) {
    while (phi >= PI) { phi -= TWOPI; }
    while (phi < -PI) { phi += TWOPI; }
    return phi;
}

double deltaPhi(const FourMomentum &p1, const FourMomentum &p2) {
    return phi_mpi_pi(p1.phi() - p2.phi());
}

double deltaR(const FourMomentum &p1, const FourMomentum &p2) {
    const double deta = p1.eta() - p2.eta();
    const double dphi = deltaPhi(p1, p2);
    return std::sqrt(deta * deta + dphi * dphi);
}

void CM22::init() {
    const double e = 2.0 * std::sqrt(s_);
    pin_ = (s_ - mqin2_) / e;
    pfin_ = lambda12(s_, mh2_, mqout2_) / e;
}

FourMomentum CM22::p1() const {
    const double energy = (s_ + mqin2_) / (2.0 * std::sqrt(s_));
    return {energy, 0, 0, pin_};
}

FourMomentum CM22::k1() const {
    const double energy = (s_ - mqin2_) / (2.0 * std::sqrt(s_));
    return {energy, 0, 0, -pin_};
}

FourMomentum CM22::p2() const {
    const double energy = (s_ - mh2_ + mqout2_) / (2.0 * std::sqrt(s_));
    return {energy, -pfin_ * sinth_ * std::cos(phi_),
            -pfin_ * sinth_ * std::sin(phi_), -pfin_ * costh_};
}

FourMomentum CM22::k2() const {
    const double energy = (s_ + mh2_ - mqout2_) / (2.0 * std::sqrt(s_));
    return {energy, pfin_ * sinth_ * std::cos(phi_),
            pfin_ * sinth_ * std::sin(phi_), pfin_ * costh_};
}
}  // namespace fchiggs
