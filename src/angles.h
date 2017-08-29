/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_ANGLES_H_
#define FCHIGGS_SRC_ANGLES_H_

#include <cmath>

namespace fchiggs {
class Angles {
private:
    double sin_alpha_, cos_alpha_;
    double sin_beta_, cos_beta_;
    double sin_alpha_beta_;

public:
    Angles() = delete;
    Angles(const double sinalpha, const double tanbeta)
        : sin_alpha_(sinalpha),
          cos_alpha_(std::sqrt(1.0 - sin_alpha_ * sin_alpha_)) {
        initBetas(tanbeta);
    }

    double sin_alpha() const { return sin_alpha_; }
    double cos_alpha() const { return cos_alpha_; }
    double sin_beta() const { return sin_beta_; }
    double cos_beta() const { return cos_beta_; }
    double sin_alpha_beta() const { return sin_alpha_beta_; }

private:
    void initBetas(const double tanbeta) {
        const double beta = std::atan(tanbeta);
        cos_beta_ = std::cos(beta);
        sin_beta_ = cos_beta_ * tanbeta;
        sin_alpha_beta_ = sin_alpha_ * cos_beta_ - cos_alpha_ * sin_beta_;
    }
};
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_ANGLES_H_
