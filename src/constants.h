/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_CONSTANTS_H_
#define FCHIGGS_SRC_CONSTANTS_H_

#include <cmath>

constexpr double PI = 3.141592653589793;
constexpr double TWOPI = 2.0 * PI;
constexpr double DELTA = 2;

/** top quark mass */
constexpr double MT = 173.0;
constexpr double MT2 = MT * MT;

/** bottom quark mass */
constexpr double MB = 4.7;
constexpr double MB2 = MB * MB;

/** Z boson mass */
constexpr double MZ = 91.188;

constexpr double GF = 1.1663787e-5;

constexpr double SQRT2 = 1.4142135623730951;

constexpr double VEW2 = 1.0 / (SQRT2 * GF);

const double VEW = std::sqrt(VEW2);

/** alpha_s (MZ) */
constexpr double ALPHAS = 0.118;

constexpr int NC = 3;

/** conversion factor GeV^-2 -> pb */
constexpr double PBCONV = 3.893793656e8;

constexpr double VUD = 0.97434;
constexpr double VUS = 0.22506;
constexpr double VUB = 0.00357;
constexpr double VCD = 0.22492;
constexpr double VCS = 0.97351;
constexpr double VCB = 0.0411;

#endif  // FCHIGGS_SRC_CONSTANTS_H_
