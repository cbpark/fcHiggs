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
constexpr double PI3 = 31.006276680299816;
constexpr double TWOPI = 2.0 * PI;
constexpr double DELTA = 2;

/** top quark mass */
constexpr double MT = 173.0;
constexpr double MT2 = MT * MT;

/** bottom quark mass */
constexpr double MB = 4.7;
constexpr double MB2 = MB * MB;

/** charm quark mass */
constexpr double MC = 1.27;
constexpr double MC2 = MC * MC;

/** muon mass */
constexpr double MMU = 105.658e-3;
constexpr double MMU2 = MMU * MMU;

/** tau lepton mass */
constexpr double MTAU = 1.777;
constexpr double MTAU2 = MTAU * MTAU;

/** W boson mass */
constexpr double MW = 80.419;
constexpr double MW2 = MW * MW;

/** Z boson mass */
constexpr double MZ = 91.188;
constexpr double MZ2 = MZ * MZ;

constexpr double GF = 1.1663787e-5;

constexpr double SQRT2 = 1.4142135623730951;

constexpr double VEW2 = 1.0 / (SQRT2 * GF);

const double VEW = std::sqrt(VEW2);

/** alpha_s (MZ) */
constexpr double ALPHAS = 0.118;

/** alpha (MW) */
constexpr double ALPHA = 1.0 / 128;

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
