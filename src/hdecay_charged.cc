/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "angles.h"
#include "constants.h"
#include "couplings.h"
#include "gamma_h_charged.h"
#include "user_interface.h"

using std::to_string;

constexpr char appname[] = "hdecay_charged";
constexpr double MHSM = 125.0;
const double Y33U = SQRT2 * MT / VEW;

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << appname
                  << "<m_H (GeV>) <tan(beta)> <cos(alpha-beta)> [output]\n";
        return 1;
    }

    const double mh = std::atof(argv[1]);
    message(appname, "m_{H^\\pm} = " + to_string(mh) +
                         " GeV, m_H(SM) = " + to_string(MHSM) + " GeV");
    const double tan_beta = std::atof(argv[2]);
    const double cos_alpha_beta = std::atof(argv[3]);
    message(appname, "tan(beta) = " + to_string(tan_beta) +
                         ", cos(alpha-beta) = " + to_string(cos_alpha_beta));
    fchiggs::Angles ang{tan_beta, cos_alpha_beta};
    const fchiggs::Hup cup{ang, Y33U};
    const fchiggs::Hdown cdown{ang};
    const fchiggs::VHd vhd{cdown};

    fchiggs::ChargedHiggsDecayWidth hdecay{mh, MHSM, cup, vhd, ang};

    if (argc == 4) { hdecay.printBR(); }

    if (argc == 5) {
        std::ofstream fout;
        fout.open(argv[4], std::ios_base::app);
        fout << hdecay << '\n';
        message(appname,
                "the output has been saved to `" + std::string(argv[4]) + "'.");
    }
}
