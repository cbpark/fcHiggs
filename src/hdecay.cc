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
#include "gamma_h.h"
#include "pdf.h"
#include "user_interface.h"

using std::to_string;

constexpr char appname[] = "hdecay";

constexpr char PDFNAME[] = "NNPDF23_lo_as_0130_qed";
constexpr double MHSM = 125.0;
constexpr double MZP = 400.0;
constexpr double GZPX = 0.01;
constexpr double MU = 200.0;
constexpr double VS = 500.0;
const double Y33U = SQRT2 * MT / VEW;

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << appname
                  << "<m_H (GeV>) <tan(beta)> <cos(alpha-beta)> [output]\n";
        return 1;
    }

    const double mh = std::atof(argv[1]);
    message(appname,
            "m_H = " + to_string(mh) + " GeV, m_H(SM) = " + to_string(MHSM) +
                " GeV");
    message(appname,
            "mu = " + to_string(MU) + " GeV, v_s = " + to_string(VS) + "GeV");
    const double tan_beta = std::atof(argv[2]);
    const double cos_alpha_beta = std::atof(argv[3]);
    message(appname,
            "tan(beta) = " + to_string(tan_beta) +
                ", cos(alpha-beta) = " + to_string(cos_alpha_beta));
    fchiggs::Angles ang{tan_beta, cos_alpha_beta};
    const fchiggs::Hup cup{ang, Y33U};
    const fchiggs::Hdown cdown{ang};

    auto pdf = fchiggs::mkPdf(PDFNAME);
    const double alpha_s = pdf->alphasQ(mh);
    message(appname, "alpha_s(m_H) = " + to_string(alpha_s));

    fchiggs::HQuartic lambda_h{MHSM, mh, fchiggs::Mu(MU), fchiggs::Vs(VS), ang};
    const double ghhh = lambda_h.trilinear();

    fchiggs::HiggsDecayWidth hdecay{
        mh,  MHSM,  MZP, alpha_s, fchiggs::GZPX(GZPX), fchiggs::GH3(ghhh),
        cup, cdown, ang};

    if (argc == 4) { hdecay.printBR(); }

    if (argc == 5) {
        std::ofstream fout;
        fout.open(argv[4], std::ios_base::app);
        fout << hdecay << '\n';
        message(appname,
                "the output has been saved to `" + std::string(argv[4]) + "'.");
    }
}
