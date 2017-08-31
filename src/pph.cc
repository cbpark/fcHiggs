/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "angles.h"
#include "constants.h"
#include "initial_states.h"
#include "pdf.h"
#include "sigma_pph.h"
#include "user_interface.h"
#include "utils.h"

using std::to_string;

constexpr char appname[] = "pph";

constexpr double ECM = 14000.0;
constexpr double SBEAM = ECM * ECM;
constexpr char PDFNAME[] = "NNPDF23_lo_as_0130_qed";
constexpr unsigned int N = 5000000;
const double Y33U = SQRT2 * MT / VEW;

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage: " << appname
                  << " <m_H (GeV)> <tan(beta)> <cos(alpha-beta)> [output]\n";
        return 1;
    }
    message(appname, "p p --> H");

    message(appname, "E_{CM} = " + to_string(ECM / 1000.0) + " TeV");
    const double mh = std::atof(argv[1]);
    message(appname, "m_H = " + to_string(mh) + " GeV");
    const double gammaH = mh / 10000.0;
    const double qmin = mh / 2.0, qmax = std::sqrt(SBEAM), mtr = mh,
                 gtr = mh / 2.0;
    const double mu = mh;
    const fchiggs::Rho rho{qmin, qmax, mtr, gtr, SBEAM};
    double val = fchiggs::rhoValue(rho);

    auto pdf = fchiggs::mkPdf(PDFNAME);
    const double alpha_s = pdf->alphasQ(mu);

    const double tan_beta = std::atof(argv[2]);
    const double cos_alpha_beta = std::atof(argv[3]);
    message(appname,
            "tan(beta) = " + to_string(tan_beta) +
                ", cos(alpha-beta) = " + to_string(cos_alpha_beta));

    const fchiggs::Angles ang{tan_beta, cos_alpha_beta};
    const fchiggs::Hup hu{ang, Y33U};
    const fchiggs::Hdown hd{ang};

    message(appname, "integrating for cross section ...");
    double sum_w = 0, sum_w_sq = 0;
    for (auto itry = 0; itry != N; ++itry, val = fchiggs::rhoValue(rho)) {
        const double shat = rho.shat(val);
        const fchiggs::InitPartons p{SBEAM, shat};
        const double w =
            fchiggs::dsigma_h(pdf, p, mu, mh, gammaH, alpha_s, hu, hd, ang) *
            rho.delta() * p.delta_y() * rho.jacobian(val);
        sum_w += w;
        sum_w_sq += w * w;
    }

    auto result = fchiggs::sigma(sum_w, sum_w_sq, N);
    const double sigma = result.first, err = result.second;
    message(appname, "... done.");
    message(appname,
            "total cross section = " + to_string(sigma) + " +- " +
                to_string(err) + " pb");

    if (argc == 5) {
        std::ofstream fout;
        fout.open(argv[4], std::ios_base::app);
        fout << std::right << std::fixed << std::setw(7) << std::setprecision(2)
             << mh << std::setw(14) << std::setprecision(9) << sigma
             << std::setw(14) << err << '\n';
        message(appname,
                "the output has been saved to `" + std::string(argv[4]) + "'.");
    }
}
