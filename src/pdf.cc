/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#include "pdf.h"
#include <memory>
#include <string>
#include <vector>
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/Info.h"
#include "LHAPDF/LHAPDF.h"
#include "constants.h"

namespace fchiggs {
std::shared_ptr<LHAPDF::PDF> mkPdf(const std::string &pdfname) {
    LHAPDF::Info &cfg{LHAPDF::getConfig()};
    cfg.set_entry("Verbosity", 0);  // make lhapdf quiet
    std::vector<int> flavors = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21};
    cfg.set_entry("Flavors", flavors);

    LHAPDF::AlphaS *alphas{new LHAPDF::AlphaS_ODE()};
    alphas->setQuarkMass(5, MB);
    alphas->setQuarkMass(6, MT);
    alphas->setMZ(MZ);
    alphas->setAlphaSMZ(ALPHAS);
    std::shared_ptr<LHAPDF::PDF> pdf{LHAPDF::mkPDF(pdfname)};
    pdf->setAlphaS(alphas);
    return pdf;
}
}  // namespace fchiggs
