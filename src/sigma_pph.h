/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_SIGMA_PPH_H_
#define FCHIGGS_SRC_SIGMA_PPH_H_

#include <memory>
#include "LHAPDF/LHAPDF.h"
#include "angles.h"
#include "couplings.h"
#include "initial_states.h"

namespace fchiggs {
double dsigma_h(std::shared_ptr<LHAPDF::PDF> pdf, const InitPartons &p,
                const double mu, const double mh, const double gammah,
                const double alpha_s, const Hup &hu, const Hdown &hd,
                const Angles &ang, const double kgg);
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_SIGMA_PPH_H_
