/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_SIGMA_PPHQ_H_
#define FCHIGGS_SRC_SIGMA_PPHQ_H_

namespace fchiggs {
/**
 * differential cross secion for qin(p1) g(k1) --> qout(p2) H(k2) process.
 */
double dsigma_dcos(const double shat, const double mh, const double mqin,
                   const double mqout, const double alpha_s, const double g,
                   const double gtilde);
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_SIGMA_PPHQ_H_
