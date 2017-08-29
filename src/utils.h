/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_UTILS_H_
#define FCHIGGS_SRC_UTILS_H_

namespace fchiggs {
double getRandom();

inline double costh(const double delta) { return -1.0 + getRandom() * delta; }
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_UTILS_H_