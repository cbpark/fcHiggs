/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_PDF_H_
#define FCHIGGS_SRC_PDF_H_

#include <memory>
#include <string>
#include "LHAPDF/LHAPDF.h"

namespace fchiggs {
std::shared_ptr<LHAPDF::PDF> mkPdf(const std::string &pdfname);
}  // namespace fchiggs

#endif  // FCHIGGS_SRC_PDF_H_
