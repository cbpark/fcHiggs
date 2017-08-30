/*
 *  Copyright (C) 2017 Chan Beom Park <cbpark@gmail.com>
 *
 *  This file is part of fcHiggs, which is released under the GNU General
 *  Public License. See file LICENSE in the top directory of this project
 *  or go to <http://www.gnu.org/licenses/> for full license details.
 */

#ifndef FCHIGGS_SRC_USER_INTERFACE_H_
#define FCHIGGS_SRC_USER_INTERFACE_H_

#include <iostream>
#include <string>

void message(const std::string &appname, const std::string &msg) {
    std::cout << appname << ": " << msg << '\n';
}

#endif  // FCHIGGS_SRC_USER_INTERFACE_H_
