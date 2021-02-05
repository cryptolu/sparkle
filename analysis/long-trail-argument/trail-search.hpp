/** Time-stamp: <2019-04-18 08:36:57 leo>
 * The general header of the trail search program.
 *
 * Copyright (C) 2019 SPARKLEgrupp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#ifndef TRAIL_SEARCH_HPP
#define TRAIL_SEARCH_HPP

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstdlib>

#include <vector>
#include <array>
#include <map>
#include <string>

// Parameters of the search. Change those to specify the search you
// want to perform!
// ================================================================

#define B 3
#define R_LEFT 3
#define R_RIGHT 0
#define ROTATION_AMOUNT 1
#define DIFFERENTIAL 1
#define USE_MATSUI 1
#define PRINT_TRAILS 1
#define RATE_TO_RATE 0



// Constants, types and global variables
// =====================================

#define B_SPACE (1 << B)
#define R_LEFT_SPACE (1 << R_LEFT)
#define R_RIGHT_SPACE (1 << R_RIGHT)
typedef std::array<unsigned int, B> Word; 


// STARTING_PROBA contains -\log_2(p) where p is the probability of a
// trail.
#if DIFFERENTIAL == 1
#define N_KNOWN 12
#else
#define N_KNOWN 8
#endif

extern std::array<unsigned int,N_KNOWN> STARTING_PROBA;

// This array is initialized by the generate_PROBA() function which
// has to be called at the beginning of the main function.
//
// Unlike STARTING_PROBA, we have that PROBA[r] is the opposite of the
// base-2 logarithm of the probability of an r-round trail. In
// particular, it has to hold that PROBA[0] = 0.
typedef std::array<unsigned int, 30> FullProba;

extern FullProba PROBA ;




#include "general_utils.hpp"
#include "linear_layer.hpp"
#include "truncated-trail.hpp"

#endif
