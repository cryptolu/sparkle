/** Time-stamp: <2019-02-03 15:04:54 leo>
 * The general header of the functions dealing with the linear layer.
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


#ifndef LINEAR_LAYER_HPP
#define LINEAR_LAYER_HPP

#include "trail-search.hpp"

std::vector<Word> XOR_combine(Word x, Word y);

// Returns a vector containing all the possible values of L(x) where L
// has the properties in Theorem 3.1.
std::vector<Word> possible_L_output(Word x);


#endif
