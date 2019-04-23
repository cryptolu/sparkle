/** Time-stamp: <2019-02-03 15:07:55 leo>
 * The header declaring basic functions.
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

#ifndef GENERAL_UTILS_HPP
#define GENERAL_UTILS_HPP

#include "trail-search.hpp"

void print_word(Word x);

Word word_from_int(unsigned int c);

void generate_PROBA();

Word rotate_word(Word x, int r);

// Returns a list of all words that are unique up to rotation. For
// example, it cannot contain both [0,0,1] and [0,1,0] as their a
// rotation of one another.
std::vector<unsigned int> unique_up_to_rotation();

// Returns the number of entries in x that set to 1.
unsigned int hamming_weight(Word x);


#endif
