/** Time-stamp: <2019-02-03 15:05:42 leo>
 * The body of basic functions.
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

#include "trail-search.hpp"


void print_word(Word x)
{
    std::cout << "[" ;
    for(unsigned int i=0; i<x.size(); i++)
        std::cout << " " << (unsigned int)x[i];
    std::cout << "]" ;
}


Word word_from_int(unsigned int c)
{
    Word result{{0}};
    for(unsigned int i=0; i<B; i++)
        result[i] = (c >> i) & 1;
    return result;
}


void generate_PROBA()
{
    for(unsigned int i=0; i<STARTING_PROBA.size(); i++) // probabilities
        // from Matsui
    {
        PROBA[i+1] = STARTING_PROBA[i] ;
    }
    for(unsigned int i=STARTING_PROBA.size()+1; i<PROBA.size(); i++)
    {
        unsigned int best_proba = 0;
        for(unsigned int j=1; j<i; j++)
        {
            unsigned int new_proba = PROBA[j] + PROBA[i-j];
            if (new_proba > best_proba)
                best_proba = new_proba;
        }
        PROBA[i] = best_proba;
    }
}


Word rotate_word(Word x, int r)
{
    Word result{{0}};
    for(unsigned int i=0; i<B; i++)
        result[i] = x[((int)i+B+r) % B];
    return result;
}


std::vector<unsigned int> unique_up_to_rotation()
{
    unsigned int mask = 0;
    for(unsigned int i=0; i<B; i++)
        mask = (mask << 1) | 1;
    std::vector<unsigned int> result(1, 0);
    for(unsigned int c=1; c<B_SPACE; c++)
    {
        bool c_is_smallest = true;
        for(unsigned int i=1; i<B; i++)
        {
            unsigned int rotated = ((c << i) | (c >> (B-i))) & mask;
            if (rotated < c)
                c_is_smallest = false;
        }
        if (c_is_smallest)
            result.push_back(c);
    }
    return result;
}


unsigned int hamming_weight(Word x)
{
    unsigned int result = 0;
    for (unsigned int i=0; i<B; i++)
        result += (unsigned int)x[i];
    return result;
}
