/** Time-stamp: <2019-02-03 14:58:33 leo>
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


#include "trail-search.hpp"



std::vector<Word> XOR_combine(Word x, Word y)
{
    std::vector<Word> result(1, Word{{0}});
    for(unsigned int i=0; i<B; i++)
    {
        if (x[i] == 1 and y[i] == 1)
        {
            std::vector<Word> added(result);
            for(unsigned int j=0; j<result.size(); j++)
            {
                // case where the sum cancels out
                result[j][i] = 0;
                // case where it doesn't
                added[j][i] = 1;
            }
            result.insert(result.end(), added.begin(), added.end());
        }
        else
        {
            for(unsigned int j=0; j<result.size(); j++)
                result[j][i] = (x[i] | y[i]) ;
        }
    }
    return result;
}


#if B>2
std::vector<Word> possible_L_output(Word x)
{
    unsigned int w = hamming_weight(x);
    std::vector<Word> result;
    if (w == 0)
    {
        result.push_back(Word{{0}});
    }
    else if (w == 1)
    {
        Word all_1;
        for (unsigned int i=0; i<B; i++)
            all_1[i] = 1;
        result.push_back(all_1);
    }
    else if (w == 2)
    {
        Word 
            identical,              // a copy of x
            all_1,                  // all outputs are set
            all_but_first,          // all outputs are set but the first
            // non-zero in x
            all_but_second          // all outputs are set but the second
            // non-zero in x
            ;
        bool first_seen = false;
        // rotation of input
        for (unsigned int i=0; i<x.size(); i++)
        {
            all_1[i] = 1;
            identical[i] = x[i];
            if (x[i] == 1)
            {
                if (not first_seen)
                {
                    all_but_first[i] = 0;
                    all_but_second[i] = 1;
                    first_seen = true;
                }
                else
                {
                    all_but_first[i] = 1;
                    all_but_second[i] = 0;
                }
            }
            else
            {
                all_but_first[i] = 1;
                all_but_second[i] = 1;
            }            
        }
        result.push_back(identical);
        result.push_back(all_1);
        result.push_back(all_but_first);
        result.push_back(all_but_second);
    }
    else if (w == B)
        // all outputs are possible in this case
    {
        for (unsigned int c=1; c<B_SPACE; c++)
        {
            Word l{{0}};
            for (unsigned int i=0; i<B; i++)
                l[i] = (c >> i) & 1;
            result.push_back(l);            
        }
    }
    else
        // the branching number and the impossibility of a weight 1 output
        // constrain this case
    {
        for (unsigned int c=1; c<B_SPACE; c++)
        {
            Word l{{0}};
            for (unsigned int i=0; i<B; i++)
                l[i] = (c >> i) & 1;
            unsigned int h = hamming_weight(l);
            if ((h > 1) and (h + w >= 4))
                result.push_back(l);            
        }
    }
    return result;
}

#else
std::vector<Word> possible_L_output(Word x)
{
    unsigned int w = hamming_weight(x);
    std::vector<Word> result;
    if (w == 0)
    {
        result.push_back(Word{{0}});
    }
    else
    {
        for (unsigned int c=1; c<B_SPACE; c++)
        {
            Word l{{0}};
            for (unsigned int i=0; i<B; i++)
                l[i] = (c >> i) & 1;
            unsigned int h = hamming_weight(l);
            if (h + w >= 3)
                result.push_back(l);            
        }
    }
    return result;
}

#endif

