/** Time-stamp: <2019-03-13 16:53:59 leo>
 * The general header of the class handling truncated trails.
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


// returns the number of steps in the trail considered
size_t TruncatedTrail::size()
{
    return b_left.size()-1;
}

  
// returns the score of the current trail 
unsigned int TruncatedTrail::score(
    unsigned int r_p_s,                  // rounds per step
    FullProba pr // probabilities
    )
{
    unsigned int result = 0;
    // adding finished long trails
    for(auto &t : ltd)
        result += pr[t.first * r_p_s] * t.second;
    // adding ongoing long trails
    for(unsigned int i=0; i<B; i++)
    {
        if (counters_left[i] != 0)
            result += pr[counters_left[i]  * r_p_s];
        if (counters_right[i] != 0)
            result += pr[counters_right[i] * r_p_s];
    }
    return result;
}

  
// returns all possible results of absorbing a_left and a_right
std::vector<TruncatedTrail> TruncatedTrail::absorb(Word a_left, Word a_right)
{
    // updating the counters where the lt:s are cut
    for(unsigned int i=0; i<B; i++)
    {
        if (a_left[i] == 1 and b_left.back()[i] == 1)
        {
            ltd[counters_left[i]] += 1; 
            counters_left[i] = 0;
        }
        if (a_right[i] == 1 and b_right.back()[i] == 1)
        {
            ltd[counters_right[i]] += 1; 
            counters_right[i] = 0;
        }          
    }
    // generating all possible successors
    std::vector<TruncatedTrail> result;
    Word 
        original_b_left(b_left.back()),
        original_b_right(b_right.back());
    std::vector<Word>
        new_b_feistel(b_feistel),
        new_absorbed_left(absorbed_left),
        new_absorbed_right(absorbed_right);
    std::vector<unsigned int> new_absorption_indices(absorption_indices);
    new_absorption_indices.push_back(size());
    new_absorbed_left.push_back(a_left);
    new_absorbed_right.push_back(a_right);
    for(auto & new_left : XOR_combine(a_left, original_b_left))
        for(auto & new_right : XOR_combine(a_right, original_b_right))
        {
            std::vector<Word>
                new_b_left(b_left),
                new_b_right(b_right);
            new_b_left.back() = new_left;
            new_b_right.back() = new_right;
            result.emplace_back(new_b_left,
                                new_b_right,
                                b_feistel,
                                counters_left,
                                counters_right,
                                ltd,
                                new_absorbed_left,
                                new_absorbed_right,
                                new_absorption_indices);
        }
    return result;
}


  
// returns all the possible trails after 1 step of the permutation
// is applied


std::vector<TruncatedTrail> TruncatedTrail::apply_step()
#if DIFFERENTIAL == 1
{
    std::vector<TruncatedTrail> result;
    // we loop over all possible outputs of the L function and
    // generate one successor for each
    std::vector<Word> all_L = possible_L_output(b_left.back()) ;
    for (Word & L : all_L)
    {
        for (Word & new_left : XOR_combine(L, b_right.back()))
        {
            new_left = rotate_word(new_left, ROTATION_AMOUNT); // rotation
            // initializing content of successor
            std::vector<Word>
                new_b_left   (b_left.begin(),    b_left.end()   ),
                new_b_right  (b_right.begin(),   b_right.end()  ),
                new_b_feistel(b_feistel.begin(), b_feistel.end());
            new_b_right.push_back(b_left.back());
            new_b_left.push_back(new_left);
            new_b_feistel.push_back(L);
            std::array<unsigned int, B> 
                new_counters_left {{0}},
                new_counters_right{{0}};
            std::map<unsigned int, unsigned int> new_ltd(ltd);
            // handling LTs on the arriving in the right
            for(unsigned int i=0; i<B; i++)
            {
                if (b_right.back()[i] == 1)
                    // case where a difference entered the step on the right
                {
                    if (L[i] == 1)
                        // lt cutting
                    {
                        new_ltd[counters_right[i]+1] += 1;
                        new_counters_left[(B+i-ROTATION_AMOUNT) % B] = 0;
                    }
                    else
                        // lt propagation
                        new_counters_left[(B+i-ROTATION_AMOUNT) % B] = counters_right[i]+1;
                }
            }
            // handling LTs arriving in the left
            for(unsigned int i=0; i<B; i++)
            {
                if (b_left.back()[i] == 0)
                    new_counters_right[i] = 0;
                else
                    new_counters_right[i] = counters_left[i]+1 ;
            }
            // adding trail to the output
            result.emplace_back(
                new_b_left,
                new_b_right,
                new_b_feistel,
                new_counters_left,
                new_counters_right,
                new_ltd,
                absorbed_left,
                absorbed_right,
                absorption_indices
                );
        }
    }
    return result;
}
#else  // If we are looking at LINEAR trails then the definition of
       // the linear layer changes; it is replaced by the transpose of
       // its inverse.
{
    std::vector<TruncatedTrail> result;
    // we loop over all possible outputs of the L function and
    // generate one successor for each
    std::vector<Word> all_L = possible_L_output(b_right.back()) ;
    for (Word & L : all_L)
    {
        for (Word & new_right : XOR_combine(L, b_left.back()))
        {
            // initializing content of successor
            std::vector<Word>
                new_b_left   (b_left.begin(),    b_left.end()   ),
                new_b_right  (b_right.begin(),   b_right.end()  ),
                new_b_feistel(b_feistel.begin(), b_feistel.end());
            new_b_right.push_back(new_right);
            new_b_left.push_back(rotate_word(b_right.back(), ROTATION_AMOUNT));
            new_b_feistel.push_back(L);
            std::array<unsigned int, B> 
                new_counters_left {{0}},
                new_counters_right{{0}};
            std::map<unsigned int, unsigned int> new_ltd(ltd);
            // handling LTs on the arriving in the right
            for(unsigned int i=0; i<B; i++)
            {
                if (b_left.back()[i] == 1)
                    // case where a difference entered the step on the left
                {
                    if (L[i] == 1)
                        // lt cutting
                    {
                        new_ltd[counters_left[i]+1] += 1;
                        new_counters_right[i] = 0;
                    }
                    else
                        // lt propagation
                        new_counters_right[i] = counters_left[i]+1;
                }
            }
            // handling LTs arriving in the right
            for(unsigned int i=0; i<B; i++)
            {
                if (b_right.back()[i] == 0)
                    new_counters_left[(B+i-ROTATION_AMOUNT) % B] = 0;
                else
                    new_counters_left[(B+i-ROTATION_AMOUNT) % B] = counters_right[i]+1 ;
            }
            // adding trail to the output
            result.emplace_back(
                new_b_left,
                new_b_right,
                new_b_feistel,
                new_counters_left,
                new_counters_right,
                new_ltd,
                absorbed_left,
                absorbed_right,
                absorption_indices
                );
        }
    }
    return result;
}
#endif



// Returns true if both left and right are fully inactive in the
// last step.
bool TruncatedTrail::finished()
{
    for(unsigned int i=0; i<B; i++)
    {
        if (b_left.back()[i] == 1)
            return false;
        if (b_right.back()[i] == 1)
            return false;
    }
    return true;
}


bool TruncatedTrail::is_empty()
{
    if (b_left.size() > 1)
        return false;
    else
    {
        for(unsigned int i=0; i<B; i++)
            if (b_left[0][i] == 1 or b_right[0][i] == 1)
                return false;
    }
    return true;
}


Word TruncatedTrail::left_branches()
{
    return b_left.back();
}


Word TruncatedTrail::right_branches()
{
    return b_right.back();
}
    

// pretty prints the truncated trail 
void TruncatedTrail::print()
{
    unsigned int a = 0;
    for(unsigned int r=0; r<size(); r++)
    {
        if (a < absorption_indices.size() and absorption_indices[a] == r)
            // case of a message absorption
        {
            for(unsigned int i=0; i<8*B+6; i++)
                std::cout << " ";
            std::cout << "<--- ";
            print_word(absorbed_left[a]);
            print_word(absorbed_right[a]);
            std::cout << std::endl;
            a++;            
        }
        std::cout << std::setw(2) << r << "  ";
        print_word(b_left[r]);
        for(unsigned int i=0; i<2*B+2; i++)
            std::cout<< " ";
        print_word(b_right[r]);
        std::cout << std::endl ;
        for(unsigned int i=0; i<2*B+6; i++)
            std::cout<< " ";
        print_word(b_feistel[r]);
        std::cout << std::endl;
    }
    std::cout << std::setw(2) << size() << "  ";
    print_word(b_left.back());
    for(unsigned int i=0; i<2*B+2; i++)
        std::cout<< " ";
    print_word(b_right.back());
    std::cout << "\nLTD: " ;
    for (auto & l : ltd)
        if (l.first != 0)
        {
            std::cout << l.first << ": " << l.second << "; ";
        }
    std::cout << std::endl ;
}

