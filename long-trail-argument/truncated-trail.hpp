/** Time-stamp: <2019-03-11 17:28:16 leo>
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


#ifndef TRUNCATED_TRAIL_HPP
#define TRUNCATED_TRAIL_HPP

#include "trail-search.hpp"



class TruncatedTrail 
{
private:
    std::vector<Word> b_left;  // branches going in the Feistel function
    // at round i
    std::vector<Word> b_right; // branches receiving the Feistel function
    // during round i
    std::vector<Word> b_feistel;  // outputs of the Feistel function at
    // round i

    std::vector<unsigned int> absorption_indices; // the steps at which
    // absorptions occured
    std::vector<Word> absorbed_left;   // differences absorbed on the left
    std::vector<Word> absorbed_right;  // differences absorbed on the right

    std::array <unsigned int, B> counters_left;   // counters for the LTS
    std::array <unsigned int, B>  counters_right; // counters for the LTS
    std::map<unsigned int, unsigned int> ltd; // current lt decomposition

public:

    TruncatedTrail(
        std::vector<Word> _b_left, 
        std::vector<Word> _b_right,
        std::vector<Word> _b_feistel,
        std::array<unsigned int, B> _counters_left,
        std::array<unsigned int, B> _counters_right,
        std::map<unsigned int, unsigned int> _ltd,
        std::vector<Word> _absorbed_left,
        std::vector<Word> _absorbed_right,
        std::vector<unsigned int> _absorption_indices
        ) :
        b_left(_b_left),
        b_right(_b_right),
        b_feistel(_b_feistel),
        absorption_indices(_absorption_indices),
        absorbed_left(_absorbed_left),
        absorbed_right(_absorbed_right),
        counters_left(_counters_left),
        counters_right(_counters_right),
        ltd(_ltd) {} ;
  

    TruncatedTrail(
        const Word left_branches, 
        const Word right_branches
        ) :
        TruncatedTrail(
            std::vector<Word>(1, left_branches),
            std::vector<Word>(1, right_branches),
            std::vector<Word>(0),
            std::array<unsigned int, B>{{0}},
            std::array<unsigned int, B>{{0}},
            std::map<unsigned int, unsigned int>(),
            std::vector<Word>(0),
            std::vector<Word>(0),
            std::vector<unsigned int>(0)
            ) {} ;


    // returns the number of steps in the trail considered
    size_t size();
  
    // returns the score of the current trail 
    unsigned int score(
        unsigned int r_p_s,                  // rounds per step
        FullProba pr // probabilities
        );
  
    // returns all possible results of absorbing a_left and a_right
    std::vector<TruncatedTrail> absorb(Word a_left, Word a_right);
  
    // returns all the possible trails after 1 step of the permutation
    // is applied
    std::vector<TruncatedTrail> apply_step();

    // Returns true if both left and right are fully inactive in the
    // last step.
    bool finished();

    bool is_empty();

    Word left_branches();
    Word right_branches();
    
    // pretty prints the truncated trail 
    void print();
};




#endif
