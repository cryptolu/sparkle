/** Time-stamp: <2019-04-18 08:34:33 leo>
 *
 * Implementation of the truncated trail search and the long trail
 * argument for SPARX-like step functions.
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

// Intializing global variables

#if (DIFFERENTIAL == 1)         // rounds:  1 2 3 4  5  6  7  8
std::array<unsigned int, N_KNOWN> STARTING_PROBA{{0,1,2,6,10,18,24,32,36,42,46,52}};
#else
std::array<unsigned int, N_KNOWN> STARTING_PROBA{{0,0,1,2, 5, 8,13,17}};
#endif

FullProba PROBA{{0}} ;

// below are estimated probabilities for Matsui's search. Not really needed if the number of steps is low.

#if USE_MATSUI == 1
// DIFFERENTIAL estimates           0 1   2   3    4    5    6    7    8    9   10   11   12
#if   (DIFFERENTIAL == 1) && (B==2)
std::vector<unsigned int> ESTIMATE{{0,5, 25, 40,  80, 100, 120, 160, 180, 220, 240, 280, 300}};
#elif (DIFFERENTIAL == 1) && (B==3)
std::vector<unsigned int> ESTIMATE{{0,5, 25, 65,  95, 170, 195, 220, 255, 255,255, 255, 255}};
#elif (DIFFERENTIAL == 1) && (B==4)
std::vector<unsigned int> ESTIMATE{{0,5, 30, 70, 110, 200, 230, 260, 270, 290, 420, 430, 490, 490}};

// LINEAR estimates                 0 1  2  3  4   5   6   7   8   9  10  11  12
#elif (DIFFERENTIAL == 0) && (B==2) 
std::vector<unsigned int> ESTIMATE{{0,2,15,20,40, 50, 70, 90,100,100,100,100,100}};
#elif (DIFFERENTIAL == 0) && (B==3)
std::vector<unsigned int> ESTIMATE{{0,2,13,20,40,60,75,95,115,140,150,150,150,150,150}};
#elif (DIFFERENTIAL == 0) && (B==4)// 1  2  3  4  5   6   7   8   9  10  11  12
std::vector<unsigned int> ESTIMATE{{0,2,15,25,45,90,100,125,150,190,205,225,250,250,250,250}};
#endif

#else
std::vector<unsigned int> ESTIMATE{{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
#endif



// !SECTION! Studying simple trails
// ================================

TruncatedTrail best_extension(
    TruncatedTrail t,
    unsigned int n_steps,
    unsigned int target_final_proba,
    unsigned int r_p_s,
    FullProba pr
    )
{
    if (n_steps == 0)
        return t;
    else 
    {
        TruncatedTrail result(Word{{0}}, Word{{0}});
        if (t.score(r_p_s, pr) + ESTIMATE[n_steps] <= target_final_proba and not t.is_empty())
        {
            unsigned int best_score = target_final_proba;
            for (TruncatedTrail & s : t.apply_step())
            {
                unsigned int local_score = s.score(r_p_s, pr);
                if ((local_score > 0) and (local_score+ESTIMATE[n_steps-1] <= best_score))
                {
                    TruncatedTrail candidate = best_extension(s, 
                                                              n_steps - 1, 
                                                              best_score,
                                                              r_p_s,
                                                              pr);
                    unsigned int candidate_score = candidate.score(r_p_s, pr);
                    if (candidate_score > 0 and candidate_score < best_score)
                    {
                        best_score = candidate_score;
                        result = candidate;
                    }
                }
            }
        }
        return result;
    }
}


void print_best_simple_trail(unsigned int target_proba,
                             unsigned int n_steps,
                             unsigned int r_p_s,
                             FullProba pr)
{
    TruncatedTrail best_trail(Word{{0}}, Word{{0}});
    unsigned int best_score = target_proba;
    std::vector<unsigned int> unique = unique_up_to_rotation();
    for(unsigned int & c_r : unique)
        for(unsigned int c_l=0; c_l<B_SPACE; c_l++)
            if (c_l != 0 or c_r != 0)
            {
                Word 
                    w_l = word_from_int(c_l),
                    w_r = word_from_int(c_r);
                // print_word(w_l);
                // print_word(w_r);
                // std::cout << std::endl;
                TruncatedTrail t = best_extension(
                    TruncatedTrail(w_l, w_r),
                    n_steps,
                    best_score,
                    r_p_s,
                    pr
                    );
                if (not t.is_empty())
                {
                    best_score = t.score(r_p_s, pr) ;
                    best_trail = t;
                }          
            }  
    std::cout << n_steps << "\t" << best_score << std::endl;
    
#if PRINT_TRAILS == 1
    best_trail.print();
    std::cout << "\n----------------------\n" << std::endl;
#endif
}



// !SECTION! Studying absorbed trails
// ==================================


bool is_output_cancelled(TruncatedTrail t)
{
#if RATE_TO_RATE == 1
  for(unsigned int i=R_LEFT; i<B; i++)
    if (t.left_branches()[i] != 0)
      return false;
  for(unsigned int i=R_RIGHT; i<B; i++)
    if (t.right_branches()[i] != 0)
      return false;
#endif
  return true;
}



TruncatedTrail best_extension_that_cancels(
    TruncatedTrail t, 
    unsigned int n_steps,
    unsigned int threshold,
    unsigned int r_p_s,
    FullProba pr)
{
    if (n_steps == 0)
    {
        unsigned int score = t.score(r_p_s, pr);
        if (
            (not t.is_empty())
            and (score > 0)
            and is_output_cancelled(t)
            and (score <= threshold)
            )
            return t;
        else
            return TruncatedTrail(Word{{0}}, Word{{0}});
    }
    else
    {
        TruncatedTrail result(Word{{0}}, Word{{0}});
        unsigned int best_score = threshold;
        for (TruncatedTrail & s : t.apply_step())
        {
            unsigned int local_score = s.score(r_p_s, pr);
            if (
                (local_score > 0) 
                and (local_score + ESTIMATE[n_steps-1] <= best_score)
                )
            {
                TruncatedTrail c = best_extension_that_cancels(s, 
                                                               n_steps - 1, 
                                                               threshold,
                                                               r_p_s,
                                                               pr);
                {
                    unsigned int candidate_score = c.score(r_p_s, pr);
                    if (
                        (not c.is_empty()) 
                        and (candidate_score > 0)
                        and (candidate_score < best_score)
                        )
                    {
                        result = c;
                        best_score = candidate_score;
                    }
                }
            }
        }
        return result;
    }
}


void print_best_absorption(unsigned int target_proba,
                           unsigned int r_p_s,
                           unsigned int s_p_a,
                           FullProba pr)
{
    TruncatedTrail best_trail(Word{{0}}, Word{{0}});  
    unsigned int best_score = target_proba;
    
    for(unsigned int b=0; b<R_RIGHT_SPACE; b++)
        for(unsigned int a=0; a<R_LEFT_SPACE; a++)
            if ((a != 0) or (b != 0)) 
            {
                Word 
                    a_left = word_from_int(a), 
                    a_right= word_from_int(b);
                TruncatedTrail t(a_left, a_right);
                TruncatedTrail u = best_extension_that_cancels(t,
                                                               s_p_a,
                                                               best_score,
                                                               r_p_s,
                                                               pr);
                unsigned int local_score = u.score(r_p_s, pr);
                if (
                    (not u.is_empty()) 
                    and (local_score > 0)
                    and (local_score <= best_score)
                    )
                {
                  best_score = local_score;
                  best_trail = u;
                  // std::cout << "--- " << best_score << "\n";
                  // t.print();
                }
            }
    std::cout << s_p_a << "\t"
              << best_score << "\t"
              << ESTIMATE[s_p_a] << "\t"
              << target_proba << "\t"
              << best_trail.score(r_p_s, pr)
              << std::endl;
#if PRINT_TRAILS == 1
    std::cout << "\n" << std::endl;
    best_trail.print();
    std::cout << "\n----------------------\n" << std::endl;
#endif
}


// !SECTION! Main Function
// =======================

void print_proba()
{
    for(unsigned int i=0; i<PROBA.size(); i++)
        std::cout << PROBA[i] << ", ";
    std::cout << std::endl;
}



int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "Usage: ./Trail-Search <type> <min steps> <max steps>\n"
                  << "- <type> must be either 'permutation' or 'absorption'\n"
                  << "- <min steps> is the minimum number of steps considered\n"
                  << "- <max steps> is the maximum number of steps considered\n"
                  << std::endl;
        exit(0);
    }
    else
    {
        if (USE_MATSUI == 1)
            std::cout << "using MATSUI" << std::endl;
        else
            std::cout << "*NOT* using MATSUI" << std::endl;
        if (DIFFERENTIAL == 1)
            std::cout << "DIFFERENTIAL" << std::endl;
        else
            std::cout << "LINEAR" << std::endl;

        generate_PROBA();

        // for(unsigned int i=0; i<PROBA.size(); i+=4)
        //   std::cout << i << "\t" << PROBA[i] << std::endl;

        std::cout << "B=" << B
                  << "  R=" << R_LEFT + R_RIGHT
                  << "=(" << R_LEFT << "+" << R_RIGHT << ")"
                  << "  Rotation=" << ROTATION_AMOUNT 
                  << std::endl;

        unsigned int
            min_steps = atoi(argv[2]),
            max_steps = atoi(argv[3]),
            r_p_s = 4;
        if (std::string(argv[1]) == "permutation")
        {
            std::cout << "PERMUTATION" << std::endl;
            for(unsigned int j=min_steps; j<max_steps; j++)
            {
                unsigned int security = 64*(2*B);
#if DIFFERENTIAL == 0
                security = 32*2*B;
#endif
                print_best_simple_trail(security, j, r_p_s, PROBA) ;
            }
        }
        else
        {
            std::cout << "ABSORPTION" << std::endl;
            if (RATE_TO_RATE == 1)
                std::cout << "rate to rate" << std::endl;
            else
                std::cout << "rate to anything" << std::endl;
            for (unsigned int s_p_a=min_steps; s_p_a<max_steps; s_p_a++)
            {
                unsigned int target_proba = 380; //64*(2*B);
#if DIFFERENTIAL == 0
                target_proba = 32*2*B;
#endif                
                
                print_best_absorption(target_proba,
                                      r_p_s,
                                      s_p_a,
                                      PROBA);
            }
        }
    }
}
