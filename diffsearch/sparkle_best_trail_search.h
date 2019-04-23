/*
 * Search for the optimal differential trails of the ARXbox of SPARKLE
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
/******************************************************************************
 *
 * sparkle_best_trail_search
 *
 ******************************************************************************/

#ifndef __SPARKLE_BEST_TRAIL_SEARCH_H__
#define __SPARKLE_BEST_TRAIL_SEARCH_H__

#include "defs.h"

extern int g_Bn;
extern Trail trail;
extern uint64_t nNodes;

#ifdef STATS
extern uint64_t nNodesPerRound[NROUNDS];
extern uint64_t nNodesPerRoundPerIbit[NROUNDS][WORD_SIZE + 1];
extern uint64_t nPathsPerRoundPerIbit[NROUNDS][WORD_SIZE + 1];
#endif

bool sparkle_best_trail_search_i(unsigned int iround, unsigned int ibit, WORD_T alpha, WORD_T beta, WORD_T gamma);
bool sparkle_best_trail_search(void);
bool sparkle_best_trail_search_split(WORD_T alpha_next, WORD_T beta_next, WORD_T gamma_next);

#endif
