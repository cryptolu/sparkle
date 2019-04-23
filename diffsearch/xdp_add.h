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
 * Estimation of "add" operator diff. probabilities
 *
 ******************************************************************************/

#ifndef __XDP_ADD_H__
#define __XDP_ADD_H__

#include <cstdint>
#include "defs.h"
#include <cassert>

inline WORD_T eq(const WORD_T x, const WORD_T y, const WORD_T z)
{
	WORD_T e = MASK & ~((x ^ y) | (x ^ z));
	return e;
}


inline int xdp_add_lm(const WORD_T da, const WORD_T db, const WORD_T dc)
{
	int p;
	const WORD_T eq_d = eq(da, db, dc);
	const WORD_T eq_d_sl_1 = ((eq_d << 1) | (WORD_T)1) & MASK;
	const WORD_T b_is_possible_if_zero = (eq_d_sl_1 & (da ^ db ^ dc ^ (da << 1)));
	if (b_is_possible_if_zero == 0)
	{
		const WORD_T neq = ~eq_d & MASK_NO_MSB;     // positions at which da,db and dc are not equal
		#if (WORD_SIZE <= 32)
		const int w = __builtin_popcount(neq);
		#else
		const int w = __builtin_popcountll(neq);
		#endif
		p = -w;
	}
	else
	{
		p = LOG0;
	}
	return p;
}


inline int xdp_add_lm(const WORD_T da, const WORD_T db, const WORD_T dc, unsigned int word_size)
{
	assert (word_size != 0);
	int p;
	if (word_size > 1)
	{
		const WORD_T mask = GEN_MASK(word_size);
		const WORD_T eq_d = eq(da, db, dc);
		const WORD_T eq_d_sl_1 = ((eq_d << 1) | (WORD_T)1UL) & mask;
		const WORD_T b_is_possible_if_zero = (eq_d_sl_1 & (da ^ db ^ dc ^ (da << 1)));
		if (b_is_possible_if_zero == 0)
		{
			const WORD_T neq = ~eq_d & (mask >> 1);   // positions at which da,db and dc are not equal
			#if (WORD_SIZE <= 32)
			const int w = __builtin_popcount(neq);
			#else
			const int w = __builtin_popcountll(neq);
			#endif
			p = -w;
		}
		else
		{
			p = LOG0;
		}
	}
	else
	{
		if (((da ^ db ^ dc) & (WORD_T)1) == 0)
		{
			p = 0;
		}
		else
		{
			p = LOG0; // prob = 0!
		}
	}
	return p;
}

#endif
