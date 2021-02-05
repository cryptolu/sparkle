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
 * Definitions
 *
 ******************************************************************************/

#ifndef __DEFS_H__
#define __DEFS_H__

#include <cstdint>

/* undef STATS if you don't want execution statistics (for faster runs) */
//#define STATS

#define ALL_TRAILS 0 // if 0 find all trails; if 1 find single trail

#define POW(n) ((WORD_T)1 << (n))

#define WORD_SIZE 32

//#define NROUNDS 7 // vv
#ifdef __TRAIL_H__
	extern unsigned NROUNDS;
	extern uint32_t gconst_r[100];
	extern uint32_t gconst_s[100];
#else
	#ifdef TRAIL_SEARCH
		extern unsigned NROUNDS;
		extern uint32_t gconst_r[100];
		extern uint32_t gconst_s[100];
	#else
	#endif
#endif

#ifdef MAIN
	unsigned NROUNDS;
	uint32_t gconst_r[100] = {7,1,15,7,1,15,7,1,15,7,1,15};
	uint32_t gconst_s[100] = {15,7,8,15,7,8,15,7,8,15,7,8};
	int g_best_B[100] = {
	     // 0,                         // 1
	    // -1,                         // 2
	    // -2,                         // 3
	    // -6,                         // 4
	    // -10,                         // 5
	    // -18,                         // 6
	    // -18,                         // 7
	};
#else
	extern int g_best_B[100];
	extern unsigned NROUNDS;
	extern uint32_t gconst_r[100];
	extern uint32_t gconst_s[100];
#endif


#define WORD_T uint32_t

#if (WORD_SIZE != 32)
#error("WORD_SIZE must be 32")
#endif


#define MASK ((uint32_t)0xffffffffUL >> (32 - WORD_SIZE))	/* mask for the WORD_SIZE LS bits of a 32-bit word. */
#define MASK_NO_MSB ((uint32_t)0xffffffffUL >> (32 - (WORD_SIZE - 1)))

#define GEN_MASK(word_size) ((uint32_t)0xffffffffUL >> (32 - (word_size)))


#define XOR(x,y) ((x ^ y) & MASK) /* The XOR operation on words of size \ref WORD_SIZE */
#define LROT(x,r) (((x << r) | (x >> (WORD_SIZE - r))) & MASK) /* Rotate \p x by \p r positions to the left; \p x is of size \ref WORD_SIZE */
#define RROT(x,r) (((x >> r) | (x << (WORD_SIZE - r))) & MASK) /* Rotate \p x by \p r positions to the right; \p x is of size \ref WORD_SIZE */

#define LOG0 -10000

/* global array of bounds */
#define NROUNDS_MAX 100

//const uint32_t gconst_r[12] = {0,31,0,31,0,31,0,31,0,31,0,31};
//const uint32_t gconst_s[12] = {29,24,29,24,29,24,29,24,29,24,29,24};


#endif


