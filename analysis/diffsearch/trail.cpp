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
 * Trail utilities
 *
 ******************************************************************************/

#include <cstdio>
#include "trail.h"
#include "defs.h"

void fprintTrail(FILE *fh, Trail trail)
{
	int p = 0;
	for (uint32_t i = 0; i < NROUNDS; i++)
	{
		#if (WORD_SIZE <=32)
		fprintf(fh, "# %2d: 0x%08x 0x%08x -> 0x%08x %+d <%+d>\n", i, trail[i].dx, trail[i].dy, trail[i].dz, trail[i].p, trail[i].cp);
		#else
		fprintf(fh, "# %2d: 0x%016lx 0x%016lx -> 0x%016lx %+d <%+d>\n", i, trail[i].dx, trail[i].dy, trail[i].dz, trail[i].p, trail[i].cp);
		#endif
		p += trail[i].p;
	}
	fprintf(fh, "# p_trail %d\n", p);
}

void fprintDifferential(FILE *fh, Trail trail)
{
	int p = 0;
	WORD_T DA[2] = {0}; // input difference
	WORD_T DB[2] = {0}; // output difference

	// input
	DA[0] = trail[0].dx; // left
	DA[1] = LROT(trail[0].dy, gconst_r[0]); // right
	// output
	DB[0] = trail[NROUNDS-1].dz;  
	DB[1] = XOR(LROT(trail[NROUNDS-1].dy,gconst_r[NROUNDS-1]),RROT(trail[NROUNDS-1].dz,gconst_s[NROUNDS-1])); // right

	// compute trail probability
	for (uint32_t i = 0; i < NROUNDS; i++)
	{
		p += trail[i].p;
	}
#if (WORD_SIZE <=32)
	fprintf(fh, "{%2d, 0x%08x, 0x%08x, 0x%08x, 0x%08x, %d},\n", NROUNDS, DA[0], DA[1], DB[0], DB[1], abs(p));
#else
	fprintf(fh, "{%2d, 0x%016lx, 0x%016lx, 0x%016lx, 0x%016lx, %d},\n", NROUNDS, DA[0], DA[1], DB[0], DB[1], abs(p));
#endif
	//	fprintf(fh, "R# p_trail %d\n", p);
}

void clearTrail(Trail &trail)
{
	for (unsigned int i = 0; i < NROUNDS; i++)
	{
		trail[i] = {0, 0, 0, LOG0, 0};
	}
}
