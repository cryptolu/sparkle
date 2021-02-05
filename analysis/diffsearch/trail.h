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

#ifndef __TRAIL_H__
#define __TRAIL_H__

#include "defs.h"
#include <cstring>
#include <array>


typedef struct
{
    WORD_T dx;
    WORD_T dy;
    WORD_T dz;
    int p;
    int cp; /* cumulated probability */
} Differential;

typedef std::array<Differential, 100> Trail;

/* iroundm1 is iround - 1 */
inline void pushTrail(Trail &trail, Differential differential, unsigned int iroundm1)
{
    std::memcpy(
        &(trail[iroundm1]),
        &differential,
        sizeof(Differential)
    );
}

/* iroundm1 is iround - 1 */
inline void popTrail(Trail &trail, unsigned int iroundm1)
{
    trail[iroundm1].p = LOG0;
}

void fprintTrail(FILE *fh, Trail trail);
void fprintDifferential(FILE *fh, Trail trail);
void clearTrail(Trail &trail);

#endif
