///////////////////////////////////////////////////////////////////////////////
// esch.h: Optimized C implementation of the hash function ESCH.             //
// This file is part of the SPARKLE package that was sent to NIST during the //
// 3rd round of the Lightweight Cryptography (LWC) standardization project.  //
// Version 1.2.1 (18-Oct-21), see <http://github.com/cryptolu/> for updates. //
// Authors: The SPARKLE Group (Christof Beierle, Alex Biryukov, Luan Cardoso //
// dos Santos, Johann Groszschaedl, Amir Moradi, Leo Perrin, Aein Rezaei     //
// Shahmirzadi, Aleksei Udovenko, Vesselin Velichkov, and Qingju Wang).      //
// License: GPLv3 (see LICENSE file), other licenses available upon request. //
// Copyright (C) 2019-2021 University of Luxembourg <http://www.uni.lu/>.    //
// ------------------------------------------------------------------------- //
// This program is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by the     //
// Free Software Foundation, either version 3 of the License, or (at your    //
// option) any later version. This program is distributed in the hope that   //
// it will be useful, but WITHOUT ANY WARRANTY; without even the implied     //
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  //
// GNU General Public License for more details. You should have received a   //
// copy of the GNU General Public License along with this program. If not,   //
// see <http://www.gnu.org/licenses/>.                                       //
///////////////////////////////////////////////////////////////////////////////

#ifndef ESCH_H
#define ESCH_H

#include <stddef.h>  // for size_t
#include <stdint.h>  // for uint8_t, uint32_t

typedef unsigned char UChar;
typedef unsigned long long int ULLInt;

// Define the ESCH instance here (api.h has to match when using SUPERCOP!). The
// main instance is ESCH256, which has a digest size of 256 bits. A second
// instance of ESCH is ESCH384. Hence, valid values for ESCH_INST are 256 and
// 384.

#ifndef ESCH_INST
#define ESCH_INST 256
#endif

// The identifier SPARKLE_ASSEMBLER determines whether the low-level functions 
// in esch.c use the C implementation or an assembler implementation of the
// SPARKLE permutation. Currently, assembler code for SPARKLE exists for AVR,
// ARMv6-M, and ARMv7-M.

// #define SPARKLE_ASSEMBLER

// Digest-size (in bytes) of all instances of ESCH.

#if (ESCH_INST == 256)
#define ESCH_DIGEST_BYTES       32
#elif (ESCH_INST == 384)
#define ESCH_DIGEST_BYTES       48
#endif

// Prototypes of the low-level functions (for benchmarking with FELICS-AEAD).

void Initialize(uint32_t *state);
void ProcessMessage(uint32_t *state, const uint8_t *in, size_t inlen);
void Finalize(uint32_t *state, uint8_t *out);

// Prototypes of the high-level functions (for benchmarking with SUPERCOP).

int crypto_hash(UChar *out, const UChar *in, ULLInt inlen);

#endif  // ESCH_H
