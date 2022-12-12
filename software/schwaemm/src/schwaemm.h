///////////////////////////////////////////////////////////////////////////////
// schwaemm.h: Optimized C implementation of the AEAD algorithm SCHWAEMM.    //
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

#ifndef SCHWAEMM_H
#define SCHWAEMM_H

#include <stddef.h>  // for size_t
#include <stdint.h>  // for uint8_t, uint32_t

typedef unsigned char UChar;
typedef unsigned long long int ULLInt;

// Define the SCHWAEMM instance here (api.h has to match when using SUPERCOP!).
// The main instance is SCHWAEMM256_128, which has a block size of 256 bits and
// a key size of 128 bits. The other instances of SCHWAEMM are SCHWAEMM128_128,
// SCHWAEMM192_192, and SCHWAEMM256_256. Hence, valid values for SCHWAEMM_INST
// are 256128, 128128, 192192, and 256256.

#ifndef SCHWAEMM_INST
#define SCHWAEMM_INST 256128
#endif

// The identifier SPARKLE_ASSEMBLER determines whether the low-level functions 
// in schwaemm.c use the C implementation or an assembler implementation of the
// SPARKLE permutation. Currently, assembler code for SPARKLE exists for AVR,
// ARMv6-M, and ARMv7-M.

// #define SPARKLE_ASSEMBLER

// Key-size, nonce-size, and tag-size (in bytes) of all instances of SCHWAEMM.

#if (SCHWAEMM_INST == 128128)
#define SCHWAEMM_KEY_BYTES      16
#define SCHWAEMM_NONCE_BYTES    16
#define SCHWAEMM_TAG_BYTES      16
#elif (SCHWAEMM_INST == 256128)
#define SCHWAEMM_KEY_BYTES      16
#define SCHWAEMM_NONCE_BYTES    32
#define SCHWAEMM_TAG_BYTES      16
#elif (SCHWAEMM_INST == 192192)
#define SCHWAEMM_KEY_BYTES      24
#define SCHWAEMM_NONCE_BYTES    24
#define SCHWAEMM_TAG_BYTES      24
#elif (SCHWAEMM_INST == 256256)
#define SCHWAEMM_KEY_BYTES      32
#define SCHWAEMM_NONCE_BYTES    32
#define SCHWAEMM_TAG_BYTES      32
#endif

// Prototypes of the low-level functions (for benchmarking with FELICS-AEAD).

void Initialize(uint32_t *state, const uint8_t *key, const uint8_t *nonce);
void ProcessAssocData(uint32_t *state, const uint8_t *in, size_t inlen);
void ProcessPlainText(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen);
void ProcessCipherText(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen);
void Finalize(uint32_t *state, const uint8_t *key);
void GenerateTag(uint32_t *state, uint8_t *tag);
int VerifyTag(uint32_t *state, const uint8_t *tag);

// Prototypes of the high-level functions (for benchmarking with SUPERCOP).

int crypto_aead_encrypt(UChar *c, ULLInt *clen, const UChar *m, ULLInt mlen, \
  const UChar *ad, ULLInt adlen, const UChar *nsec, const UChar *npub,       \
  const UChar *k);
int crypto_aead_decrypt(UChar *m, ULLInt *mlen, UChar *nsec, const UChar *c, \
  ULLInt clen, const UChar *ad, ULLInt adlen, const UChar *npub,             \
  const UChar *k);

#endif  // SCHWAEMM_H
