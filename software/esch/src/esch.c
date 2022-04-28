///////////////////////////////////////////////////////////////////////////////
// esch.c: Optimized C implementation of the hash function ESCH.             //
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


// This source code file should be compiled with the following set of flags:
// -std=c99 -Wall -Wextra -Wshadow -fsanitize=address,undefined -O2


// gencat_hash.c shall be used to generate the test vector output file. The
// test vector output file shall be provided in the corresponding 
// crypto_hash/[algorithm]/ directory


#include <string.h>  // for memcpy, memset
#include "esch.h"
#include "sparkle.h"


///////////////////////////////////////////////////////////////////////////////
//////////// SPARKLE CONFIGURATIONS FOR THE TWO INSTANCES OF ESCH /////////////
///////////////////////////////////////////////////////////////////////////////

#if (ESCH_INST == 256)
#define STATE_BYTES     48
#define RATE_BYTES      16
#define STEPS_SLIM       7
#define STEPS_BIG       11

#elif (ESCH_INST == 384)
#define STATE_BYTES     64
#define RATE_BYTES      16
#define STEPS_SLIM       8
#define STEPS_BIG       12

#else
#error "Invalid definition of ESCH instance!"
#endif

#define ESCH_DIGEST_WORDS (ESCH_DIGEST_BYTES/4)

#define STATE_BRANS (STATE_BYTES/8)
#define STATE_WORDS (STATE_BYTES/4)
#define RATE_BRANS  (RATE_BYTES/8)
#define RATE_WORDS  (RATE_BYTES/4)
#define CAP_BYTES   (STATE_BYTES-RATE_BYTES)
#define CAP_BRANS   (CAP_BYTES/8)
#define CAP_WORDS   (CAP_BYTES/4)

#define CONST_M1 (((uint32_t) 1) << 24)
#define CONST_M2 (((uint32_t) 2) << 24)


///////////////////////////////////////////////////////////////////////////////
//// PREPROCESSOR DIRECTIVES TO REPLACE THE C CODE OF SPARKLE BY ASM CODE /////
///////////////////////////////////////////////////////////////////////////////


// When this file is compiled for an AVR microcontroller and SPARKLE_ASSEMBLER
// is defined (see esch.h), then the AVR assembler implementation of the
// SPARKLE permutation is used. On the other hand, if SPARKLE_ASSEMBLER is not
// defined, then the C version (i.e. the function `sparkle`) is used.

#if (defined(__AVR) || defined(__AVR__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_avr(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_avr((state), (brans), (steps))
#endif // if (defined(__AVR__) || ...


// When this file is compiled for a MSP430 (or a MSP430X) microcontroller and
// SPARKLE_ASSEMBLER is defined (see esch.h), then the MSP430 assembler
// implementation of the SPARKLE permutation is used. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// `sparkle`) is used.

#if (defined(MSP430) || defined(__MSP430__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_msp(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_msp((state), (brans), (steps))
#endif // if (defined(MSP430) || ...


// When this file is compiled for a 32-bit RISC-V microcontroller (e.g. RV32I)
// and SPARKLE_ASSEMBLER is defined (see esch.h), then one of the two
// branch-unrolled RV32 assembler implementations of the SPARKLE permutation is
// used, depending on the concrete ESCH instance. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// `sparkle`) is used.

#if defined(__riscv_xlen) && (__riscv_xlen == 32) && defined(SPARKLE_ASSEMBLER)
#if (STATE_BYTES == 48)
extern void sparkle384_rv32(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle384_rv32((state), (steps))
#elif (STATE_BYTES == 64)
extern void sparkle512_rv32(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle512_rv32((state), (steps))
#endif // if (STATE_BYTES == 48)
#endif // if defined(__riscv_xlen) && ...


// When this file is compiled for an ARM microcontroller and SPARKLE_ASSEMBLER
// is defined (see esch.h), then one of the two branch-unrolled ARM assembler
// implementations of the SPARKLE permutation is used, depending on the
// concrete ESCH instance. On the other hand, if SPARKLE_ASSEMBLER is not
// defined, then the C version (i.e. the function `sparkle`) is used.

#if (defined(__arm__) || defined(_M_ARM)) && defined(SPARKLE_ASSEMBLER)
#if (STATE_BYTES == 48)
extern void sparkle384_arm(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle384_arm((state), (steps))
#elif (STATE_BYTES == 64)
extern void sparkle512_arm(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle512_arm((state), (steps))
#endif // if (STATE_BYTES == 48)
#endif // if (defined(__arm__) || ...


// In all ARMv7 architectures, including ARMv7-M, unaligned data access is
// available (and is the default), whereas ARMv6-M (e.g. Cortex-M0/M0+) and
// ARMv8-M baseline (e.g. Cortex-M23) require strict alignment, which means
// the address of a 32-bit integer must be a multiple of four. The following
// preprocessor directives set the identifier ALIGN_OF_UI32 to 1 (see below
// for an explanation) when this file is compiled for an ARMv7 device.

#ifndef ALIGN_OF_UI32
#if defined(__arm__) || defined(_M_ARM)
#if ((__ARM_ARCH == 7) && (__ARM_ARCH_ISA_THUMB == 2)) ||     \
  ((__TARGET_ARCH_ARM == 7) && (__TARGET_ARCH_THUMB == 4)) || \
  ((__TARGET_ARCH_ARM == 0) && (__TARGET_ARCH_THUMB == 4))
#define ALIGN_OF_UI32 1
#endif // if ((__ARM_ARCH == 7) && ..
#endif // if defined(__arm__) || ...
#endif // ifndef ALIGN_OF_UI32


///////////////////////////////////////////////////////////////////////////////
/////// HELPER FUNCTIONS AND MACROS (INJECTION OF MESSAGE BLOCK, ETC.) ////////
///////////////////////////////////////////////////////////////////////////////


#define ROT(x, n) (((x) >> (n)) | ((x) << (32-(n))))
#define ELL(x) (ROT(((x) ^ ((x) << 16)), 16))


// The high-level Hash API specifies that the message to be hashed is stored in
// an array of type unsigned char. However, the SPARKLE permutation operates on
// 32-bit words and performs best when the data to be processed is stored in an
// uint32_t-array. Casting an unsigned-char pointer to an uint32_t-pointer
// increases the alignment requirements, i.e. the base address of the array
// must be even on 16-bit architectures and a multiple of four (i.e. 4-byte
// aligned) on 32-bit and 64-bit platforms. The preprocessor statements below
// can be used to determine the alignment requirements an unsigned-char-pointer
// has to meet to permit casting to an uint32_t-pointer.

#ifndef ALIGN_OF_UI32
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)  // C11
#include <stdalign.h>
#define ALIGN_OF_UI32 alignof(uint32_t)
#else  // C11 standard is not available
#define MIN_SIZE(a, b) ((sizeof(a) < sizeof(b)) ? sizeof(a) : sizeof(b))
#define ALIGN_OF_UI32 MIN_SIZE(uint32_t, uint_fast8_t)  // stdint.h
#endif // if defined(__STDC_VERSION__) && ...
#endif // ifndef ALIGN_OF_UI32


// Injection of a 16-byte block of the message to the state. According to the
// specification, the Feistel function is performed on a message block that is
// padded with 0-bytes to reach a length of STATE_BYTES/2 bytes (i.e. 24 bytes
// for ESCH256, 32 bytes for ESCH384). However, this padding can be omitted by
// adapting the Feistel function accordingly. The third parameter indicates
// whether the uint8_t-pointer `in` is properly aligned to permit casting to a
// uint32_t-pointer. If this is the case then array `in` is processed directly,
// otherwise it is first copied to an aligned buffer. 

static void add_msg_blk(uint32_t *state, const uint8_t *in, int aligned)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t *in32;
  uint32_t tmpx = 0, tmpy = 0;
  int i;
  
  if (aligned) {  // `in` can be casted to uint32_t pointer
    in32 = (uint32_t *) in;
  } else {  // `in` is not sufficiently aligned for casting
    memcpy(buffer, in, RATE_BYTES);
    in32 = (uint32_t *) buffer;
  }
  
  for(i = 0; i < RATE_WORDS; i += 2) {
    tmpx ^= in32[i];
    tmpy ^= in32[i+1];
  }
  tmpx = ELL(tmpx);
  tmpy = ELL(tmpy);
  for(i = 0; i < RATE_WORDS; i += 2) {
    state[i] ^= (in32[i] ^ tmpy);
    state[i+1] ^= (in32[i+1] ^ tmpx);
  }
  for(i = RATE_WORDS; i < (STATE_WORDS/2); i += 2) {
    state[i] ^= tmpy;
    state[i+1] ^= tmpx;
  }
}


// Injection of the last message block to the state. Since this last block may
// require padding, it is always copied to a buffer.

static void add_msg_blk_last(uint32_t *state, const uint8_t *in, size_t inlen)
{
  uint32_t buffer[RATE_WORDS];
  uint8_t *bufptr;
  uint32_t tmpx = 0, tmpy = 0;
  int i;
  
  memcpy(buffer, in, inlen);
  if (inlen < RATE_BYTES) {  // padding
    bufptr = ((uint8_t *) buffer) + inlen;
    memset(bufptr, 0, (RATE_BYTES - inlen));
    *bufptr = 0x80;
  }
  
  for(i = 0; i < RATE_WORDS; i += 2) {
    tmpx ^= buffer[i];
    tmpy ^= buffer[i+1];
  }
  tmpx = ELL(tmpx);
  tmpy = ELL(tmpy);
  for(i = 0; i < RATE_WORDS; i += 2) {
    state[i] ^= (buffer[i] ^ tmpy);
    state[i+1] ^= (buffer[i+1] ^ tmpx);
  }
  for(i = RATE_WORDS; i < (STATE_WORDS/2); i += 2) {
    state[i] ^= tmpy;
    state[i+1] ^= tmpx;
  }
}


///////////////////////////////////////////////////////////////////////////////
///////////// LOW-LEVEL HASH FUNCTIONS (FOR USE WITH FELICS-HASH) /////////////
///////////////////////////////////////////////////////////////////////////////


// The Initialize function sets all branches of the state to 0.

void Initialize(uint32_t *state)
{
  int i;
  
  for (i = 0; i < STATE_WORDS; i++) {
    state[i] = 0;
  }
}


// The ProcessMessage function absorbs the message into the state (in blocks of
// 16 bytes). According to the specification, the constant Const_M is first
// transformed via the inverse Feistel function, added to the (padded) message
// block, and finally injected to the state via the Feistel function. Since the
// Feistel function and the inverse Feistel function cancel out, we can simply
// inject the constant directly to the state.

void ProcessMessage(uint32_t *state, const uint8_t *in, size_t inlen)
{
  // check whether `in` can be casted to uint32_t pointer (we can use here
  // size_t instead of uintptr_t since ALIGN_OF_UI32 is either 1, 2, or 4)
  int aligned = ((size_t) in) % ALIGN_OF_UI32 == 0;
  // printf("Address of `in`: %p\n", in);
  
  // Main Hashing Loop
  
  while (inlen > RATE_BYTES) {
    // addition of a message block to the state
    add_msg_blk(state, in, aligned);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Hashing of Last Block
  
  // addition of constant M1 or M2 to the state
  state[STATE_BRANS-1] ^= ((inlen < RATE_BYTES) ? CONST_M1 : CONST_M2);
  // addition of last msg block (incl. padding)
  add_msg_blk_last(state, in, inlen);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


// The Finalize function generates the message digest by "squeezing" (i.e. by
// calling SPARKLE with a slim number of steps) until the digest has reached a
// byte-length of ESCH_DIGEST_BYTES.

void Finalize(uint32_t *state, uint8_t *out)
{
  size_t outlen;
  
  memcpy(out, state, RATE_BYTES);
  outlen = RATE_BYTES;
  out += RATE_BYTES;
  while (outlen < ESCH_DIGEST_BYTES) {
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    memcpy(out, state, RATE_BYTES);
    outlen += RATE_BYTES;
    out += RATE_BYTES;
  }
}


///////////////////////////////////////////////////////////////////////////////
////////////// HIGH-LEVEL HASH FUNCTIONS (FOR USE WITH SUPERCOP) //////////////
///////////////////////////////////////////////////////////////////////////////


// To ensure compatibility with the SUPERCOP, the below implementation of 
// crypto_hash can handle overlapping input and output buffers.

int crypto_hash(UChar *out, const UChar *in, ULLInt inlen)
{
  uint32_t state[STATE_WORDS];
  size_t insize = (size_t) inlen;
  
  Initialize(state);
  ProcessMessage(state, in, insize);
  Finalize(state, out);
  
  return 0;
}
