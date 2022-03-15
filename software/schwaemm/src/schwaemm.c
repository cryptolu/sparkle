///////////////////////////////////////////////////////////////////////////////
// schwaemm.c: Optimized C implementation of the AEAD algorithm SCHWAEMM.    //
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

// gencat_aead.c shall be used to generate the test vector output file. The
// test vector output file shall be provided in the corresponding
// crypto_aead/[algorithm]/ directory


#include <string.h>  // for memcpy, memset
#include "schwaemm.h"
#include "sparkle.h"


///////////////////////////////////////////////////////////////////////////////
////////// SPARKLE CONFIGURATIONS FOR THE FOUR INSTANCES OF SCHWAEMM //////////
///////////////////////////////////////////////////////////////////////////////


#if (SCHWAEMM_INST == 128128)
#define STATE_BYTES     32
#define RATE_BYTES      16
#define STEPS_SLIM       7
#define STEPS_BIG       10

#elif (SCHWAEMM_INST == 256128)
#define STATE_BYTES     48
#define RATE_BYTES      32
#define STEPS_SLIM       7
#define STEPS_BIG       11

#elif (SCHWAEMM_INST == 192192)
#define STATE_BYTES     48
#define RATE_BYTES      24
#define STEPS_SLIM       7
#define STEPS_BIG       11

#elif (SCHWAEMM_INST == 256256)
#define STATE_BYTES     64
#define RATE_BYTES      32
#define STEPS_SLIM       8
#define STEPS_BIG       12

#else
#error "Invalid definition of SCHWAEMM instance!"
#endif

#define SCHWAEMM_KEY_WORDS   (SCHWAEMM_KEY_BYTES/4)
#define SCHWAEMM_NONCE_WORDS (SCHWAEMM_NONCE_BYTES/4)
#define SCHWAEMM_TAG_WORDS   (SCHWAEMM_TAG_BYTES/4)

#define STATE_BRANS (STATE_BYTES/8)
#define STATE_WORDS (STATE_BYTES/4)
#define RATE_BRANS  (RATE_BYTES/8)
#define RATE_WORDS  (RATE_BYTES/4)
#define CAP_BYTES   (STATE_BYTES-RATE_BYTES)
#define CAP_BRANS   (CAP_BYTES/8)
#define CAP_WORDS   (CAP_BYTES/4)

#define CONST_A0 (((uint32_t) (0 ^ (1 << CAP_BRANS))) << 24)
#define CONST_A1 (((uint32_t) (1 ^ (1 << CAP_BRANS))) << 24)
#define CONST_M2 (((uint32_t) (2 ^ (1 << CAP_BRANS))) << 24)
#define CONST_M3 (((uint32_t) (3 ^ (1 << CAP_BRANS))) << 24)


///////////////////////////////////////////////////////////////////////////////
//// PREPROCESSOR DIRECTIVES TO REPLACE THE C CODE OF SPARKLE BY ASM CODE /////
///////////////////////////////////////////////////////////////////////////////


// When this file is compiled for an AVR microcontroller and SPARKLE_ASSEMBLER
// is defined (see schwaemm.h), then the AVR assembler implementation of the
// SPARKLE permutation is used. On the other hand, if SPARKLE_ASSEMBLER is not
// defined, then the C version (i.e. the function sparkle) is used.

#if (defined(__AVR) || defined(__AVR__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_avr(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_avr((state), (brans), (steps))
#endif // if (defined(__AVR__) || ...


// When this file is compiled for a MSP430 (or a MSP430X) microcontroller and
// SPARKLE_ASSEMBLER is defined (see schwaemm.h), then the MSP430 assembler
// implementation of the SPARKLE permutation is used. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// sparkle) is used.

#if (defined(MSP430) || defined(__MSP430__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_msp(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_msp((state), (brans), (steps))
#endif // if (defined(MSP430) || ...


// When this file is compiled for an ARM microcontroller and SPARKLE_ASSEMBLER
// is defined (see schwaemm.h), then one of the three branch-unrolled ARM
// assembler implementations of the SPARKLE permutation is used, depending on
// the concrete SCHWAEMM instance. On the other hand, if SPARKLE_ASSEMBLER is
// not defined, then the C version (i.e. the function sparkle) is used.

#if (defined(__arm__) || defined(_M_ARM)) && defined(SPARKLE_ASSEMBLER)
#if (STATE_BYTES == 32)
extern void sparkle256_arm(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle256_arm((state), (steps))
#elif (STATE_BYTES == 48)
extern void sparkle384_arm(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle384_arm((state), (steps))
#elif (STATE_BYTES == 64)
extern void sparkle512_arm(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle512_arm((state), (steps))
#endif // if (STATE_BYTES == 32)
#endif // if (defined(__arm__) || ...


// When this file is compiled for a 32-bit RISC-V microcontroller (e.g. RV32I)
// and SPARKLE_ASSEMBLER is defined (see schwaemm.h), then one of the three
// branch-unrolled RV32 assembler implementations of the SPARKLE permutation is
// used, depending on the concrete SCHWAEMM instance. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// sparkle) is used.

#if defined(__riscv_xlen) && (__riscv_xlen == 32) && defined(SPARKLE_ASSEMBLER)
#if (STATE_BYTES == 32)
extern void sparkle256_rv32(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle256_rv32((state), (steps))
#elif (STATE_BYTES == 48)
extern void sparkle384_rv32(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle384_rv32((state), (steps))
#elif (STATE_BYTES == 64)
extern void sparkle512_rv32(uint32_t *state, int steps);
#define sparkle(state, brans, steps) sparkle512_rv32((state), (steps))
#endif // if (STATE_BYTES == 32)
#endif // if defined(__riscv_xlen) && ...


///////////////////////////////////////////////////////////////////////////////
/////// HELPER FUNCTIONS AND MACROS (RHO1, RHO2, RATE-WHITENING, ETC.) ////////
///////////////////////////////////////////////////////////////////////////////


// The high-level AEAD API specifies that the associated data, plaintext, and
// ciphertext are stored in arrays of type unsigned char. However, the SPARKLE
// permutation operates on 32-bit words and performs best when the data to be
// processed is stored in an uint32_t-array. Casting an unsigned-char pointer
// to an uint32_t-pointer increases the alignment requirements, i.e. the base
// address of the array must be even on 16-bit architectures and a multiple of
// four (i.e. 4-byte aligned) on 32-bit and 64-bit platforms. The preprocessor
// statements below can be used to determine the alignment requirements an
// unsigned-char-pointer has to meet to permit casting to an uint32_t-pointer.

#if defined(_MSC_VER) && (_MSC_VER < 1600)
#define ALIGN_OF_UI32 4
#else  // compiler is not ancient MSVC
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)  // C11
#include <stdalign.h>
#define ALIGN_OF_UI32 alignof(uint32_t)
#else  // C11 standard is not available
#define MIN_SIZE(a, b) ((sizeof(a) < sizeof(b)) ? sizeof(a) : sizeof(b))
#define ALIGN_OF_UI32 MIN_SIZE(uint32_t, uint_fast8_t)  // stdint.h
#endif // if defined(__STDC_VERSION__) && ...
#endif // defined(_MSC_VER) && ...


// The rate-whitening for SCHWAEMM256_128 applies the "tweak" described in
// Section 2.3.2 of the specification. Therefore, the indices used to load the
// 32-bit words from the capacity-part of the state need to be reduced modulo
// CAP_WORDS, which the C implementation below does by ANDing the index with
// (CAP_WORDS - 1) = 3. Performing the modulo reduction in this way only works
// when CAP_WORDS is a power of 2, which is the case for SCHWAEMM256_128.

#if (RATE_WORDS > CAP_WORDS)
#define CAP_INDEX(i) ((i) & (CAP_WORDS-1))
#else  // RATE_WORDS <= CAP_WORDS
#define CAP_INDEX(i) (i)
#endif


// Rho and rate-whitening for the authentication of associated data. The third
// parameter indicates whether the uint8_t-pointer `in` is properly aligned to
// permit casting to a uint32_t-pointer. If this is the case then array `in` is
// processed directly, otherwise it is first copied to an aligned buffer.

static void rho_whi_aut(uint32_t *state, const uint8_t *in, int aligned)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t *in32;
  uint32_t tmp;
  int i, j;
  
  if (aligned) {  // `in` can be casted to uint32_t pointer
    in32 = (void *) in;  // casting in this way prevents a warning
  } else {  // `in` is not sufficiently aligned for casting
    memcpy(buffer, in, RATE_BYTES);
    in32 = (uint32_t *) buffer;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp = state[i];
    state[i] = state[j] ^ in32[i] ^ state[RATE_WORDS+i];
    state[j] ^= tmp ^ in32[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
  }
}


// Rho and rate-whitening for the authentication of the last associated-data
// block. Since this last block may require padding, it is always copied to a
// buffer.

static void rho_whi_aut_last(uint32_t *state, const uint8_t *in, size_t inlen)
{
  uint32_t buffer[RATE_WORDS];
  uint8_t *bufptr;
  uint32_t tmp;
  int i, j;
  
  memcpy(buffer, in, inlen);
  if (inlen < RATE_BYTES) {  // padding
    bufptr = ((uint8_t *) buffer) + inlen;
    memset(bufptr, 0, (RATE_BYTES - inlen));
    *bufptr = 0x80;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp = state[i];
    state[i] = state[j] ^ buffer[i] ^ state[RATE_WORDS+i];
    state[j] ^= tmp ^ buffer[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
  }
}


// Rho and rate-whitening for the encryption of plaintext. The third parameter
// indicates whether the uint8_t-pointers `in` and `out` are properly aligned
// to permit casting to uint32_t-pointers. If this is the case then array `in`
// and `out` are processed directly, otherwise `in` is copied to an aligned
// buffer.

static void rho_whi_enc(uint32_t *state, uint8_t *out, const uint8_t *in, \
  int aligned)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t *in32, *out32;
  uint32_t tmp1, tmp2;
  int i, j;
  
  if (aligned) {  // `in` and `out` can be casted to uint32_t pointer
    in32 = (void *) in;  // casting in this way prevents a warning
    out32 = (void *) out;  // casting in this way prevents a warning
  } else {  // `in` or `out` is not sufficiently aligned for casting
    memcpy(buffer, in, RATE_BYTES);
    in32 = out32 = (uint32_t *) buffer;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp1 = state[i];
    tmp2 = state[j];
    state[i] = state[j] ^ in32[i] ^ state[RATE_WORDS+i];
    state[j] ^= tmp1 ^ in32[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
    out32[i] = in32[i] ^ tmp1;
    out32[j] = in32[j] ^ tmp2;
  }
  
  if (!aligned) {
    memcpy(out, buffer, RATE_BYTES);
  }
}


// Rho and rate-whitening for the encryption of the last plaintext block. Since
// this last block may require padding, it is always copied to a buffer.

static void rho_whi_enc_last(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t tmp1, tmp2;
  uint8_t *bufptr;
  int i, j;
  
  memcpy(buffer, in, inlen);
  if (inlen < RATE_BYTES) {  // padding
    bufptr = ((uint8_t *) buffer) + inlen;
    memset(bufptr, 0, (RATE_BYTES - inlen));
    *bufptr = 0x80;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp1 = state[i];
    tmp2 = state[j];
    state[i] = state[j] ^ buffer[i] ^ state[RATE_WORDS+i];
    state[j] ^= tmp1 ^ buffer[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
    buffer[i] ^= tmp1;
    buffer[j] ^= tmp2;
  }
  memcpy(out, buffer, inlen);
}


// Rho and rate-whitening for the decryption of ciphertext. The third parameter
// indicates whether the uint8_t-pointers `in` and `out` are properly aligned
// to permit casting to uint32_t-pointers. If this is the case then array `in`
// and `out` are processed directly, otherwise `in` is copied to an aligned
// buffer.

static void rho_whi_dec(uint32_t *state, uint8_t *out, const uint8_t *in, \
  int aligned)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t *in32, *out32;
  uint32_t tmp1, tmp2;
  int i, j;
  
  if (aligned) {  // `in` and `out` can be casted to uint32_t pointer
    in32 = (void *) in;  // casting in this way prevents a warning
    out32 = (void *) out;  // casting in this way prevents a warning
  } else {  // `in` or `out` is not sufficiently aligned for casting
    memcpy(buffer, in, RATE_BYTES);
    in32 = out32 = (uint32_t *) buffer;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp1 = state[i];
    tmp2 = state[j];
    state[i] ^= state[j] ^ in32[i] ^ state[RATE_WORDS+i];
    state[j] = tmp1 ^ in32[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
    out32[i] = in32[i] ^ tmp1;
    out32[j] = in32[j] ^ tmp2;
  }
  
  if (!aligned) {
    memcpy(out, buffer, RATE_BYTES);
  }
}


// Rho and rate-whitening for the decryption of the last ciphertext block.
// Since this last block may require padding, it is always copied to a buffer.

static void rho_whi_dec_last(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen)
{
  uint32_t buffer[RATE_WORDS];
  uint32_t tmp1, tmp2;
  uint8_t *bufptr;
  int i, j;
  
  memcpy(buffer, in, inlen);
  if (inlen < RATE_BYTES) {  // padding
    bufptr = ((uint8_t *) buffer) + inlen;
    memcpy(bufptr, (((uint8_t *) state) + inlen), (RATE_BYTES - inlen));
    *bufptr ^= 0x80;
  }
  
  for (i = 0, j = RATE_WORDS/2; i < RATE_WORDS/2; i++, j++) {
    tmp1 = state[i];
    tmp2 = state[j];
    state[i] ^= state[j] ^ buffer[i] ^ state[RATE_WORDS+i];
    state[j] = tmp1 ^ buffer[j] ^ state[RATE_WORDS+CAP_INDEX(j)];
    buffer[i] ^= tmp1;
    buffer[j] ^= tmp2;
  }
  memcpy(out, buffer, inlen);
}


///////////////////////////////////////////////////////////////////////////////
///////////// LOW-LEVEL AEAD FUNCTIONS (FOR USE WITH FELICS-AEAD) /////////////
///////////////////////////////////////////////////////////////////////////////


// The Initialize function loads nonce and key into the state and executes the
// SPARKLE permutation with the big number of steps.

void Initialize(uint32_t *state, const uint8_t *key, const uint8_t *nonce)
{
  // load nonce into the rate-part of the state
  memcpy(state, nonce, SCHWAEMM_NONCE_BYTES);
   // load key into the capacity-part of the sate
  memcpy((state + RATE_WORDS), key, SCHWAEMM_KEY_BYTES);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


// The ProcessAssocData function absorbs the associated data, which becomes
// only authenticated but not encrypted, into the state (in blocks of size
// RATE_BYTES). Note that this function MUST NOT be called when the length of
// the associated data is 0.

void ProcessAssocData(uint32_t *state, const uint8_t *in, size_t inlen)
{
  // check whether `in` can be casted to uint32_t pointer (we can use here
  // size_t instead of uintptr_t since ALIGN_OF_UI32 is either 1, 2, or 4)
  int aligned = ((size_t) in) % ALIGN_OF_UI32 == 0;
  // printf("Address of `in`: %p\n", in);
  
  // Main Authentication Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_aut(state, in, aligned);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Authentication of Last Block
  
  // addition of constant A0 or A1 to the state
  state[STATE_WORDS-1] ^= ((inlen < RATE_BYTES) ? CONST_A0 : CONST_A1);
  // combined Rho and rate-whitening (incl. padding)
  rho_whi_aut_last(state, in, inlen);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


// The ProcessPlainText function encrypts the plaintext (in blocks of size
// RATE_BYTES) and generates the respective ciphertext. The uint8_t-array `in`
// contains the plaintext and the ciphertext is written to uint8_t-array `out`
// (`in` and `out` can be the same array, i.e. they can have the same start
// address). Note that this function MUST NOT be called when the length of the
// plaintext is 0.

void ProcessPlainText(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen)
{
  // check whether `in` and `out` can be casted to uint32_t pointer (we can use
  // here size_t instead of uintptr_t since ALIGN_OF_UI32 is either 1, 2, or 4)
  int aligned = (((size_t) in) | ((size_t) out)) % ALIGN_OF_UI32 == 0;
  // printf("Address of `in` and `out`: %p, %p\n", in, out);
  
  // Main Encryption Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_enc(state, out, in, aligned);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    out += RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Encryption of Last Block
  
  // addition of constant M2 or M3 to the state
  state[STATE_WORDS-1] ^= ((inlen < RATE_BYTES) ? CONST_M2 : CONST_M3);
  // combined Rho and rate-whitening (incl. padding)
  rho_whi_enc_last(state, out, in, inlen);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


// The Finalize function adds the key to the capacity part of the state.

void Finalize(uint32_t *state, const uint8_t *key)
{
  uint32_t buffer[SCHWAEMM_KEY_WORDS];
  int i;
  
  // to prevent (potentially) unaligned memory accesses
  memcpy(buffer, key, SCHWAEMM_KEY_BYTES);
  // add key to the capacity-part of the state
  for (i = 0; i < SCHWAEMM_KEY_WORDS; i++) {
    state[RATE_WORDS+i] ^= buffer[i];
  }
}


// The GenerateTag function generates an authentication tag.

void GenerateTag(uint32_t *state, uint8_t *tag)
{
  memcpy(tag, (state + RATE_WORDS), SCHWAEMM_TAG_BYTES);
}


// The VerifyTag function checks whether the given authentication tag is valid.
// It performs a simple constant-time comparison and returns 0 if the provided
// tag matches the computed tag and -1 otherwise.

int VerifyTag(uint32_t *state, const uint8_t *tag)
{
  uint32_t buffer[SCHWAEMM_TAG_WORDS], diff = 0;
  int i;
  
  // to prevent (potentially) unaligned memory accesses
  memcpy(buffer, tag, SCHWAEMM_TAG_BYTES);
  // constant-time comparison: 0 if equal, -1 otherwise
  for (i = 0; i < SCHWAEMM_TAG_WORDS; i++) {
    diff |= (state[RATE_WORDS+i] ^ buffer[i]);
  }

  return (((int) (diff == 0)) - 1);
}


// The ProcessCipherText function decrypts the ciphertext (in blocks of size
// RATE_BYTES) and generates the respective plaintext. The uint8_t-array `in`
// contains the ciphertext and the plaintext is written to uint8_t-array `out`
// (`in` and `out` can be the same array, i.e. they can have the same start
// address). Note that this function MUST NOT be called when the length of the
// ciphertext is 0.

void ProcessCipherText(uint32_t *state, uint8_t *out, const uint8_t *in, \
  size_t inlen)
{
  // check whether `in` and `out` can be casted to uint32_t pointer (we can use
  // here size_t instead of uintptr_t since ALIGN_OF_UI32 is either 1, 2, or 4)
  int aligned = (((size_t) in) | ((size_t) out)) % ALIGN_OF_UI32 == 0;
  // printf("Address of `in` and `out`: %p, %p\n", in, out);
  
  // Main Decryption Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_dec(state, out, in, aligned);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    out += RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Decryption of Last Block
  
  // addition of constant M2 or M3 to the state
  state[STATE_WORDS-1] ^= ((inlen < RATE_BYTES) ? CONST_M2 : CONST_M3);
  // combined Rho and rate-whitening (incl. padding)
  rho_whi_dec_last(state, out, in, inlen);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


///////////////////////////////////////////////////////////////////////////////
////////////// HIGH-LEVEL AEAD FUNCTIONS (FOR USE WITH SUPERCOP) //////////////
///////////////////////////////////////////////////////////////////////////////


// High-level encryption function from SUPERCOP.
// nsec is kept for compatibility with SUPERCOP, but is not used.

int crypto_aead_encrypt(UChar *c, ULLInt *clen, const UChar *m, ULLInt mlen, \
  const UChar *ad, ULLInt adlen, const UChar *nsec, const UChar *npub,       \
  const UChar *k)
{
  uint32_t state[STATE_WORDS];
  size_t msize = (size_t) mlen;
  size_t adsize = (size_t) adlen;
  
  (void) nsec;  // to get rid of a warning
  
  Initialize(state, k, npub);
  if (adsize) {
    ProcessAssocData(state, ad, adsize);
  }
  if (msize) {
    ProcessPlainText(state, c, m, msize);
  }
  Finalize(state, k);
  GenerateTag(state, (c + msize));

  *clen = msize;
  *clen += SCHWAEMM_TAG_BYTES;
  
  return 0;
}


// High-level decryption function from SUPERCOP.
// nsec is kept for compatibility with SUPERCOP, but is not used.

int crypto_aead_decrypt(UChar *m, ULLInt *mlen, UChar *nsec, const UChar *c, \
  ULLInt clen, const UChar *ad, ULLInt adlen, const UChar *npub,             \
  const UChar *k)
{
  uint32_t state[STATE_WORDS];
  size_t csize = (size_t) (clen - SCHWAEMM_TAG_BYTES);
  size_t adsize = (size_t) adlen;
  int retval;
  
  (void) nsec;  // to get rid of a warning
  
  Initialize(state, k, npub);
  if (adsize) {
    ProcessAssocData(state, ad, adsize);
  }
  if (csize) {
    ProcessCipherText(state, m, c, csize);
  }
  Finalize(state, k);
  retval = VerifyTag(state, (c + csize));
  *mlen = csize;
  
  return retval;
}
