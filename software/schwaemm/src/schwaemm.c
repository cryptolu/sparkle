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
// defined, then the C version (i.e. the function `sparkle`) is used.

#if (defined(__AVR) || defined(__AVR__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_avr(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_avr((state), (brans), (steps))
#endif // if (defined(__AVR__) || ...


// When this file is compiled for a MSP430 (or a MSP430X) microcontroller and
// SPARKLE_ASSEMBLER is defined (see schwaemm.h), then the MSP430 assembler
// implementation of the SPARKLE permutation is used. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// `sparkle`) is used.

#if (defined(__MSP430__) || defined(__ICC430__)) && defined(SPARKLE_ASSEMBLER)
extern void sparkle_msp(uint32_t *state, int brans, int steps);
#define sparkle(state, brans, steps) sparkle_msp((state), (brans), (steps))
#endif // if (defined(__MSP430__) || ...


// When this file is compiled for a 32-bit RISC-V microcontroller (e.g. RV32I)
// and SPARKLE_ASSEMBLER is defined (see schwaemm.h), then one of the three
// branch-unrolled RV32 assembler implementations of the SPARKLE permutation is
// used, depending on the concrete SCHWAEMM instance. On the other hand, if
// SPARKLE_ASSEMBLER is not defined, then the C version (i.e. the function
// `sparkle`) is used.

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


// When this file is compiled for an ARM microcontroller and SPARKLE_ASSEMBLER
// is defined (see schwaemm.h), then one of the three branch-unrolled ARM
// assembler implementations of the SPARKLE permutation is used, depending on
// the concrete SCHWAEMM instance. On the other hand, if SPARKLE_ASSEMBLER is
// not defined, then the C version (i.e. the function `sparkle`) is used.

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


///////////////////////////////////////////////////////////////////////////////
/////// HELPER FUNCTIONS AND MACROS (RHO1, RHO2, RATE-WHITENING, ETC.) ////////
///////////////////////////////////////////////////////////////////////////////


// The high-level AEAD API specifies that the associated data, plaintext, and
// ciphertext are stored in arrays of type unsigned char. However, the SPARKLE
// permutation operates on 32-bit words and performs best when the data to be
// processed is stored in an uint32_t-array. Casting an unsigned-char-pointer
// to a uint32_t-pointer is only permitted if the pointer is properly aligned,
// which means the address must be even on 16-bit architectures and a multiple
// of four (i.e. 4-byte aligned) on 32-bit and 64-bit platforms, though there
// exist also architectures that can handle unaligned RAM accesses, e.g. Intel
// X86/X64 and ARM Cortex-M3/M4. The following preprocessor statements try to
// determine the alignment requirements that an unsigned-char-pointer needs to
// satisfy to allow casting to a uint32_t-pointer. Unfortunately, there is no
// clean way to do this with pre-C11 compilers, but to be on the safe side the
// required alignment can be assumed to be the processor's word-size or to the
// size of uint32_t, whichever value is smaller. While this approach reliably
// prevents unaligned RAM accesses, it may be overly conservative and increase
// the execution time, in particular on architectures that do not have strict
// alignment requirements.

#ifndef ALIGN_OF_UI32
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)  // C11
#include <stdalign.h>
#define ALIGN_OF_UI32 alignof(uint32_t)
#else  // C11 standard is not supported
#define MIN_SIZE(a, b) ((sizeof(a) < sizeof(b)) ? sizeof(a) : sizeof(b))
#define ALIGN_OF_UI32 MIN_SIZE(uint32_t, uint_fast8_t)  // stdint.h
#endif // if defined(__STDC_VERSION__) && ...
#endif // ifndef ALIGN_OF_UI32


// Rho and rate-whitening for the authentication of associated data. The third
// parameter `bytes2buf` specifies the number of bytes of array `in` that have
// to be copied to an aligned buffer; a value of 0 indicates that the uint8_t-
// pointer `in` is properly aligned to permit casting to a uint32_t-pointer,
// in which case the bytes can be directly processed (without copying them to
// an aligned buffer). This implementation supports any SCHWAEMM instance with
// RATE_WORDS == CAP_WORDS or RATE_WORDS == 2*CAP_WORDS.

static void rho_whi_aut(uint32_t *state, const uint8_t *in, int bytes2buf)
{
  uint32_t buffer[RATE_WORDS+1];
  uint32_t rw1, rw2;  // left/right rate-word
  uint32_t dw1, dw2;  // left/right data-word
  uint32_t *in32;
  int i;
  
  if (bytes2buf == 0) {  // `in` can be directly processed
    in32 = (uint32_t *) (void *) in;  // to prevent cast-warning
  } else {  // bytes of array `in` are copied to (padded) buffer
    i = (int) (bytes2buf >> 2);  // start-index of the 0-words
    memset((buffer + i), 0, 4*(RATE_WORDS - i));
    ((uint8_t *) buffer)[bytes2buf] = 0x80;
    memcpy(buffer, in, bytes2buf);
    in32 = (uint32_t *) buffer;
  }
  
  for (i = 0; i < RATE_WORDS/2; i++) {
    rw1 = state[i];
    rw2 = state[RATE_WORDS/2+i];
    dw1 = in32[i];
    dw2 = in32[RATE_WORDS/2+i];
    dw2 ^= rw2;  // part of Rho1
    rw1 ^= dw2;  // Feistel-XOR
    rw2 ^= dw1;  // part of Rho1
    // Rate-Whitening + Feistel-Swap rw1 <-> rw2
    rw2 ^= state[RATE_WORDS+i];
    rw1 ^= state[(STATE_WORDS-RATE_WORDS/2)+i];
    state[i] = rw2;
    state[RATE_WORDS/2+i] = rw1; 
  }
}


// Rho and rate-whitening for the encryption of plaintext. The third parameter
// `bytes2buf` specifies the number of bytes of array `in` that have to be
// copied to an aligned buffer; a value of 0 indicates that the uint8_t-pointer
// `in` is properly aligned to permit casting to a uint32_t-pointer, in which
// case the bytes can be directly processed (without copying them to an aligned
// buffer). This implementation supports any SCHWAEMM instance with RATE_WORDS
// == CAP_WORDS or RATE_WORDS == 2*CAP_WORDS.

static void rho_whi_enc(uint32_t *state, uint8_t *out, const uint8_t *in, \
  int bytes2buf)
{
  uint32_t buffer[RATE_WORDS+1];
  uint32_t rw1, rw2;  // left/right rate-word
  uint32_t dw1, dw2;  // left/right ptxt-word
  uint32_t *in32, *out32;
  int i;
  
  if (bytes2buf == 0) {   // `in`, `out` can be directly processed
    in32 = (uint32_t *) (void *) in;    // to prevent cast-warning
    out32 = (uint32_t *) (void *) out;  // to prevent cast-warning
  } else {  // bytes of array `in` are copied to a (padded) buffer
    i = (int) (bytes2buf >> 2);  // start-index of the 0-words
    memset((buffer + i), 0, 4*(RATE_WORDS - i));
    ((uint8_t *) buffer)[bytes2buf] = 0x80;
    memcpy(buffer, in, bytes2buf);
    in32 = out32 = (uint32_t *) buffer;
  }
  
  for (i = 0; i < RATE_WORDS/2; i++) {
    rw1 = state[i];
    rw2 = state[RATE_WORDS/2+i];
    dw1 = in32[i];
    dw2 = in32[RATE_WORDS/2+i];
    dw2 ^= rw2;  // part of Rho2
    rw2 ^= dw1;  // part of Rho1
    dw1 ^= rw1;  // part of Rho2
    rw1 ^= dw2;  // part of Rho1 + Feistel-XOR
    out32[i] = dw1;
    out32[RATE_WORDS/2+i] = dw2;
    // Rate-Whitening + Feistel-Swap rw1 <-> rw2
    rw2 ^= state[RATE_WORDS+i];
    rw1 ^= state[(STATE_WORDS-RATE_WORDS/2)+i];
    state[i] = rw2;
    state[RATE_WORDS/2+i] = rw1;
  }
  
  if (bytes2buf) memcpy(out, buffer, bytes2buf);
}


// Rho and rate-whitening for the decryption of ciphertext. The third parameter
// `bytes2buf` specifies the number of bytes of array `in` that have to be
// copied to an aligned buffer; a value of 0 indicates that the uint8_t-pointer
// `in` is properly aligned to permit casting to a uint32_t-pointer, in which
// case the bytes can be directly processed (without copying them to an aligned
// buffer). This implementation supports any SCHWAEMM instance with RATE_WORDS
// == CAP_WORDS or RATE_WORDS == 2*CAP_WORDS. Instead of padding the last block
// with 0-bytes, we pad it with the corresponding state-bytes since, in this
// way, the rho'1 function for decryption (see specification Section 2.3.2) can
// omit the XOR of the state. Thanks to this small modification of the padding,
// the rho and rate-whitening for decryption is equally fast as for encryption.

static void rho_whi_dec(uint32_t *state, uint8_t *out, const uint8_t *in, \
  int bytes2buf)
{
  uint32_t buffer[RATE_WORDS+1];
  uint32_t rw1, rw2;  // left/right rate-word
  uint32_t dw1, dw2;  // left/right ctxt-word
  uint32_t *in32, *out32;
  int i;
  
  if (bytes2buf == 0) {   // `in`, `out` can be directly processed
    in32 = (uint32_t *) (void *) in;    // to prevent cast-warning
    out32 = (uint32_t *) (void *) out;  // to prevent cast-warning
  } else {  // bytes of array `in` are copied to a (padded) buffer
    i = (int) (bytes2buf >> 2);  // start-index of the state-words
    memcpy((buffer + i), (state + i), 4*(RATE_WORDS - i));
    ((uint8_t *) buffer)[bytes2buf] ^= 0x80;
    memcpy(buffer, in, bytes2buf);
    in32 = out32 = (uint32_t *) buffer;
  }
  
  for (i = 0; i < RATE_WORDS/2; i++) {
    rw1 = state[i];
    rw2 = state[RATE_WORDS/2+i];
    dw1 = in32[i];
    dw2 = in32[RATE_WORDS/2+i];
    dw1 ^= rw1;  // part of Rho2
    rw1 ^= dw2;  // part of Rho1 + Feistel-XOR
    dw2 ^= rw2;  // part of Rho2
    rw2 ^= dw1;  // part of Rho1
    out32[i] = dw1;
    out32[RATE_WORDS/2+i] = dw2;
    // Rate-Whitening + Feistel-Swap rw1 <-> rw2
    rw2 ^= state[RATE_WORDS+i];
    rw1 ^= state[(STATE_WORDS-RATE_WORDS/2)+i];
    state[i] = rw2;
    state[RATE_WORDS/2+i] = rw1;
  }
  
  if (bytes2buf) memcpy(out, buffer, bytes2buf);
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
  // check whether `in` can be casted to a uint32_t pointer
  int aligned = ((uintptr_t) in) % (ALIGN_OF_UI32) == 0;
  // `bytes2buf` is the number of bytes that `rho_whi_aut` has to copy to an
  // aligned buffer; it is 0 when `in` is aligned to permit direct processing
  // (w/o copying to a buffer) by casting it to a uint32_t pointer
  int bytes2buf = (aligned) ? 0 : RATE_BYTES;
  
  // Main Authentication Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_aut(state, in, bytes2buf);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Authentication of Last Block
  
  if (inlen == RATE_BYTES) {
    state[STATE_WORDS-1] ^= CONST_A1;
  } else {  // inlen < RATE_BYTES
    state[STATE_WORDS-1] ^= CONST_A0;
    bytes2buf = (int) inlen;
  }
  // combined Rho and rate-whitening operation
  rho_whi_aut(state, in, bytes2buf);
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
  // check whether `in` and `out` can be casted to uint32_t pointers
  int aligned = (((uintptr_t) in) | ((uintptr_t) out)) % (ALIGN_OF_UI32) == 0;
  // `bytes2buf` is the number of bytes that `rho_whi_enc` has to copy to an
  // aligned buffer; it is 0 when `in` and `out` are aligned to permit direct
  // processing (w/o copying to a buffer) by casting them to uint32_t pointers
  int bytes2buf = (aligned) ? 0 : RATE_BYTES;
  
  // Main Encryption Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_enc(state, out, in, bytes2buf);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    out += RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Encryption of Last Block
  
  if (inlen == RATE_BYTES) {
    state[STATE_WORDS-1] ^= CONST_M3;
  } else {  // inlen < RATE_BYTES
    state[STATE_WORDS-1] ^= CONST_M2;
    bytes2buf = (int) inlen;
  }
  // combined Rho and rate-whitening operation
  rho_whi_enc(state, out, in, bytes2buf);
  // execute SPARKLE with big number of steps
  sparkle(state, STATE_BRANS, STEPS_BIG);
}


// The Finalize function adds the key to the capacity part of the state.

void Finalize(uint32_t *state, const uint8_t *key)
{
  uint32_t buffer[SCHWAEMM_KEY_WORDS];
  uint32_t *key32;  // pointer for 32-bit access to `key`
  // check whether `key` can be casted to uint32_t pointer
  int i, aligned = ((uintptr_t) key) % (ALIGN_OF_UI32) == 0;
  
  if (aligned) {  // `key` can be casted to uint32_t pointer
    key32 = (uint32_t *) (void *) key;  // to prevent cast-warning
  } else {  // `key` is not sufficiently aligned for casting
    memcpy(buffer, key, SCHWAEMM_KEY_BYTES);
    key32 = (uint32_t *) buffer;
  }
  
  // add key to the capacity-part of the state
  for (i = 0; i < SCHWAEMM_KEY_WORDS; i++) {
    state[RATE_WORDS+i] ^= key32[i];
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
  uint32_t *tag32;  // pointer for 32-bit access to `tag`
  // check whether `tag` can be casted to uint32_t pointer
  int i, aligned = ((uintptr_t) tag) % (ALIGN_OF_UI32) == 0;
  
  if (aligned) {  // `tag` can be casted to uint32_t pointer
    tag32 = (uint32_t *) (void *) tag;  // to prevent cast-warning
  } else {  // `tag` is not sufficiently aligned for casting
    memcpy(buffer, tag, SCHWAEMM_TAG_BYTES);
    tag32 = (uint32_t *) buffer;
  }
  
  // constant-time comparison: 0 if equal, -1 otherwise
  for (i = 0; i < SCHWAEMM_TAG_WORDS; i++) {
    diff |= (state[RATE_WORDS+i] ^ tag32[i]);
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
  // check whether `in` and `out` can be casted to uint32_t pointers
  int aligned = (((uintptr_t) in) | ((uintptr_t) out)) % (ALIGN_OF_UI32) == 0;
  // `bytes2buf` is the number of bytes that `rho_whi_dec` has to copy to an
  // aligned buffer; it is 0 when `in` and `out` are aligned to permit direct
  // processing (w/o copying to a buffer) by casting them to uint32_t pointers
  int bytes2buf = (aligned) ? 0 : RATE_BYTES;
  
  // Main Decryption Loop
  
  while (inlen > RATE_BYTES) {
    // combined Rho and rate-whitening operation
    rho_whi_dec(state, out, in, bytes2buf);
    // execute SPARKLE with slim number of steps
    sparkle(state, STATE_BRANS, STEPS_SLIM);
    inlen -= RATE_BYTES;
    out += RATE_BYTES;
    in += RATE_BYTES;
  }
  
  // Decryption of Last Block
  
  if (inlen == RATE_BYTES) {
    state[STATE_WORDS-1] ^= CONST_M3;
  } else {  // inlen < RATE_BYTES
    state[STATE_WORDS-1] ^= CONST_M2;
    bytes2buf = (int) inlen;
  }
  // combined Rho and rate-whitening operation
  rho_whi_dec(state, out, in, bytes2buf);
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
  if (adsize > 0) ProcessAssocData(state, ad, adsize);
  if (msize > 0) ProcessPlainText(state, c, m, msize);
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
  if (clen < SCHWAEMM_TAG_BYTES) return -1;
  
  Initialize(state, k, npub);
  if (adsize > 0) ProcessAssocData(state, ad, adsize);
  if (csize > 0) ProcessCipherText(state, m, c, csize);
  Finalize(state, k);
  retval = VerifyTag(state, (c + csize));
  *mlen = csize;
  
  return retval;
}
