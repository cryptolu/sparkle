///////////////////////////////////////////////////////////////////////////////
// sparkle384_32i.S: RV32I assembler implementation (GCC) of SPARKLE384.     //
// This file is part of the SPARKLE package that was sent to NIST during the //
// 3rd round of the Lightweight Cryptography (LWC) standardization project.  //
// Version 1.1.2 (18-Oct-21), see <http://github.com/cryptolu/> for updates. //
// Authors: the SPARKLE group (Christof Beierle, Alex Biryukov, Luan Cardoso //
// dos Santos, Johann Groszschaedl, Amir Moradi, Leo Perrin, Aein Rezaei     //
// Shahmirzadi, Aleksei Udovenko, Vesselin Velichkov, and Qingju Wang).      //
// License: LGPLv3 (see LICENSE.md), other licenses available upon request.  //
// Copyright (C) 2019-2021 University of Luxembourg <http://www.uni.lu/>.    //
// ------------------------------------------------------------------------- //
// This source code is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU Lesser General Public License as published  //
// by the Free Software Foundation, either version 3 of the License, or (at  //
// your option) any later version. This source code is distributed in the    //
// hope that it will be useful, but WITHOUT ANY WARRANTY; without even the   //
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  //
// See the GNU Lesser General Public License for more details. You should    //
// have received a copy of the GNU Lesser General Public License along with  //
// this source code. If not, see <http://www.gnu.org/licenses/>.             //
///////////////////////////////////////////////////////////////////////////////


// Function prototype:
// -------------------
// void sparkle384_rv32(uint32_t *state, int steps)
//
// Parameters:
// -----------
// `state`: pointer to an uint32_t-array containing the 12 state words
// `steps`: number of steps (must be either 7 or 11)
//
// Return value:
// -------------
// None


// This file contains a size-optimized and a speed-optimized implementation of
// SPARKLE384. By default, the size-optimized implementation is used. In order
// to switch to the speed-optimized version, the symbol `SPARKLE_ASM_FAST` has
// to be defined, either by uncommenting the line below or by compiling this
// file with the option `-DSPARKLE_ASM_FAST`.

// #define SPARKLE_ASM_FAST


.section .text
.balign 4


///////////////////////////////////////////////////////////////////////////////
//////////////////////// REGISTER NAMES AND CONSTANTS /////////////////////////
///////////////////////////////////////////////////////////////////////////////

// register `sptr` holds the start address of array `state`
#define sptr a0
// register `steps` holds the number of steps (parameter `steps`)
#define steps a1
// registers `x0w` to `y5w` hold 12 words from array `state`
#define x0w a2
#define y0w a3
#define x1w a4
#define y1w a5
#define x2w a6
#define y2w a7
#define x3w t0
#define y3w t1
#define x4w t2
#define y4w t3
#define x5w t4
#define y5w t5
// register `tmp` holds the two temporary values `tx` and `ty`
#define tmp t6
// registers `c0w` to `c7w` hold round constants from array `RCON`
#define c0w s1
#define c1w s2
#define c2w s3
#define c3w s4
#define c4w s5
#define c5w s6
#define c6w s7
#define c7w s8
// register `cptr` holds the start address of array `RCON`
#define cptr s7
// register `scnt` holds the step-counter (for loop termination)
#define scnt s8
// register `trr` serves as temporary register for rotations
#define trr s9

// number of slim and big steps for the SPARKLE384 permutation
.equ STEPS_SLIM, 7
.equ STEPS_BIG, 11


///////////////////////////////////////////////////////////////////////////////
////// MACROS USED BY BOTH SIZE-OPTIMIZED AND SPEED-OPTIMIZED SPARKLE384 //////
///////////////////////////////////////////////////////////////////////////////


// The macro `PROLOGUE` pushes all callee-saved registers on the stack.

.macro PROLOGUE
    addi    sp, sp, -36
    sw      s1, 0(sp)
    sw      s2, 4(sp)
    sw      s3, 8(sp)
    sw      s4, 12(sp)
    sw      s5, 16(sp)
    sw      s6, 20(sp)
    sw      s7, 24(sp)
    sw      s8, 28(sp)
    sw      s9, 32(sp)
.endm


// The macro `EPILOGUE` pops all callee-saved registers from the stack and
// returns to the caller.

.macro EPILOGUE
    lw      s1, 0(sp)
    lw      s2, 4(sp)
    lw      s3, 8(sp)
    lw      s4, 12(sp)
    lw      s5, 16(sp)
    lw      s6, 20(sp)
    lw      s7, 24(sp)
    lw      s8, 28(sp)
    lw      s9, 32(sp)
    addi    sp, sp, 36
    ret
.endm


// The macro `LDSTATE` loads the 12 state-words from RAM and puts them in the
// registers `x0w`-`y5w`. Furthermore, either six or eight round constants are
// loaded from RAM or flash (depending on whether the `RCON` array was placed
// in the .data or .rodata segment) and put into the registers `c0w`-`c5w`
// (resp. `c0w`-`c7w`).

.macro LDSTATE
    // load the 12 state-words
    lw      x0w, 0(sptr)
    lw      y0w, 4(sptr)
    lw      x1w, 8(sptr)
    lw      y1w, 12(sptr)
    lw      x2w, 16(sptr)
    lw      y2w, 20(sptr)
    lw      x3w, 24(sptr)
    lw      y3w, 28(sptr)
    lw      x4w, 32(sptr)
    lw      y4w, 36(sptr)
    lw      x5w, 40(sptr)
    lw      y5w, 44(sptr)
    // load 6 (or 8) RCON-words
    la      cptr, RCON
    lw      c0w, 0(cptr)
    lw      c1w, 4(cptr)
    lw      c2w, 8(cptr)
    lw      c3w, 12(cptr)
    lw      c4w, 16(cptr)
    lw      c5w, 20(cptr)
#if defined(SPARKLE_ASM_FAST)
    lw      c7w, 28(cptr)
    lw      c6w, 24(cptr)
#endif
.endm


// The macro `RORADD` adds two 32-bit operands, whereby one operand is rotated
// `dist` bits right before the addition: bw = bw + (aw >>> dist). Operand `aw`
// is not modified.

.macro RORADD bw:req, aw:req, dist:req
#if (defined(__riscv_zbb) || defined(__riscv_zbkb))
    rori    trr, \aw, \dist
#else  // rotate instructions not supported
    srli    trr, \aw, \dist
    add     \bw, \bw, trr
    slli    trr, \aw, 32-(\dist)
#endif
    add     \bw, \bw, trr
.endm


// The macro `RORXOR` XORs two 32-bit operands, whereby one operand is rotated
// `dist` bits right before the XOR: bw = bw ^ (aw >>> dist). Operand `aw` is
// not modified.

.macro RORXOR bw:req, aw:req, dist:req
#if (defined(__riscv_zbb) || defined(__riscv_zbkb))
    rori    trr, \aw, \dist
#else  // rotate instructions not supported
    srli    trr, \aw, \dist
    xor     \bw, \bw, trr
    slli    trr, \aw, 32-(\dist)
#endif
    xor     \bw, \bw, trr
.endm


// The macro `ARXBOX` computes the ARX-box Alzette.

.macro ARXBOX xi:req, yi:req, ci:req
    RORADD  \xi, \yi, 31    // X = X + (Y >>> 31)
    RORXOR  \yi, \xi, 24    // Y = Y ^ (X >>> 24)
    xor     \xi, \xi, \ci   // X = X ^ RCON
    RORADD  \xi, \yi, 17    // X = X + (Y >>> 17)
    RORXOR  \yi, \xi, 17    // Y = Y ^ (X >>> 17)
    xor     \xi, \xi, \ci   // X = X ^ RCON
    add     \xi, \xi, \yi   // X = X + Y
    RORXOR  \yi, \xi, 31    // Y = Y ^ (X >>> 31)
    xor     \xi, \xi, \ci   // X = X ^ RCON
    RORADD  \xi, \yi, 24    // X = X + (Y >>> 24)
    RORXOR  \yi, \xi, 16    // Y = Y ^ (X >>> 16)
    xor     \xi, \xi, \ci   // X = X ^ RCON
.endm


// The macro `FEIFUNC` computes a part (actually half) of the Feistel function
// $\mathcal{M}_3$ as described in Section 2.1.2 of the SPARKLE specification.
// At first, a temporary value `tmp` is computed, which is done by XORing the
// three left-side words `l0` to `l2`, and thereafter the `ell` operation is
// performed on `tmp`. Then, `ell(tmp)` is XORed to the three right-side words
// `r0` to `r2`, whereby the result of the XOR with `r2` is actually stored in
// register `tmp` (and not in `r2`) to make the subsequent Feistel swap more
// efficient. Note that this macro has to be executed twice: once with x-words
// on the left side and y-words on the right side, and then a second time with
// y-words on the left and x-words on the right.

.macro FEIFUNC l0:req, l1:req, l2:req, r0:req, r1:req, r2:req
    xor     tmp, \l0, \l1
    xor     tmp, tmp, \l2
    slli    trr, tmp, 16
    xor     tmp, tmp, trr
#if (defined(__riscv_zbb) || defined(__riscv_zbkb))
    rori    tmp, tmp, 16
#else  // rotate instructions not supported
    srli    trr, tmp, 16
    slli    tmp, tmp, 16
    xor     tmp, tmp, trr
#endif
    xor     \r0, \r0, tmp
    xor     \r1, \r1, tmp
    xor     tmp, \r2, tmp
.endm


///////////////////////////////////////////////////////////////////////////////
#if !defined(SPARKLE_ASM_FAST) ///// SIZE-OPTIMIZED SPARKLE384 PERMUTATION ////
///////////////////////////////////////////////////////////////////////////////


// The macro `STSTATE` stores the 12 state-words in the registers `x0w`-`y5w`
// to RAM.

.macro STSTATE
    sw      x0w, 0(sptr)
    sw      y0w, 4(sptr)
    sw      x1w, 8(sptr)
    sw      y1w, 12(sptr)
    sw      x2w, 16(sptr)
    sw      y2w, 20(sptr)
    sw      x3w, 24(sptr)
    sw      y3w, 28(sptr)
    sw      x4w, 32(sptr)
    sw      y4w, 36(sptr)
    sw      x5w, 40(sptr)
    sw      y5w, 44(sptr)
.endm


// The macro `ADDRCON` XORs one of the round constants to state-word `y0` and
// the step-counter to state-word `y1`.

.macro ADDRCON
    lw      tmp, 0(cptr)
    add     cptr, cptr, 4
    xor     y1w, y1w, scnt
    xor     y0w, y0w, tmp
.endm


// The macro `ARXLAYER` computes six instances of the ARX-box Alzette.

.macro ARXLAYER
    ARXBOX  x0w, y0w, c0w
    ARXBOX  x1w, y1w, c1w
    ARXBOX  x2w, y2w, c2w
    ARXBOX  x3w, y3w, c3w
    ARXBOX  x4w, y4w, c4w
    ARXBOX  x5w, y5w, c5w
.endm


// The macro `FEISWAP` computes a part (actually half) of the Feistel swap as
// described in Section 2.1.2 of the SPARKLE specification, which means it XORs
// the three left-side words `l0` to `l2` to the three right-side words `r0` to
// `r2` and then swaps the two sides. However, before performing the swap, the
// three right-side words are rotated one word to the left. The right-side word
// `r2` is actually stored in register `tmp`, which means register `r2` is free
// and can be immediately used as target register for the swap. Note that this
// macro has to be executed twice: once with the x-words and then a second time
// with the y-words.

.macro FEISWAP l0:req, l1:req, l2:req, r0:req, r1:req, r2:req
    mv      \r2, \l2
    xor     \l2, \r0, \l0
    mv      \r0, \l0
    xor     \l0, \r1, \l1
    mv      \r1, \l1
    xor     \l1, tmp, \r2
.endm


// The macro `LINLAYER` computes the linear diffusion layer, which consists of
// the Feistel function $\mathcal{M}_3$ and a special Feistel swap (including a
// 1-branch left rotation of the right-side branches). The macro `FEIFUNC` for
// the Feistel function is contained twice: once with x-words on the left side
// and y-words on the right, and then a second time with y-words on the left
// and x-words on the right. Note that the `FEISWAP` after the first `FEIFUNC`
// swaps the left and right y-words, which means the left-side y-words for the
// second `FEIFUNC` are actually on the right side.

.macro LINLAYER
    FEIFUNC  x0w, x1w, x2w, y3w, y4w, y5w
    FEISWAP  y0w, y1w, y2w, y3w, y4w, y5w
    FEIFUNC  y3w, y4w, y5w, x3w, x4w, x5w
    FEISWAP  x0w, x1w, x2w, x3w, x4w, x5w
.endm


// This implementation of the SPARKLE384 permutation is optimized for small
// code size, which means the branch-loop is unrolled, but not the step-loop.

.global sparkle384_rv32
.type sparkle384_rv32, %function
.func sparkle384_rv32
sparkle384_rv32:
    PROLOGUE                // push callee-saved registers
    LDSTATE                 // load the 12 state-words
    mv scnt, zero           // initialize step-counter scnt
.LSTEPLOOP:                 // start of step-loop
    ADDRCON                 // macro for addition of round constant
    ARXLAYER                // macro for complete ARX-box layer
    LINLAYER                // macro for complete linear layer
    addi scnt, scnt, 1      // increment step-counter scnt
    bne scnt, steps, .LSTEPLOOP // branch to start of loop if scnt != steps
    STSTATE                 // store the 12 state-words
    EPILOGUE                // pop callee-saved registers and return
.endfunc


///////////////////////////////////////////////////////////////////////////////
#else ///////////////// SPEED-OPTIMIZED SPARKLE384 PERMUTATION ////////////////
///////////////////////////////////////////////////////////////////////////////


// The macro `STSTATE` stores the 12 state-words in the registers `x0w`-`y5w`
// to RAM. This version differs from the `STSTATE` macro of the size-optimized
// implementation because the left-side words are actually in `x3w`-`y5w` and
// the right-side words in `x0w`-`y2w` (which is due to the fact that the
// speed-optimized implementation does not explicitly swap the left-side and
// right-side words as part of the linear layer).

.macro STSTATE
    sw      x3w, 0(sptr)
    sw      y3w, 4(sptr)
    sw      x4w, 8(sptr)
    sw      y4w, 12(sptr)
    sw      x5w, 16(sptr)
    sw      y5w, 20(sptr)
    sw      x0w, 24(sptr)
    sw      y0w, 28(sptr)
    sw      x1w, 32(sptr)
    sw      y1w, 36(sptr)
    sw      x2w, 40(sptr)
    sw      y2w, 44(sptr)
.endm


// The macro `FEISWAP` computes a part (actually half) of the Feistel swap as
// described in Section 2.1.2 of the SPARKLE specification, which means it XORs
// the three left-side words `l0` to `l2` to the three right-side words `r0` to
// `r2` and rotates the three right-side words one word to the left. Unlike the
// `FEISWAP` macro for the size-optimized implementation, this version does not
// include the swap of the left-side and right-side words. The right-side word
// `r2` is actually stored in register `tmp`, which means register `r2` is free
// and can be immediately used as target register for the rotation. Note that
// this macro has to be executed twice: once with the x-words and then a second
// time with the y-words.

.macro FEISWAP l0:req, l1:req, l2:req, r0:req, r1:req, r2:req
    xor     \r2, \r0, \l0
    xor     \r0, \r1, \l1
    xor     \r1, tmp, \l2
.endm


// The macro `PRMSTEP` computes one full step of the SPARKLE384 permutation
// consisting of RCON addition, ARX-box layer, and linear layer.

.macro PRMSTEP step:req, rcon:req, x0:req, y0:req, x1:req, y1:req, x2:req, \
               y2:req, x3:req, y3:req, x4:req, y4:req, x5:req, y5:req
    // RCON addition
.if (\step != 0)
    xori    \y1, \y1, \step
.endif
    xor     \y0, \y0, \rcon
    // ARX-box Layer
    ARXBOX  \x0, \y0, c0w
    ARXBOX  \x1, \y1, c1w
    ARXBOX  \x2, \y2, c2w
    ARXBOX  \x3, \y3, c3w
    ARXBOX  \x4, \y4, c4w
    ARXBOX  \x5, \y5, c5w
    // Linear Layer
    FEIFUNC \x0, \x1, \x2, \y3, \y4, \y5
    FEISWAP \y0, \y1, \y2, \y3, \y4, \y5
    FEIFUNC \y0, \y1, \y2, \x3, \x4, \x5
    FEISWAP \x0, \x1, \x2, \x3, \x4, \x5
.endm


// The macro `RETSLIM` checks whether the number of steps (passed as parameter
// `steps` to the function) equals `STEPS_SLIM`. If this is the case then the
// 12 state-words are stored to RAM and the function returns to the caller.

.macro RETSLIM
    addi    steps, steps, -(STEPS_SLIM)
    bne     steps, zero, .+96
    // store the 12 state-words
    STSTATE
    // pop callee-saved registers
    EPILOGUE
.endm


// This implementation of the SPARKLE384 permutation is optimized for high
// speed, which means both the branch-loop and the step-loop are unrolled.
// Unlike the size-optimized variant, this speed-optimized implementation does
// not explicitly perform the swap of the left and right-side branches as part
// of the linear layer. Instead, this swap is carried out by adapting the order
// of the 12 state-word-registers that are passed as 3rd to 14th argument to
// the `PRMSTEP` macro.

.global sparkle384_rv32
.type sparkle384_rv32, %function
.func sparkle384_rv32
sparkle384_rv32:
    PROLOGUE                // push callee-saved registers
    LDSTATE                 // load the 12 state-words
    PRMSTEP  0, c0w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    PRMSTEP  1, c1w, x3w, y3w, x4w, y4w, x5w, y5w, x0w, y0w, x1w, y1w, x2w, y2w
    PRMSTEP  2, c2w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    PRMSTEP  3, c3w, x3w, y3w, x4w, y4w, x5w, y5w, x0w, y0w, x1w, y1w, x2w, y2w
    PRMSTEP  4, c4w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    PRMSTEP  5, c5w, x3w, y3w, x4w, y4w, x5w, y5w, x0w, y0w, x1w, y1w, x2w, y2w
    PRMSTEP  6, c6w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    RETSLIM                 // return when the number of steps is slim
    PRMSTEP  7, c7w, x3w, y3w, x4w, y4w, x5w, y5w, x0w, y0w, x1w, y1w, x2w, y2w
    PRMSTEP  8, c0w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    PRMSTEP  9, c1w, x3w, y3w, x4w, y4w, x5w, y5w, x0w, y0w, x1w, y1w, x2w, y2w
    PRMSTEP 10, c2w, x0w, y0w, x1w, y1w, x2w, y2w, x3w, y3w, x4w, y4w, x5w, y5w
    STSTATE                 // store the 12 state-words
    EPILOGUE                // pop callee-saved registers and return
.endfunc


#endif  // #if !defined(SPARKLE_ASM_FAST)


.size sparkle384_rv32, .-sparkle384_rv32


///////////////////////////////////////////////////////////////////////////////
/////////////////////////// SPARKLE ROUND CONSTANTS ///////////////////////////
///////////////////////////////////////////////////////////////////////////////

// This implementation places the round constants in the .data segment, which
// means they are loaded from RAM during the computation of the ARX-boxes. It
// would also be possible to place them in the .rodata segment (by replacing
// the ".section .data" directive below by ".section .rodata") so that they are
// loaded from flash, which reduces the RAM consumption by 44 bytes, but may
// increase the execution time on devices with a high number of flash wait
// states.

.section .data
.balign 4


.type RCON, %object
.size RCON, 44
RCON:
.word 0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738
.word 0xBB1185EB, 0x4F7C7B57, 0xCFBFA1C8, 0xC2B3293D
.word 0xB7E15162, 0xBF715880, 0x38B4DA56


.end
