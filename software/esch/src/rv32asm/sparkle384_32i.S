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

// #define SPARKLE_ASM_FAST --> NOT YET IMPLEMENTED !!!


.section .text
.balign 4


///////////////////////////////////////////////////////////////////////////////
//////////////////////// REGISTER NAMES AND CONSTANTS /////////////////////////
///////////////////////////////////////////////////////////////////////////////

// register sptr holds the start address of array `state`
#define sptr a0
// register steps holds the number of steps (parameter `steps`)
#define steps a1
// registers x0w to y5w hold 12 words from array `state`
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
// register txy holds the temporary values `tx` and `ty`
#define txy t6
// registers c0w to c5w hold round constants from array `rcon`
#define c0w s1
#define c1w s2
#define c2w s3
#define c3w s4
#define c4w s5
#define c5w s6
// register cptr holds the start address of array `rcon`
#define cptr s7
// register scnt holds the step-counter (for loop termination)
#define scnt s8
// register trr serves as temporary register for rotations
#define trr s9


///////////////////////////////////////////////////////////////////////////////
//////////////////////////// MACROS FOR SPARKLE384 ////////////////////////////
///////////////////////////////////////////////////////////////////////////////


.macro PROLOGUE_384
    // allocate space on stack
    addi    sp, sp, -36  
    // push nine callee-saved registers
    sw      s1, 0(sp)
    sw      s2, 4(sp)
    sw      s3, 8(sp)
    sw      s4, 12(sp)
    sw      s5, 16(sp)
    sw      s6, 20(sp)
    sw      s7, 24(sp)
    sw      s8, 28(sp)
    sw      s9, 32(sp)
    // load the 12 state-words x0-y5
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
    // load the 6 round constants c0-c5
    la      cptr, RCON 
    lw      c0w, 0(cptr)
    lw      c1w, 4(cptr)
    lw      c2w, 8(cptr)
    lw      c3w, 12(cptr)
    lw      c4w, 16(cptr)
    lw      c5w, 20(cptr)
.endm


.macro EPILOGUE_384
    // store the 12 state-words x0-y5
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
    // pop nine callee-saved registers
    lw      s1, 0(sp)
    lw      s2, 4(sp)
    lw      s3, 8(sp)
    lw      s4, 12(sp)
    lw      s5, 16(sp)
    lw      s6, 20(sp)
    lw      s7, 24(sp)
    lw      s8, 28(sp)
    lw      s9, 32(sp)
    // deallocate space on stack
    addi    sp, sp, 36
    ret
.endm


// addition of step-counter
.macro ADD_STEP_CNT_384
    lw      txy, 0(cptr)
    add     cptr, cptr, 4
    xor     y1w, y1w, scnt
    xor     y0w, y0w, txy
.endm


// operation B = B + (A >>> dist)
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


// operation B = B ^ (A >>> dist)
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


// complete ARX-box (Alzette)
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


// complete non-linear layer
.macro ARXBOX_LAYER_384
    ARXBOX  x0w, y0w, c0w
    ARXBOX  x1w, y1w, c1w
    ARXBOX  x2w, y2w, c2w
    ARXBOX  x3w, y3w, c3w
    ARXBOX  x4w, y4w, c4w
    ARXBOX  x5w, y5w, c5w
.endm


// half of Feistel function
.macro HFFUNC l0:req, l1:req, l2:req, r0:req, r1:req, r2:req
    xor     txy, \l0, \l1
    xor     txy, txy, \l2
    slli    trr, txy, 16
    xor     txy, txy, trr
#if (defined(__riscv_zbb) || defined(__riscv_zbkb))
    rori    txy, txy, 16
#else  // rotate instructions not supported
    srli    trr, txy, 16
    slli    txy, txy, 16
    xor     txy, txy, trr
#endif
    xor     \r0, \r0, txy
    xor     \r1, \r1, txy
    xor     txy, \r2, txy
.endm


// half of Feistel XOR+swap
.macro HFSWAP l0:req, l1:req, l2:req, r0:req, r1:req, r2:req
    mv      \r2, \l2
    xor     \l2, \r0, \l0
    mv      \r0, \l0
    xor     \l0, \r1, \l1
    mv      \r1, \l1
    xor     \l1, txy, \r2
.endm


// complete linear layer
.macro LINEAR_LAYER_384
    HFFUNC  x0w, x1w, x2w, y3w, y4w, y5w
    HFSWAP  y0w, y1w, y2w, y3w, y4w, y5w
    HFFUNC  y3w, y4w, y5w, x3w, x4w, x5w
    HFSWAP  x0w, x1w, x2w, x3w, x4w, x5w
.endm


///////////////////////////////////////////////////////////////////////////////
////////////////// SPARKLE384 PERMUTATION (BRANCH-UNROLLED) ///////////////////
///////////////////////////////////////////////////////////////////////////////


.global sparkle384_rv32
.type sparkle384_rv32, %function
.func sparkle384_rv32
sparkle384_rv32:
    PROLOGUE_384            // push callee-saved registers and load state
    mv scnt, zero           // initialize step-counter scnt
.LSTEPLOOP:                 // start of step-loop
    ADD_STEP_CNT_384        // macro to add step-counter and rcon to state
    ARXBOX_LAYER_384        // macro for the ARXBOX layer
    LINEAR_LAYER_384        // macro for the linear layer
    addi scnt, scnt, 1      // increment step-counter scnt
    bne scnt, steps, .LSTEPLOOP // branch to start of loop if scnt != steps
    EPILOGUE_384            // store state and pop callee-saved registers
.endfunc


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
