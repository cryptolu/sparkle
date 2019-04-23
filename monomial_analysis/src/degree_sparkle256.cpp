/*
 * Experiments on the monomial distribution in SPARKLE.
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
 * degree_sparkle256.c
 *
 * Compile as: g++ -std=c++11 -Wall degree_sparkle256.cpp -o degree_sparkle256
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include <cassert>
#include <string.h>
#include <iostream>

using namespace std;

#define N_IN 256
#define N_OUT 256

// here, define the degree of the monomials to be checked and the number of samples 
// for degree 1 and 2, all monomials will be checked and there is no sampling
#define DEGREE 3
#define SAMPLES 10000

// define the number of steps to check
#define STEPS 4

uint8_t input[N_IN];
uint8_t output[N_OUT];
uint8_t monomial[N_IN];
uint64_t sol[N_OUT];

// random number generator to generate integer uniformly in [0,N_IN)
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_int_distribution<int> dis(0,N_IN-1);

// round constants
uint64_t c0 = 0xB7E1516200000000;
uint64_t c1 = 0xBF71588000000000;
uint64_t c2 = 0x38B4DA5600000000;
uint64_t c3 = 0x324E773800000000;
uint64_t c4 = 0xBB1185EB00000000;
uint64_t c5 = 0x4F7C7B5700000000;
uint64_t c6 = 0xCFBFA1C800000000;
uint64_t c7 = 0xC2B3293D00000000;

// randomly sets the monomial array at degree positions
int random_monomial(int degree) {
    memset(monomial,0,sizeof(monomial));
    int var[degree];
    memset(var,-1,sizeof(var));
    int i = 0;
    int v = 0;
    int stop;
    while(i < degree) {
        v = dis(gen);
        // check if v already in var
        stop = 0;
        for (int j=0; j<degree; j++) {
            if (v==var[j])
                stop = 1;
        }
        if (stop==0) {
            var[i] = v;
            i++;
        }
    }
    for (int j=0; j<degree; j++) {
        monomial[var[j]] = 1;
    }
    return 0;
}

uint32_t rotr(uint32_t x, int r) {
    return (x >> r % 32) | (x << (32-r) % 32);
}

uint64_t ARXround(uint64_t x, int r, int s, uint64_t c) {
    uint32_t xr = x & 0xFFFFFFFF;
    uint32_t xl = (x >> 32) & 0xFFFFFFFF;
    xl = xl + rotr(xr,r);
    xr = xr ^ rotr(xl,s);
    return (((uint64_t)(xl) << 32) ^ (uint64_t)(xr) ^ c);
}

uint64_t ARXbox(uint64_t x, uint64_t c) {
  uint64_t y = x;
  y = ARXround(y, 31, 24, c);
  y = ARXround(y, 17, 17, c);
  y = ARXround(y, 0, 31, c);
  y = ARXround(y, 24, 16, c);
  return y;
}



// the function for which the ANF is checked (here Sparkle256)
int F() {
	uint64_t x0 = 0;
	uint64_t x1 = 0;
	uint64_t x2 = 0;
	uint64_t x3 = 0;
	uint32_t x0r = 0;
    uint32_t x0l = 0;
    uint32_t x1r = 0;
    uint32_t x1l = 0;
    uint32_t x2r = 0;
    uint32_t x2l = 0;
    uint32_t x3r = 0;
    uint32_t x3l = 0;
	uint32_t tx = 0;
	uint32_t ty = 0;
	uint32_t temp0,temp1;
	for (int i=0; i<64; i++) {
		x0 = x0 ^ ((uint64_t)input[i]<<(63-i));
		x1 = x1 ^ ((uint64_t)input[i+64]<<(63-i));
		x2 = x2 ^ ((uint64_t)input[i+2*64]<<(63-i));
		x3 = x3 ^ ((uint64_t)input[i+3*64]<<(63-i));
	}
	for (int s=0; s<STEPS; s++)
	{
        if ((s==0) % 8)
            x0 = x0^(c0>>32);
        if ((s==1) % 8)
            x0 = x0^(c1>>32);
        if ((s==2) % 8)
            x0 = x0^(c2>>32);
        if ((s==3) % 8)
            x0 = x0^(c3>>32);
        if ((s==4) % 8)
            x0 = x0^(c4>>32);
        if ((s==5) % 8)
            x0 = x0^(c5>>32);
        if ((s==6) % 8)
            x0 = x0^(c6>>32);
        if ((s==7) % 8)
            x0 = x0^(c7>>32);

        x1 = x1^s;

        x0 = ARXbox(x0,c0);
        x1 = ARXbox(x1,c1);
        x2 = ARXbox(x2,c2);
        x3 = ARXbox(x3,c3);

        // linear layer
        x0r = x0 & 0xFFFFFFFF;
        x0l = (x0 >> 32) & 0xFFFFFFFF;
        x1r = x1 & 0xFFFFFFFF;
        x1l = (x1 >> 32) & 0xFFFFFFFF;
        x2r = x2 & 0xFFFFFFFF;
        x2l = (x2 >> 32) & 0xFFFFFFFF;
        x3r = x3 & 0xFFFFFFFF;
        x3l = (x3 >> 32) & 0xFFFFFFFF;

        tx = x0l^x1l;
        ty = x0r^x1r;
        tx = rotr(tx^(tx<<16),16);
        ty = rotr(ty^(ty<<16),16);

        x2r=x2r^x0r^tx;
        x3r=x3r^x1r^tx;
        x2l=x2l^x0l^ty;
        x3l=x3l^x1l^ty;

        temp0=x0l;
        temp1=x1l;
        x0l=x3l;
        x1l=x2l;
        x2l=temp0;
        x3l=temp1;

        temp0=x0r;
        temp1=x1r;
        x0r=x3r;
        x1r=x2r;
        x2r=temp0;
        x3r=temp1;

        x0= ((uint64_t)(x0l) << 32) ^ (uint64_t)(x0r);
        x1= ((uint64_t)(x1l) << 32) ^ (uint64_t)(x1r);
        x2= ((uint64_t)(x2l) << 32) ^ (uint64_t)(x2r);
        x3= ((uint64_t)(x3l) << 32) ^ (uint64_t)(x3r);
    }

	for (int i=0; i<64; i++) {
		output[i] = (x0 >> (63-i))&0x1;
		output[i+64] = (x1 >> (63-i))&0x1;
		output[i+2*64] = (x2 >> (63-i))&0x1;
		output[i+3*64] = (x3 >> (63-i))&0x1;
	}

	return 0;
}

// checks if the monomial set in array monomial occurs in the ANF (for all output bits). sol[l] will be set to 1 iff monomial occurs in the l-th output bit
// sol[l] is 0 otherwise. It uses the Möbius transform to check the occurence of the monomial
int check_monomial() {
    int j=0;
    memset(sol, 0, sizeof(sol));
    for (int i=0; i<N_IN; i++) {
        if (monomial[i]==1) {
            printf("x%d", i);
            j = j+1;
        }
    }
    printf("\n");

    int pos[j];
    int k=0;
    for (int i=0; i<N_IN; i++) {
        if (monomial[i]==1) {
            pos[k] = i;
            k=k+1;
        }
    }

    for (uint64_t t=0; t<((uint64_t)1<<j); t++)
    {
        memset(input, 0, sizeof(input));
        for (int i=0; i<j; i++) {
            input[pos[i]] = (t>>i)&(0x1);
        }
        F();
        for (int i=0; i<N_OUT; i++) {
            if (output[i] == 1) {
                sol[i]++;
            }
        }
    }
    for (int i=0; i<N_OUT; i++) {
        sol[i] = (sol[i] % 2);
    }

    return 0;

}

int main(int argc, char* argv[])
{
    int c_min,c_max;

   char filename[400];
   FILE *fp;
   sprintf(filename, "sparkle256_monomials%d_%d.dat", DEGREE, STEPS);
   fp = fopen(filename, "w");
   int nb[N_OUT];
   memset(nb,0,sizeof(nb));


   #if DEGREE==1
   memset(nb,0,sizeof(nb));
   for (int i=0; i<N_IN; i++) {
        memset(monomial,0,sizeof(monomial));
        monomial[i] = 1;
        check_monomial();
        for (int l=0; l<N_OUT; l++)
        {
            fprintf(fp, "%d", (int)sol[l]);
            if (sol[l]==1)
                nb[l]++;
        }
        fprintf(fp,"\n");
    }
    c_min = 200000000;
    c_max = 0;
    for (int l=0; l<N_OUT;l++) {
        if (nb[l] < c_min)
            c_min = nb[l];
        if (nb[l] > c_max)
            c_max = nb[l];
        fprintf(fp, "%d\t", nb[l]);
    }
    fprintf(fp, "\nmin: %d", c_min);
    fprintf(fp, "\nmax: %d", c_max);
    #endif

    #if DEGREE==2
    memset(nb,0,sizeof(nb));
   for (int i=0; i<N_IN-1; i++) {
    for (int j=i+1; j<N_IN; j++) {
        memset(monomial,0,sizeof(monomial));
        monomial[i] = 1;
        monomial[j] = 1;
        check_monomial();
        for (int l=0; l<N_OUT; l++)
        {
            fprintf(fp, "%d", (int)sol[l]);
            if (sol[l]==1)
                nb[l]++;
        }
        fprintf(fp,"\n");
        }
    }
    c_min = 200000000;
    c_max = 0;
    for (int l=0; l<N_OUT;l++) {
        if (nb[l] < c_min)
            c_min = nb[l];
        if (nb[l] > c_max)
            c_max = nb[l];
        fprintf(fp, "%d\t", nb[l]);
    }
    fprintf(fp, "\nmin: %d", c_min);
    fprintf(fp, "\nmax: %d", c_max);
   #endif


   #if DEGREE >2
   for (int i=0; i<SAMPLES; i++) {
        random_monomial(DEGREE);
        for (int k=0; k<N_IN; k++) {
            if (monomial[k]==1) {
                fprintf(fp, "x%d", k);
            }
        }
        fprintf(fp,"\n");
        check_monomial();
        for (int l=0; l<N_OUT; l++)
        {
            fprintf(fp, "%d", (int)sol[l]);
            if (sol[l]==1)
                nb[l]++;
        }
        fprintf(fp,"\n");
        for (int l=0; l<N_OUT;l++) {
            printf("%d\t", nb[l]);
        }
        printf("\n");
    }
    c_min = 200000000;
    c_max = 0;
    for (int l=0; l<N_OUT;l++) {
        if (nb[l] < c_min)
            c_min = nb[l];
        if (nb[l] > c_max)
            c_max = nb[l];
        fprintf(fp, "%d\t", nb[l]);
    }
    fprintf(fp, "\nmin: %d", c_min);
    fprintf(fp, "\nmax: %d", c_max);
   #endif

   fclose(fp);

    return 0;
}
