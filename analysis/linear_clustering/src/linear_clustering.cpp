/*
 * Experiments on linear clustering of the double ARXbox of SPARKLE.
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
 * linear_clustering.c
 *
 * Compile as: g++ -O3 -std=c++11 -Wall linear_clustering.cpp -o linear_clustering -lpthread
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include <string.h>
#include "approx_8rnds.cpp"
#include <pthread.h>

#define SAMPLES (1ull << 40) // 2^40
#define NTHREADS 50


char in_parity[N_IN] = {};
char out_parity[N_OUT] = {};
uint64_t hits[N_IN][N_OUT] = {};


uint32_t rotr(uint32_t x, int r) {
    return (x >> r) | (x << (32-r));
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
    y = ARXround(y,  0, 31, c);
    y = ARXround(y, 24, 16, c);
    return y;
}

// parity of x
int parity(uint64_t x) {
  return __builtin_popcountll(x)&1;
}

void account(uint64_t x, uint64_t y) {
    for (int in=0; in<N_IN; in++) {
            in_parity[in] = parity(IN_MASKS[in]&x);
    }
    for (int out=0; out<N_OUT; out++) {
            out_parity[out] = parity(OUT_MASKS[out]&y);
    }
    for (int in=0; in<N_IN; in++) {
        for (int out=0; out<N_OUT; out++) {
            hits[in][out] += (in_parity[in] == out_parity[out]);
        }
    }
}

uint64_t c_a;
uint64_t c_b;

void* worker(void *arg) {
    // the random number generator for uniform 64bit inputs
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    uint64_t x, y;
    uint64_t thread_hits[N_IN][N_OUT] = {};
    char in_parity[N_IN];
    char out_parity[N_OUT];
    for (uint64_t i=0; i<SAMPLES / NTHREADS; i++) {
        x = dis(gen);
        y = ARXbox(x, c_a);
        y = ARXbox(y, c_b);

        for (int in=0; in<N_IN; in++) {
            in_parity[in] = parity(IN_MASKS[in]&x);
        }
        for (int out=0; out<N_OUT; out++) {
            out_parity[out] = parity(OUT_MASKS[out]&y);
        }
        for (int in=0; in<N_IN; in++) {
            for (int out=0; out<N_OUT; out++) {
                thread_hits[in][out] += (in_parity[in] == out_parity[out]);
            }
        }
    } // end of main loop

    // merge thread results
    for (int in=0; in<N_IN; in++) {
        for (int out=0; out<N_OUT; out++) {
            hits[in][out] += thread_hits[in][out];
        }
    }
    return NULL;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        printf("Usage: %s <first_constant> <second_constant> ...\n", argv[0]);
        return -1;
    }
    int first = atoi(argv[1]);
    int second = atoi(argv[2]);


    // round constants for the ARX-boxes
    uint64_t c[8] = {
    0xB7E1516200000000,
    0xBF71588000000000,
    0x38B4DA5600000000,
    0x324E773800000000,
    0xBB1185EB00000000,
    0x4F7C7B5700000000,
    0xCFBFA1C800000000,
    0xC2B3293D00000000};

    c_a = c[first];
    c_b = c[second];

    printf("Samples: 0x%llx\n", SAMPLES);

    printf("Round constants:\n");
    printf("c_a: 0x%jx\n", c_a);
    printf("c_b: 0x%jx\n", c_b);

    memset(hits,0,sizeof(hits));

    int rc;
    pthread_t threads[NTHREADS];
    pthread_attr_t attr;
    void *status;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int i = 0; i < NTHREADS; i++) {
        printf("creating thread %d\n", i);
        rc = pthread_create(&threads[i], &attr, worker, NULL);
        if (rc) {
            printf("error create\n");
            exit(-1);
        }
    }
    for (int i = 0; i < NTHREADS; i++) {
        rc = pthread_join(threads[i], &status);
        if (rc) {
            printf("error join\n");
            exit(-1);
        }
        printf("thread %d finished\n", i);
    }

    // write correlations to file
    char filename[400];
    FILE *fp;
    sprintf(filename, "linear_clustering_8rnds_Ac%dAc%d.dat", first,second);
    fp = fopen(filename, "w");
    int current_input, current_output;
    double cor = 0.0;

    for (int ilin=0; ilin<NLIN; ilin++) {
        fprintf(fp,"%2d: R#%1d 0x%016lx -> 0x%016lx 2^%d\n", ilin, (int)APPROXIMATIONS[ilin][0], APPROXIMATIONS[ilin][1], APPROXIMATIONS[ilin][2], -(int)APPROXIMATIONS[ilin][3]);
        // find approximations in table
        current_input = -1;
        current_output = -1;
        for (int in=0; in<N_IN; in++) {
                if (IN_MASKS[in]==APPROXIMATIONS[ilin][1])
                        current_input = in;
        }
        for (int out=0; out<N_OUT; out++) {
                if (OUT_MASKS[out]==APPROXIMATIONS[ilin][2])
                        current_output = out;
        }
        // compute correlation
        cor = 2.0*((double)hits[current_input][current_output]/SAMPLES) - 1.0;
        fprintf(fp,"cor: %4.10f\n",cor);
        if (cor < 0) {
                fprintf(fp, "-log2 (abs): %4.2f\n", -log2(-cor));
        }
        else {
                fprintf(fp, "-log2 (abs): %4.2f\n", -log2(cor));
        }

    }

    fclose(fp);

    return 0;
}
