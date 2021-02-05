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
 * sparkle_best_trail_search (serial version)
 *
 * Example command-line execution: rounds, rot const, bounds
 *
 * make && ./main_serial 6 31 17 0 24 24 17 31 16 0 -1 -2 -6 -10 -18
 *
 * or
 * 
 * g++ -std=c++11 -o main_serial main_serial.cpp trail.cpp sparkle_best_trail_search.cpp -I.
 *
 ******************************************************************************/

#include <cstdio>
#include <cmath>
#include <csignal>
#include <chrono>
#include <unistd.h>
#include <cassert>

#define MAIN 1
#include "defs.h"
#include "trail.h"
#include "sparkle_best_trail_search.h"


/******************************************************************************
 * Signal handlers to catch ctrl-c
 ******************************************************************************/
void signal_handler(int signal)
{
	exit(-1);
}

/******************************************************************************
 * Helper
 ******************************************************************************/
void finalize(std::chrono::microseconds startTime)
{
	std::chrono::microseconds stopTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch());
	unsigned long int duration = stopTime.count() - startTime.count();
	double speed;
	if (duration != 0)
	{
		speed = double(nNodes)/double(duration);
	}
	else
	{
		speed = nan("");
	}
	printf("[%g s] {%lu nodes -> %g Mnodes/s}\n", (double)duration/1e6, nNodes, speed);
	#ifdef STATS
	std::string filename = "nodes_" + std::to_string(g_Bn) + ".stats";
	FILE *fh = fopen(filename.c_str(), "w");
	fprintf(fh, "# Nodes per round:\n");
	for (unsigned int i = 0; i < NROUNDS; i++)
	{
		fprintf(fh, "%u,%lu\n", i, nNodesPerRound[i]);
	}
	fprintf(fh, "# Nodes per round per ibit:\n");
	for (unsigned int i = 0; i < NROUNDS; i++)
	{
		for (unsigned int j = 0; j < WORD_SIZE + 1; j++)
		{
			fprintf(fh, "%u,%u,%lu\n", i, j, nNodesPerRoundPerIbit[i][j]);
		}
	}
	fprintf(fh, "# Paths per round per ibit:\n");
	for (unsigned int i = 0; i < NROUNDS; i++)
	{
		for (unsigned int j = 0; j < WORD_SIZE + 1; j++)
		{
			fprintf(fh, "%u,%u,%lu\n", i, j, nPathsPerRoundPerIbit[i][j]);
		}
	}
	fclose(fh);
	#endif
}

/******************************************************************************
 * Main
 ******************************************************************************/
int main(int argc, char *argv[])
{
	std::signal(SIGINT, signal_handler);

	
	bool benchMode = false;

	if (argc <= 2) {
		printf("Usage: %s <NROUNDS> <R1> <R2> ... <S1> <S2> ...\n", argv[0]);
		return -1;
	}
	NROUNDS = atoi(argv[1]);
	printf("NROUNDS %d\n", NROUNDS);

	for(int i = 0; i < NROUNDS; i++) {
		gconst_r[i] = atoi(argv[2+i%4]);
		gconst_s[i] = atoi(argv[2+4+i%4]);
		printf("r%d = %d, s%d = %d\n", i, gconst_r[i], i, gconst_s[i]);
	}
	for(int i = 0; i < NROUNDS-1; i++) {
		g_best_B[i] = atoi(argv[2+8+i]);
		printf("g_best_B[%d] = %d\n", i, g_best_B[i]);
	}

	g_Bn = g_best_B[NROUNDS - 2];
	printf("-- g_Bn forced to %d\n", g_Bn);
#if 0 // DEBUG
	if(NROUNDS == 6) { // for 6R there are no trails with p > 2^-18
	  g_Bn = -18;
	  printf("-- g_Bn forced to %d\n", g_Bn);
	}
#endif // #if 1 // DEBUG

	printf("-- Searching for %d rounds (WORD_SIZE %d bits):\n", NROUNDS, WORD_SIZE);

	#ifdef STATS
	printf("-- compiled with STATS\n");
	#endif

	std::chrono::microseconds initialTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch());
	while (1)
	{
		fprintf(stdout, "-- g_Bn = %+3d ... ", g_Bn);
		fflush(stdout);
		std::chrono::microseconds startTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch());
		bool res = sparkle_best_trail_search();
		if (res == false)
		{
			printf("no trail found ");
		}
		else
		{
			printf("trail found ");
		}
		finalize(startTime);
		if (res == true)
		{
			fprintTrail(stdout, trail);
			break; // job done!
		}
		if (benchMode == true)
		{
			break;
		}
		g_Bn--;
#if 0 // DEBUG
		if(g_Bn < -19) { // find only trails of prob>= 2^-19
		  break;
		}
#endif // #if 1 // DEBUG
	}

	//std::chrono::microseconds finalTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch());
	//printf("== Total duration: %g s\n", (double)(finalTime.count() - initialTime.count())/1.0e6);

	return 0;
}
