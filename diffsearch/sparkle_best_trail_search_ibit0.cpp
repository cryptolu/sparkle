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
bool sparkle_best_trail_search_ibit0(unsigned int iround, WORD_T alpha, WORD_T beta)
{
  /* NB: This function is called when iround != NROUNDS && ibit == 0. Under those conditions, we have gamma = 0 */
  /* Iteration #0 */
  // printf("%d %08x %08x\n", iround, alpha, beta);
  nNodes++;
#ifdef STATS
  nNodesPerRound[iround - 1]++;
  nNodesPerRoundPerIbit[iround - 1][0]++;
  bool gotChildren_0 = false;
#endif
  int p_trail = g_best_B[NROUNDS - iround - 1];
  p_trail += trail[iround - 2].cp;
  for (unsigned int w_0 = 0; w_0 < 2; w_0++)
    {
      const WORD_T gamma_part_0 = (w_0 << 0);
      int p_part_0 = xdp_add_lm(alpha, beta, gamma_part_0, 1);
      int p_est_0 = p_trail + p_part_0;
      if (p_est_0 >= g_Bn)
	{
	  /* Iteration #1 */
	  nNodes++;
#ifdef STATS
	  gotChildren_0 = true;
	  nNodesPerRound[iround - 1]++;
	  nNodesPerRoundPerIbit[iround - 1][1]++;
	  bool gotChildren_1 = false;
#endif
	  for (unsigned int w_1 = 0; w_1 < 3; w_1 += 2)
	    {
	      const WORD_T gamma_part_1 = gamma_part_0 | w_1;
	      int p_part_1 = xdp_add_lm(alpha, beta, gamma_part_1, 2);
	      int p_est_1 = p_trail + p_part_1;
	      if (p_est_1 >= g_Bn)
		{
		  /* Iteration #2 */
		  nNodes++;
#ifdef STATS
		  gotChildren_1 = true;
		  nNodesPerRound[iround - 1]++;
		  nNodesPerRoundPerIbit[iround - 1][2]++;
		  bool gotChildren_2 = false;
#endif
		  for (unsigned int w_2 = 0; w_2 < 7; w_2 += 4)
		    {
		      const WORD_T gamma_part_2 = gamma_part_1 | w_2;
		      int p_part_2 = xdp_add_lm(alpha, beta, gamma_part_2, 3);
		      int p_est_2 = p_trail + p_part_2;
		      if (p_est_2 >= g_Bn)
			{
			  /* Iteration #3 */
			  nNodes++;
#ifdef STATS
			  gotChildren_2 = true;
			  nNodesPerRound[iround - 1]++;
			  nNodesPerRoundPerIbit[iround - 1][3]++;
			  bool gotChildren_3 = false;
#endif
			  for (unsigned int w_3 = 0; w_3 < 15; w_3 += 8)
			    {
			      const WORD_T gamma_part_3 = gamma_part_2 | w_3;
			      int p_part_3 = xdp_add_lm(alpha, beta, gamma_part_3, 4);
			      int p_est_3 = p_trail + p_part_3;
			      if (p_est_3 >= g_Bn)
				{
				  /* Iteration #4 */
				  nNodes++;
#ifdef STATS
				  gotChildren_3 = true;
				  nNodesPerRound[iround - 1]++;
				  nNodesPerRoundPerIbit[iround - 1][4]++;
				  bool gotChildren_4 = false;
#endif
				  for (unsigned int w_4 = 0; w_4 < 31; w_4 += 16)
				    {
				      const WORD_T gamma_part_4 = gamma_part_3 | w_4;
				      int p_part_4 = xdp_add_lm(alpha, beta, gamma_part_4, 5);
				      int p_est_4 = p_trail + p_part_4;
				      if (p_est_4 >= g_Bn)
					{
					  /* Iteration #5 */
					  nNodes++;
#ifdef STATS
					  gotChildren_4 = true;
					  nNodesPerRound[iround - 1]++;
					  nNodesPerRoundPerIbit[iround - 1][5]++;
					  bool gotChildren_5 = false;
#endif
					  for (unsigned int w_5 = 0; w_5 < 63; w_5 += 32)
					    {
					      const WORD_T gamma_part_5 = gamma_part_4 | w_5;
					      int p_part_5 = xdp_add_lm(alpha, beta, gamma_part_5, 6);
					      int p_est_5 = p_trail + p_part_5;
					      if (p_est_5 >= g_Bn)
						{
						  /* Iteration #6 */
						  nNodes++;
#ifdef STATS
						  gotChildren_5 = true;
						  nNodesPerRound[iround - 1]++;
						  nNodesPerRoundPerIbit[iround - 1][6]++;
						  bool gotChildren_6 = false;
#endif
						  for (unsigned int w_6 = 0; w_6 < 127; w_6 += 64)
						    {
						      const WORD_T gamma_part_6 = gamma_part_5 | w_6;
						      int p_part_6 = xdp_add_lm(alpha, beta, gamma_part_6, 7);
						      int p_est_6 = p_trail + p_part_6;
						      if (p_est_6 >= g_Bn)
							{
							  /* Iteration #7 */
							  nNodes++;
#ifdef STATS
							  gotChildren_6 = true;
							  nNodesPerRound[iround - 1]++;
							  nNodesPerRoundPerIbit[iround - 1][7]++;
							  bool gotChildren_7 = false;
#endif
							  for (unsigned int w_7 = 0; w_7 < 255; w_7 += 128)
							    {
							      const WORD_T gamma_part_7 = gamma_part_6 | w_7;
							      int p_part_7 = xdp_add_lm(alpha, beta, gamma_part_7, 8);
							      int p_est_7 = p_trail + p_part_7;
							      if (p_est_7 >= g_Bn)
								{
#ifdef STATS
								  gotChildren_7 = true;
#endif
								  bool res = sparkle_best_trail_search_i(iround, 8, alpha, beta, gamma_part_7);
								  if (res == true)
								    {
								      return true;
								    }
								}
							    }
#ifdef STATS
							  if (gotChildren_7 == false)
							    {
							      nPathsPerRoundPerIbit[iround -1][7]++;
							    }
#endif
							}
						    }
#ifdef STATS
						  if (gotChildren_6 == false)
						    {
						      nPathsPerRoundPerIbit[iround -1][6]++;
						    }
#endif
						}
					    }
#ifdef STATS
					  if (gotChildren_5 == false)
					    {
					      nPathsPerRoundPerIbit[iround -1][5]++;
					    }
#endif
					}
				    }
#ifdef STATS
				  if (gotChildren_4 == false)
				    {
				      nPathsPerRoundPerIbit[iround -1][4]++;
				    }
#endif
				}
			    }
#ifdef STATS
			  if (gotChildren_3 == false)
			    {
			      nPathsPerRoundPerIbit[iround -1][3]++;
			    }
#endif
			}
		    }
#ifdef STATS
		  if (gotChildren_2 == false)
		    {
		      nPathsPerRoundPerIbit[iround -1][2]++;
		    }
#endif
		}
	    }
#ifdef STATS
	  if (gotChildren_1 == false)
	    {
	      nPathsPerRoundPerIbit[iround -1][1]++;
	    }
#endif
	}
    }
#ifdef STATS
  if (gotChildren_0 == false)
    {
      nPathsPerRoundPerIbit[iround -1][0]++;
    }
#endif
  return false;
}
