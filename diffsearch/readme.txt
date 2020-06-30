Search for the optimal differential trails of the ARXbox of SPARKLE

- Compile:

make

or

g++ -std=c++11 -O3 -o main_serial main_serial.cpp trail.cpp sparkle_best_trail_search.cpp -I.

- Clean object files:

make clean

- Execute:

The command-line syntax is as follows

./diffsearch <number of rounds R> <the 8 rotation constants of the ARXbox> <log2 probabilities of the best trails for rounds 1, 2, ... R-1>

- Examples:

-- The following command searches for the best differential trails for 4 rounds

./diffsearch 4 31 17 0 24 24 17 31 16 0 -1 -2

where  the command-line arguments are to be interpreted as follows:

4 = number of rounds

31 17 0 24 24 17 31 16 = the 8 rotation constants of the ARXbox

0 -1 -2 = the log2 probabilities of the best trails for 1 round (2^0), 2 rounds (2^-1) and 3 rounds (2^-2)

-- The following command searches for the best differential trails for 5 rounds

./diffsearch 5 31 17 0 24 24 17 31 16 0 -1 -2 -6

- Note:

By default the program finds a single trail with optimal probability. To find all trails, set the macro ALL_TRAILS in defs.h to 1.

- IMPORTANT! information on bounds and offsets (cf. Table 3.2 of the specification)

Denote by

B^r_i: best probability for i rounds, starting from round r: 1 <= r <= 4 (the ARXBox has 4 rounds)

The specification refers to "r" as an "offset". As i increases from 1 to 2,3,4,5,6,7,..., r cycles 1,2,3,4,1,2,3,4,1,...

The bounds B^1_i, B^2_i, B^3_i, B^4_i are given respectively in the four rows of Table 3.2 of the specification.

To find the best probability B^1_6 for 6 rounds (starting from round 1) one has to give as input to the algorithm the following bounds for rounds 1,2,3,4,5 starting at different offsets r:

B^2_1, B^1_2, B^4_3, B^3_4, B^2_5

B^2_1: bound for 1 rounds starting from round 2
B^1_2: bound for 2 rounds starting from round 1
B^4_3: bound for 3 rounds starting from round 4
B^3_4: bound for 4 rounds starting from round 3
B^2_5: bound for 5 rounds starting from round 2

The above bounds stem from the following sequence through which the diffsearch algorithms proceeds:

B^1_6
<= p1 B^2_5
<= p1 p2 B^3_4
<= p1 p2 p3 B^4_3
<= p1 p2 p3 p4 B^1_2
<= p1 p2 p3 p4 p5 B^2_1
<= p1 p2 p3 p4 p5 p6

From Table 3.2 we get the following values for the above bounds

B^2_1 =  0
B^1_2 =  1
B^4_3 =  2
B^3_4 =  6
B^2_5 = 10

Therefore the comamand line parameters for 6 rounds (starting from round 1) would be

./diffsearch 6 31 17 0 24 24 17 31 16 0 -1 -2 -6 -10

For 7 rounds we'll have the following sequence

B^1_7
<= p1 B^2_6
<= p1 p2 B^3_5
<= p1 p2 p3 B^4_4
<= p1 p2 p3 p4 B^1_3
<= p1 p2 p3 p4 p5 B^2_2
<= p1 p2 p3 p4 p5 p6 B^3_1
<= p1 p2 p3 p4 p5 p6 p7

From Table 3.2 we get the following values for the above bounds

B^2_6 = 17
B^3_5 = 10 
B^4_4 =  6
B^1_3 =  2
B^2_2 =  1
B^3_1 =  0

The comamand line parameters for 7 rounds (starting from round 1) would be

./diffsearch 7 31 17 0 24 24 17 31 16 0 -1 -2 -6 -10 -17

Notice how r starts always at 2 at round n-1 (n=6 above) and cycles backwards to 3,4,1,2,3,4,1,... In general, the command line parameters for n rounds, starting from round 1 would be

./diffsearch n 31 17 0 24 24 17 31 16 B^{?}_{1} ... B^{1}_{n-4} B^{4}_{n-3} B^{3}_{n-2} B^{2}_{n-1}

where ? is the value the results from r cycling from 2 at round n-1 to ? at round 1.





