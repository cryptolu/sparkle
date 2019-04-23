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
