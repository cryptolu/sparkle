LONG TRAIL ARGUMENT
===================


# Usage

This program uses `cmake` to generate the relevant makefile. To
compile it, uses the following command from this folder:

    mkdir bin
    cd bin/
    cmake ..
    make
    
Then, to run the trail search, go in the `bin` folder you just created
and run
    
    ./Trail-Search <type> <min steps> <max steps>

where
- `<type>` must be either 'permutation' or 'absorption'
- `<min steps>` is the minimum number of steps considered
- `<max steps>` is the maximum number of steps considered


# Options

The parameters from the trail search are specified as constants in the
`trail-search.hpp` file as follows:

* `B` is the number of branches on each side, so that the SPARKLE
  instance considered has an internal state of 128*B bits.
  
* `R_LEFT` and `R_RIGHT` are the size of the rate (measured in
  branches) on the left and on the right side. For example, for the
  hash function, use `R_LEFT=2` and `R_RIGHT=0`. For c=128 and r=256
  in SPARKLE384, use `B=R_LEFT=3` and `R_RIGHT=1`.

* `ROTATION_AMOUNT` describes the rotation that is performed just
  after the XOR of the output of M_B in the linear layer of the
  SPARKLE permutation. Its value in SPARKLE is always 1 but feel free
  to modify it to study its effect.

* `DIFFERENTIAL` decides if the search is for differential
  (`DIFFERENTIAL=1`) or for linear (`DIFFERENTIAL=0`) trails. This
  changes both the initial value of the probabilities and the
  definition of the linear layer.
  
To speed up the search, you can uncomment the probability estimates at
the top of the `trail-search.cpp` file. That will allow a proper
Matsui search where more branches are cut in the tree search. In
practice, it is not really needed unless you look at `B=4` and many
steps.

* `USE_MATSUI` defines whether the search algorithm will use knowledge
  of the bounds for fewer steps to speed up the search over a given
  number of steps. If set to 1, it uses the bounds defined at the top
  of `trail-search.cpp` in the ESTIMATE vector.

* `PRINT_TRAILS` decides whether the actual truncated trail are
  printed or not once the one yielding the best differential (or
  linear) bound is found.

* `RATE_TO_RATE` defines whether the pattern in the output of the
  permutation is constrained so as to fit in the rate or not. Does
  nothing when the program is looking for trails in the "permutation"
  case; it is only relevant in the "absorption" case.
