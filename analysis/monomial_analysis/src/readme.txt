This is to perform checks on the distribution of low-weight monomials in (step-reduced) versions of Sparkle. In each of the files "degree_sparkleN.cpp", you should adjust DEGREE and SAMPLES as you need. In case DEGREE is > 2, random monomials will be sampled and their occurence in each output bits of the ANF are checked. In case DEGREE is <= 2, all monomials will be checked and no sampling will be done.

The program will output a file "sparkleN_monomials_DEGREE_STEPS.dat" for the solution. The vector below each monomial denotes whether the monomial occurs in the ANF of the corresponding output bits. In particular, the i-th entry is set to 1 if the i-th output bit contains this monomial, and set to 0 otherwise.

At the end of each file it is shown, for each index, how many of the vectors contain a 1, together with the minimum and maxium of the number of occurences over all indices.  

This program thus allows to check whether there are output bits that contain less monomials of a particular degree than expected (expected value is SAMPLES/2).
