# An optimized implementation of the dynamic wheel sieve of Pritchard

This repository concerns the dynamic wheel sieve of Pritchard as described on its [Wikipedia page](https://en.wikipedia.org/wiki/Sieve_of_Pritchard).
It is an algorithm for finding the primes up to a given limit N that takes sublinear time O(N/log log N).
The repository gives an optimized implementation of the algorithm, aiming to minimize the execution time for N around 10^9.

The original implementation is described in the paper
[Paul Pritchard, "A Sublinear Additive Sieve for Finding Prime Numbers", *Communications of the ACM*, vol. 24, no. 1, pp. 18–23](https://dl.acm.org/doi/10.1145/358527.358540).
A detailed code animation with a limit of N=150 is given in [this video](https://www.youtube.com/watch?v=GxgGMwLfTjE).

Since the dynamic wheel sieve was presented as a contribution to computational complexity (and as a beautiful algorithm),
the original implementation was designed to demonstrate the sublinear running time as simply and clearly as possible.
Nevertheless, it is interesting to ask just how fast (or not) this algorithm can be.
The implementation herein is our best attempt at answering that question.

## The approach taken

The classic approach to optimizing the sieve of Eratosthenes is to use a bitset, i.e.,
a sequence of bits representing the sifted status of successive prime candidates:
0 for known to be composite (sifted by a prime already considered), and 1 for possibly prime (not sifted by the primes already considered).
The simplest way is to use 1 bit for each number up to the limit N.
Using 1 bit for each odd number is a common improvement, equivalent to using wheel 1 rolled up to N.
We use a compression scheme that exploits the fact that wheel 3 (for primes 2,3,5) has 8 members, which fit nicely in one byte.
See [this explanation](http://tverniquet.com/hprime/#p1:i2).

The only other important technique used is the intrinsic *__builtin_ctzll* in the g++ compiler,
which enables iterating over set bits quickly (especially on modern x64 processors).
See [here](https://lemire.me/blog/2018/03/08/iterating-over-set-bits-quickly-simd-edition/) for a discussion which also presents
an alternative SIMD approach (i.e. operating on short vectors).
The intrinsic *__builtin_popcountll* is also for the final counting of set bits.

Note that the code does not use explicit SIMD processing (though optimizing compilers may do so),
or multiprocessing, extensive loop unrolling, or extensive precomputation.

### wheel_sieve_optimized.cpp

This is a self-contained C++ program.
It essentially implements the abstract algorithm,
differing only in starting with wheel 3 (and therefore taking special action for N \< 30).

Provided *__builtin_ctzll* finds each successive set bit in O(1) time (which is the case on mopdern processors),
the resulting program runs in O(N/log log N) time, and the implicit constant factor is very small.
It should be competitive with optimized single-threaded non-segmented implementations of the classic sieve of Eratosthenes.