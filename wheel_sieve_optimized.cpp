/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard
   arguments: N [-p]
      N: finds primes up to N
     -p: (optional) print the primes found
   optimized single-threaded implementation using a bitset compressed with wheel 3
   2 <= N <= 2000000000
   (like the classic Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements) */

// g++ -O3 -o wheel_sieve_optimized wheel_sieve_optimized.cpp -march=native -lm

#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include <ctime>
#include <math.h>
#include <x86intrin.h>

// constants for wheel 3 compression
const int8_t mod30_to_bit8[] = {-1,0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7};
const uint32_t bit64toval240[] = {
    1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89,91,97,101,103,107,109,113,119,121,127,
    131,133,137,139,143,149,151,157,161,163,167,169,173,179,181,187,191,193,197,199,203,209,211,217,221,223,227,229,233,239};
const char diff[] = {6,4,2,4,2,4,6,2};

#define marked(x,bitset) (bitset[x / 30] & 1 << ((x % 30) * 8 / 30))
#define unmark(x,bitset) bitset[x / 30] &= ~(1 << ((x % 30) * 8 / 30))

uint32_t print(uint64_t *bitmap,  uint32_t bitmapsize) {
    // prints the set bits in W_3-compressed bitmap of bitmapsize 64-bit words, one per line, and returns count
    uint32_t pos = 0, base = 0;
    for (uint32_t k = 0; k < bitmapsize; ++k) {
        uint64_t bitset = bitmap[k];
        while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            uint32_t r = __builtin_ctzll(bitset);
            uint32_t x = base + bit64toval240[r];
            printf("%d\n", x); pos++;
            bitset ^= t;
        }
        base += 240;
    }
    return(pos);
}

uint32_t count(char *bitmap, uint32_t size) {
    // returns number of set bits in bitmap of size bytes
    uint32_t sum = 0, i = size/8;
    uint64_t *p = (uint64_t *) bitmap;
    while (i--) sum += __builtin_popcountll(*p++);
    return(sum);
}

void Extend (char* &bitmap, uint32_t &length, uint32_t n) {
    // Rolls full wheel with given length in W_3 compressed bitmap up to n, and sets length=n
    uint32_t offset = 0;
    for (uint32_t i=1; i < n/length; i++) {
        offset += length/30;
        memcpy(bitmap+offset, bitmap, length/30); // full copy
    }
    offset += length/30;
    uint32_t rem = n % length;
    memcpy(bitmap+offset, bitmap, rem/30); // partial copy
    offset += rem/30;
    uint32_t rem30 = rem % 30;
    if (rem30 > 0) bitmap[offset] = bitmap[rem/30] & (0xFF >> (7-mod30_to_bit8[rem30]));
    length = n;
}

void Delete (char* &bitmap, uint32_t p, uint32_t length) {
    // Deletes multiples of p in W_3 compressed bitmap that are <= length
    /*uint32_t pr240[64];*/
    uint32_t pr240On30[64];
    uint32_t pr240bit8[64];
    for (char r=0; r < 64; r++) { uint32_t t = p*bit64toval240[r];/* pr240[r] = t;*/ pr240On30[r] = t/30; pr240bit8[r] = t%30*8/30; }
    uint64_t *bitmap64 = (uint64_t *)bitmap;
    uint32_t lengthOnp = length/p; if (lengthOnp % 2 == 0) lengthOnp -= 1;
    /*uint32_t p240 = p*240;*/
    uint32_t kmin = p/240, kmax = lengthOnp/240;
    uint32_t bit64 = lengthOnp%240/30*8 + mod30_to_bit8[lengthOnp%240%30];
    /*uint32_t base = kmin*p240;*/
    uint32_t baseon30 = kmin*p*8;
    uint32_t lengthto1on3 = cbrt(length);
    if (p > lengthto1on3) { // no need to stack composites c
        for (uint32_t k=kmin; k < kmax; k++) {
            uint64_t bitset = bitmap64[k];
            while (bitset != 0) {
                uint64_t t = bitset & -bitset;
                uint32_t r = __builtin_ctzll(bitset);
                /*uint32_t c = base + pr240[r];*/
                uint32_t cOn30 = baseon30 + pr240On30[r];
                uint32_t cbit8 = pr240bit8[r];
                bitmap[cOn30] &= ~(0x1 << cbit8); // unmark(c,bitmap)
                bitset ^= t;
            }
            baseon30 += p*8;/* base += p240;*/
        }
        uint64_t bitset = bitmap64[kmax] & (0xFFFFFFFFFFFFFFFF >> (63-bit64));
        while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            uint32_t r = __builtin_ctzll(bitset);
            /*uint32_t c = base + pr240[r];*/
            uint32_t cOn30 = baseon30 + pr240On30[r];
            uint32_t cbit8 = pr240bit8[r];
            bitmap[cOn30] &= ~(0x1 << cbit8); // unmark(c,bitmap)
            bitset ^= t;
        }
    } else { // p <= lengthto1on3
        uint32_t lengthOnp2 = lengthOnp/p; if (lengthOnp2 % 2 == 0) lengthOnp2 -= 1;
        uint32_t kmid = lengthOnp2/240;
        uint32_t bit64mid = lengthOnp2%240/30*8 + mod30_to_bit8[lengthOnp2%240%30];
        uint64_t bitset;
        // stack composites c <= lengthOnp2
        uint32_t* cstack = new uint32_t[((lengthOnp-1)/30+1)*8]; uint32_t ics = 0;
        for (uint32_t k=kmin; k < kmid; k++) {
            bitset = bitmap64[k];
            while (bitset != 0) {
                uint64_t t = bitset & -bitset;
                uint32_t r = __builtin_ctzll(bitset);
                /*uint32_t c = base + pr240[r];*/
                uint32_t cOn30 = baseon30 + pr240On30[r];
                uint32_t cbit8 = pr240bit8[r];
                cstack[ics++] = (cOn30 << 3) | cbit8; // packed
                bitset ^= t;
            }
            baseon30 += p*8;/* base += p240;*/
        }
        bitset = bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF >> (63-bit64mid));
        while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            uint32_t r = __builtin_ctzll(bitset);
            /*uint32_t c = base + pr240[r];*/
            uint32_t cOn30 = baseon30 + pr240On30[r];
            uint32_t cbit8 = pr240bit8[r];
            cstack[ics++] = (cOn30 << 3) | cbit8; // packed
            bitset ^= t;
        }
        // process composites c > lengthOnp2
        bitset = (bit64mid == 63 ? 0 : bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF << bit64mid+1));
        if (kmax > kmid) {
            while (bitset != 0) {
                uint64_t t = bitset & -bitset;
                uint32_t r = __builtin_ctzll(bitset);
                /*uint32_t c = base + pr240[r];*/
                uint32_t cOn30 = baseon30 + pr240On30[r];
                uint32_t cbit8 = pr240bit8[r];
                bitmap[cOn30] &= ~(0x1 << cbit8); // unmark(c,bitmap)  
                bitset ^= t;
            }
            baseon30 += p*8;/* base += p240;*/
            for (uint32_t k=kmid+1; k < kmax; k++) {
                uint64_t bitset = bitmap64[k];
                while (bitset != 0) {
                    uint64_t t = bitset & -bitset;
                    uint32_t r = __builtin_ctzll(bitset);
                    /*uint32_t c = base + pr240[r];*/
                    uint32_t cOn30 = baseon30 + pr240On30[r];
                    uint32_t cbit8 = pr240bit8[r];
                    bitmap[cOn30] &= ~(0x1 << cbit8); // unmark(c,bitmap)  
                    bitset ^= t;
                }
                baseon30 += p*8;/* base += p240;*/
            }
            bitset = bitmap64[kmax] & (0xFFFFFFFFFFFFFFFF >> (63-bit64));
        } else bitset &= (bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF >> (63-bit64)));
        while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            uint32_t r = __builtin_ctzll(bitset);
            /*uint32_t c = base + pr240[r];*/
            uint32_t cOn30 = baseon30 + pr240On30[r];
            uint32_t cbit8 = pr240bit8[r];
            bitmap[cOn30] &= ~(0x1 << cbit8); // unmark(c,bitmap)  
            bitset ^= t;
        }
        // process c <= lengthOnp2 in reverse order
        while (ics > 0) { uint32_t t = cstack[--ics]; bitmap[t >> 3] &= ~(1 << (t & 0x7)); }
        delete[] cstack;
    }
    unmark(p,bitmap);
}

uint32_t Sift(uint32_t N, bool printPrimes) { 
    // counts the primes up to N, printing them if printPrimes
    uint32_t bitmapsize = ((N/30+1 - 1)/8+1)*8; // round to 8 bytes
    char *bitmap = (char*)calloc(bitmapsize, 1);
    /* representation invariant (for the main loop): bitmap of size bitmapsize is the ordered set W compressed with W_3 */
    uint32_t k = 1; if (printPrimes) printf("%d\n", 2);
    if (N >= 3) { k++; printPrimes && printf("%d\n", 3); }
    if (N >= 5) { k++; printPrimes && printf("%d\n", 5); }
    if (N < 7) return(k);
    /* W,k,length = {1,7,11,13,17,19,23,29},4,30: */
    bitmap[0] = 0xFF;
    uint32_t length = 30, p = 7, p_index = 1, p2 = 49;
    if (N < 30) { bitmap[0] &= (0xFF >> (7-mod30_to_bit8[N])); length = N; }
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr(inted) = the first k primes */
    /* (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p2 <= N) {
        if (length < N) { Extend (bitmap, length, std::min(p*length,N)); if (length == N) unmark(1,bitmap); }
        Delete(bitmap, p, length); // (deletes p)
        k++; printPrimes && printf("%d\n", p);
        while (1) { p += diff[p_index++ % 8]; if marked(p,bitmap) break; } /* p = next(W, 1): */
        p2 = p*p;
    }
    if (length < N) Extend (bitmap, length, N);
    unmark(1,bitmap);
    uint32_t nr;
    if (printPrimes) { // print remaining primes: */
        nr = print((uint64_t *)bitmap, bitmapsize/8);
    } else {
        nr = count(bitmap, bitmapsize);
    }
    free(bitmap);
    return(k+nr);
}

int main (int argc, char *argw[]) {
    bool error = false; bool printPrimes = false; uint32_t max = 2000000000;
    uint32_t N, piN;
    if (argc == 3) {
        if (strcmp(argw[2], "-p") == 0) {
            printPrimes = true;
            argc--;
        } else {
            error = true;
        }
    }
    if (argc == 2) {
        N = atoi(argw[1]);
        if (N < 2 || N > max) error = true;
    } else {
        error = true;
    }
    if (error) {
        printf("call with: %s N -p where 2 <= N <= %d and -p to print the primes is optional \n", argw[0], max);
        exit(1);
    }
    int start_s = clock();
    piN = Sift(N, printPrimes);
    int stop_s=clock();
    float duration = (stop_s-start_s)*1E3/double(CLOCKS_PER_SEC);
    printf("%d primes up to %d found in %.2f ms\n", piN, N, duration);
}