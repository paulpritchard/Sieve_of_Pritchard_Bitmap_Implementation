/* Sieve of Pritchard in C++, as described at https://en.wikipedia.org/wiki/Sieve_of_Pritchard
   arguments: N [-p]
      N: finds primes up to N
     -p: (optional) print the primes found
   optimized single-threaded implementation using a bitset compressed with wheel 3
   2 <= N <= 4000000000
   (like the classic Sieve of Eratosthenes, this algorithm is not suitable for very large N due to memory requirements) */

// g++ -O3 -o wheel_sieve_bitmap wheel_sieve_bitmap.cpp -march=native -lm

#include <cstring>
#include <stdio.h>
#include <ctime>
#include <math.h>
#include <inttypes.h>

// constants for wheel 3 compression
const int mod30_to_bit8[] = {-1,0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7};
const uint64_t bit64toval240[] = {
    1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89,91,97,101,103,107,109,113,119,121,127,
    131,133,137,139,143,149,151,157,161,163,167,169,173,179,181,187,191,193,197,199,203,209,211,217,221,223,227,229,233,239};
const uint64_t diff[] = {6,4,2,4,2,4,6,2};

#define marked(x,bitset) (bitset[x / 30] & 1 << ((x % 30) * 8 / 30))
#define unmark(x,bitset) bitset[x / 30] &= ~(1 << ((x % 30) * 8 / 30))

uint64_t print(uint64_t *bitmap, uint64_t bitmapsize) {
    // prints the set bits in W_3-compressed bitmap of bitmapsize 64-bit words, one per line, and returns count
    uint64_t pos = 0, base = 0;
    for (uint64_t k = 0; k < bitmapsize; ++k) {
        uint64_t bitset = bitmap[k];
        while (bitset != 0) {
            uint64_t t = bitset & -bitset;
            uint64_t r = __builtin_ctzll(bitset);
            uint64_t x = base + bit64toval240[r];
            printf("%lu\n", x); pos++;
            bitset ^= t;
        }
        base += 240;
    }
    return(pos);
}

uint64_t count(char *bitmap, uint64_t size) {
    // returns number of set bits in bitmap of size bytes
    uint64_t sum = 0, i = size/8;
    uint64_t *p = (uint64_t *) bitmap;
    while (i--) sum += __builtin_popcountll(*p++);
    return(sum);
}

void Extend (char* &bitmap, uint64_t &length, uint64_t n) {
    // Rolls full wheel in bitmap with given length up to n, and sets length=n
    uint64_t offset = 0;
    for (uint64_t i=1; i < n/length; i++) {
        offset += length/30;
        memcpy(bitmap+offset, bitmap, length/30); // full copy
    }
    offset += length/30; uint64_t rem = n % length;
    memcpy(bitmap+offset, bitmap, rem/30); // partial copy
    offset += rem/30; uint64_t rem30 = rem % 30;
    if (rem30 > 0) bitmap[offset] = bitmap[rem/30] & (0xFF >> (7-mod30_to_bit8[rem30]));
    length = n;
}

// Can be done much faster with AVX-512. See:
// https://lemire.me/blog/2022/05/10/faster-bitset-decoding-using-intel-avx-512/
#define doLeastSetBit(bitset,cOn30,cbit8mask) uint64_t t = bitset & -bitset;\
    uint64_t r = __builtin_ctzll(bitset);\
    /*uint64_t c = base + pr240[r];*/\
    uint64_t cOn30 = baseon30 + pr240On30[r];\
    uint8_t cbit8mask = pr240bit8mask[r];\
    bitset ^= t;

void Delete (char* &bitmap, uint64_t p, uint64_t length) {
    // Deletes multiples of p in bitmap that are <= length
    uint32_t /*pr240[64],*/pr240On30[64]; // < p*8
    uint8_t pr240bit8mask[64];
    for (uint32_t r=0; r < 64; r++) {
        uint64_t t = p*bit64toval240[r];/* pr240[r] = t;*/ pr240On30[r] = t/30; pr240bit8mask[r] = ~(0x1 << t%30*8/30);
    }
    uint64_t *bitmap64 = (uint64_t *)bitmap;
    uint64_t maxf = length/p; if (maxf % 2 == 0) maxf -= 1;
    /*uint64_t p240 = p*240;*/
    uint64_t kmin = p/240, kmax = maxf/240;
    uint64_t bit64 = maxf%240/30*8 + mod30_to_bit8[maxf%240%30];
    uint64_t baseon30 = kmin*p*8/*, base = kmin*p240*/;
    uint64_t lengthto1on3 = cbrt(length);
    if (p > lengthto1on3) { // no need to stack composites c
        for (uint64_t k=kmin; k < kmax; k++) {
            uint64_t bitset = bitmap64[k];
            while (bitset != 0) {
                doLeastSetBit(bitset,cOn30,cbit8mask);
                bitmap[cOn30] &= cbit8mask; // unmark(c,bitmap)
            }
            baseon30 += p*8;/* base += p240;*/
        }
        uint64_t bitset = bitmap64[kmax] & (0xFFFFFFFFFFFFFFFF >> (63-bit64));
        while (bitset != 0) {
            doLeastSetBit(bitset,cOn30,cbit8mask);
            bitmap[cOn30] &= cbit8mask; // unmark(c,bitmap)
        }
        unmark(p,bitmap);
        return;
    }
    // p <= lengthto1on3, so both types of factors
    uint64_t maxfOnp = maxf/p; if (maxfOnp % 2 == 0) maxfOnp -= 1;
    uint64_t kmid = maxfOnp/240;
    uint64_t bit64mid = maxfOnp%240/30*8 + mod30_to_bit8[maxfOnp%240%30];
    // stack composites c <= maxfOnp
    uint64_t* cstack = new uint64_t[((maxf-1)/30+1)*8]; uint64_t ics = 0;
    for (uint64_t k=kmin; k < kmid; k++) {
        uint64_t bitset = bitmap64[k];
        while (bitset != 0) {
            doLeastSetBit(bitset,cOn30,cbit8mask);
            cstack[ics++] = (cOn30 << 8) | cbit8mask; // packed
        }
        baseon30 += p*8;/* base += p240;*/
    }
    uint64_t bitset = bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF >> (63-bit64mid));
    while (bitset != 0) {
        doLeastSetBit(bitset,cOn30,cbit8mask);
        cstack[ics++] = (cOn30 << 8) | cbit8mask; // packed
    }
    // process composites c > maxfOnp
    bitset = (bit64mid == 63 ? 0 : bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF << bit64mid+1));
    if (kmax > kmid) {
        while (bitset != 0) {
            doLeastSetBit(bitset,cOn30,cbit8mask);
            bitmap[cOn30] &= cbit8mask; // unmark(c,bitmap)  
        }
        baseon30 += p*8;/* base += p240;*/
        for (uint64_t k=kmid+1; k < kmax; k++) {
            uint64_t bitset = bitmap64[k];
            while (bitset != 0) {
                doLeastSetBit(bitset,cOn30,cbit8mask);
                bitmap[cOn30] &= cbit8mask; // unmark(c,bitmap)  
            }
            baseon30 += p*8;/* base += p240;*/
        }
        bitset = bitmap64[kmax] & (0xFFFFFFFFFFFFFFFF >> (63-bit64));
    } else bitset &= (bitmap64[kmid] & (0xFFFFFFFFFFFFFFFF >> (63-bit64)));
    while (bitset != 0) {
        doLeastSetBit(bitset,cOn30,cbit8mask);
        bitmap[cOn30] &= cbit8mask; // unmark(c,bitmap)  
    }
    // process c <= maxfOnp in reverse order
    while (ics > 0) { uint64_t t = cstack[--ics]; bitmap[t >> 8] &= (t & 0xFF); }
    delete[] cstack;
    unmark(p,bitmap);
}

uint64_t Sift(uint64_t N, bool printPrimes) { 
    // counts the primes up to N, printing them if printPrimes
    uint64_t bitmapsize = ((N/30+1 - 1)/8+1)*8; // round to 8 bytes
    char *bitmap = (char*)calloc(bitmapsize, 1);
    /* representation invariant (for the main loop): bitmap of size bitmapsize is the ordered set W compressed with W_3 */
    uint64_t k = 1; if (printPrimes) printf("%u\n", 2);
    if (N >= 3) { k++; printPrimes && printf("%u\n", 3); }
    if (N >= 5) { k++; printPrimes && printf("%u\n", 5); }
    if (N < 7) return(k);
    /* W,k,length = {1,7,11,13,17,19,23,29},4,30: */
    bitmap[0] = 0xFF;
    uint64_t length = 30, p = 7, p2 = 49;
    uint32_t p_index = 1;
    if (N < 30) { bitmap[0] &= (0xFF >> (7-mod30_to_bit8[N])); length = N; }
    /* invariant: p = p_(k+1) and W = W_k inter {1,...,N} and length = min(P_k,N) and Pr(inted) = the first k primes */
    /* (where p_i denotes the i'th prime, W_i denotes the i'th wheel, P_i denotes the product of the first i primes) */
    while (p2 <= N) {
        if (length < N) { Extend (bitmap, length, std::min(p*length,N)); if (length == N) unmark(1,bitmap); }
        Delete(bitmap, p, length); // (deletes p)
        k++; printPrimes && printf("%lu\n", p);
        while (1) { p += diff[p_index++ % 8]; if marked(p,bitmap) break; } /* p = next(W, 1): */
        p2 = p*p;
    }
    if (length < N) Extend (bitmap, length, N);
    unmark(1,bitmap);
    uint64_t nr = (printPrimes ? print((uint64_t *)bitmap, bitmapsize/8) : count(bitmap, bitmapsize));
    free(bitmap);
    return(k+nr);
}

int main (int argc, char *argw[]) {
    #if __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
    printf("sorry: this code relies on little endian byte order\n"); exit(1);
    #endif
    bool error = false; bool printPrimes = false; uint64_t max = 100000000000;
    uint64_t N, piN;
    if (argc == 3) {
        if (strcmp(argw[2], "-p") == 0) { printPrimes = true; argc--; } else error = true;
    }
    if (argc == 2) { N = strtol(argw[1], NULL, 10); if (N < 2 || N > max) error = true; }
    else error = true;
    if (error) { printf("call with: %s N -p where 2 <= N <= %lu and -p to print the primes is optional \n", argw[0], max); exit(1); }
    int start_s = clock();
    piN = Sift(N, printPrimes);
    float duration = (clock()-start_s)*1E3/double(CLOCKS_PER_SEC);
    printf("%lu primes up to %lu found in %.2f ms\n", piN, N, duration);
}