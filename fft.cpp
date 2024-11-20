#include "fft_finite.h"

using namespace std;

// replaced by example's rev construction
// u32 find_ind(u32 w_exp, u32 modulus) {
//     u32 ind = 0;
//     for (u32 i = (modulus >> 1); w_exp && i > 0; i >>= 1) {
//         if (w_exp % 2) ind |= i;
//         w_exp >>= 1;
//     }
//     return ind;
// }

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft(const vector<F>& w, vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Field& field, Fint& mul_counter) {
    u32 len = 1ULL << logn;
    assert(arr.size() == len);

    // compute all values of w, manually
    
    F u,v;
    u.setField(&field);
    v.setField(&field);
    Fint m = field.getModulus();

    u32 j_rev;
    for (u32 i = len; i > 1; i >>= 1) { // i: gap between next coefficient 
        for (u32 j = 0; j < len; j += i) { // j: margin between initial point and point of interest
            j_rev = rev[j/i];
            for (u32 k = 0; k < (i >> 1); ++k) { // k: initial starting point
                u = arr[j + k]; // P_e(X^2)
                v = arr[j + k + (i >> 1)] * w[j_rev >> 1]; // x*P_o(X^2)
                arr[j + k] = u + v; // P(X) = P_e(X^2) + x*P_o(X^2)
                arr[j + k + (i >> 1)] = u - v; // P(-X) = P_e(X^2) - x*P_0(X^2)
                mul_counter++;
            }
        }
    }
    for (u32 i = 0; i < len; ++i) {
        if (i < rev[i]) {
            u = arr[i];
            arr[i] = arr[rev[i]];
            arr[rev[i]] = u;
        }
    }
}