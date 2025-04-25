#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include "fft_multivar.h"
using namespace std;

struct FieldData {
    u32 prime;                   // prime
    Fint prime_mp;               // prime_mp
    Fint rou;                    // root of unity
    vector<Fint> w;              // root of unity power table
    vector<Fint> dlog;           // discrete logarithm
    vector<u32> rev;             // reverse order of coefficients
    vector<vector<u32>> rev_rev; // inverse of rev, for multiple depth
    vector<u32> radix_vec;       // factorization
    vector<u32> dist;

    FieldData() {}; // default constructor
    FieldData(const vector<Fint>& w, const vector<Fint>& dlog, const vector<u32>& rev, const vector<vector<u32>>& rev_rev, const vector<u32>& radix_vec, const vector<u32>& dist, u32 prime) 
        : prime(prime), prime_mp(prime), w(w), dlog(dlog), rev(rev), rev_rev(rev_rev), radix_vec(radix_vec), dist(dist) {
        }
};

#endif // FIELD_DATA_H