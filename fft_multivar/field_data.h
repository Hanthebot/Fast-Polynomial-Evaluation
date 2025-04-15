#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include "fft_multivar.h"
using namespace std;

struct FieldData {
    u32 prime;                   // prime
    Field field;                 // field element
    vector<F> w;                 // root of unity
    vector<F> dlog;              // discrete logarithm
    vector<u32> rev;             // reverse order of coefficients
    vector<vector<u32>> rev_rev; // inverse of rev, for multiple depth
    vector<u32> radix_vec;           // factorization
    vector<u32> dist;
    F zero_F;                    // zero element of field

    FieldData(): field(2, 1, 1) {}; // default constructor
    FieldData(Field& field, const vector<F>& w, const vector<F>& dlog, const vector<u32>& rev, const vector<vector<u32>>& rev_rev, const vector<u32>& radix_vec, const vector<u32>& dist, u32 prime) 
        : prime(prime), field(field), w(w), dlog(dlog), rev(rev), rev_rev(rev_rev), radix_vec(radix_vec), dist(dist), zero_F(&field, 0) {}
};

#endif // FIELD_DATA_H