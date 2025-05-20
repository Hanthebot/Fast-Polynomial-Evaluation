#ifndef CRT_DATA_H
#define CRT_DATA_H

#include "fft_multivar.h"
#include "field_data.h"
#include "nd_fft.h"
#include "util.h"
using namespace std;

struct CRTData {
    u32 modulo;                   // prime
    bool isPrime;
    Fint prime_mp;               // prime_mp
    vector<u32> fields_used;
    CRTData** subField = nullptr; // subfield data (array)
    FieldData* FD = nullptr; // subfield data (array)
    Fint M;
    Fint buffer_size = 0;
    vector<u32> interpolate;
    vector<size_t> shape_vec;
    vector<size_t> unit_vec;
    Fint common_prod;
    nd_vector<Fint>* coeff = nullptr; // pointer to coefficient data
    nd_vector<Fint>* preprocessed = nullptr;

    CRTData() = default; // default constructor

    // Construct for a single field (leaf node)
    CRTData(nd_vector<Fint>& arr, u32 prime);

    // Recursively count tables
    u32 getTableCount() const;

    // Recursively count table size
    Fint getTableSize() const;

    // Direct access for leaf node
    Fint access(const vector<Fint>& eval_point);

    // Recursive evaluation
    Fint evaluate(const vector<Fint>& eval_point);
    void evaluate_stats(Fint& mul_count, Fint& access_count);

    void checkSpaceOptimum(const span<size_t>& shape);
    void findFieldsCRT();
    void preprocess();
    void printTableSize(string st = "") const;
};

#endif // CRT_DATA_H