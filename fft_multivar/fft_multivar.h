#ifndef FFT_MULTIVAR_H
#define FFT_MULTIVAR_H

#include <gmpxx.h>
#include "../nd_vector/nd_vector.h"
#include <map>
#include <vector>
#include <span>
#include <chrono>

using namespace std;

// Fint: arbitrary size non-negative integer
// F: Field element of finite field of MODULUS

typedef unsigned int u32;
typedef unsigned long long u64;
typedef mpz_class Fint;

struct TimeRecord {
    long long int precomp = 0;
    long long int mem_clean = 0;
    long long int fft = 0;
    long long int retrieval = 0;
    long long int crt = 0;
};

struct EvalIO {
    TimeRecord times;
    bool print_time = false;
    bool debug = false;
};

// changed
int evaluate_all_point(const nd_vector<Fint>& arr, const Fint& modulo, const Fint& capital_M, Fint& mul_counter, const vector<vector<Fint>>& points_to_evaluate, vector<Fint>& results, EvalIO& meta);
Fint evaluate_brutal(const nd_vector<Fint>& coeff, const vector<Fint>& points, Fint& mul_counter);

#endif // FFT_MULTIVAR_H