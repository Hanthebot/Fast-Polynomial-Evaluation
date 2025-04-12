#ifndef FFT_MULTIVAR_H
#define FFT_MULTIVAR_H

#include "../galois/gfelement.h"
#include "../nd_vector/nd_vector.h"
#include <map>
#include <span>
#include <chrono>

using namespace shk_galoiscpp;

// Fint: arbitrary size non-negative integer
// F: Field element of finite field of MODULUS

typedef unsigned int u32;
typedef unsigned long long u64;
typedef GFelement F;
typedef GaloisField Field;

// changed
int evaluate_all_point(const nd_vector<Fint>& arr, const Fint& modulo, const Fint& capital_M, Fint& mul_counter, const vector<vector<Fint>>& points_to_evaluate, vector<Fint>& results);
Fint evaluate_brutal(const nd_vector<Fint>& coeff, const vector<Fint>& points, Fint& mul_counter);

const bool DEBUG = false;

template class nd_vector<F>;

#endif // FFT_MULTIVAR_H