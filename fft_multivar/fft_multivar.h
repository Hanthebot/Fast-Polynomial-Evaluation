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

void print_dlog(const nd_vector<F>& coeff, F* dlog, const F& zero_F, const u32& prime, const string& prefix);

bool getRootOfUnity(const Fint& modulus, const u64& len, const u32& prime, Fint& prim_root);
u64 pow(u64 base, u32 exp);
Fint pow(Fint base, Fint exp);

void init_setup(Fint& modulo, vector<u32>& degs_vec, u64& total_len, u32& m);
Field* field_setup(const u32& m, const u32& prime);
void init_evaluation_points(vector<vector<Fint>>& evaluation_points, const u32& m, const Fint& modulo);
void rou_init(F*& w, F*& dlog, const F& zero_F, const u32& prime);
void coeff_init(const nd_vector<Fint>& coeff, const u32& m, const vector<u32>& degs_vec, const Fint& modulo);
void rev_init(u32*& rev, u32 logn);

void print_stats(Fint& mul_counter, std::chrono::microseconds duration);
void compute_M(Fint& capital_M, const vector<u32>& degs_vec, const Fint& modulo, const u64& total_len, const u32& m);

u64 mpz2ull(Fint& mp);
void ull2mpz(Fint& mp, const u64& input);

const bool DEBUG = false;

template class nd_vector<F>;

#endif // FFT_MULTIVAR_H