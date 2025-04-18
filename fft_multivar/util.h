#ifndef UTIL_H
#define UTIL_H

#include "fft_multivar.h"

using namespace shk_galoiscpp;

// changed
void init_setup(Fint& modulo, vector<u32>& degs_vec, u64& total_len, u32& m);
Field* field_setup(const u32& m, const u32& prime);
void init_evaluation_points(vector<vector<Fint>>& evaluation_points, const u32& m, const Fint& modulo);
void rou_init(F*& w, F*& dlog, const F& zero_F, const u32& prime);
void coeff_init(const nd_vector<Fint>& coeff, const u32& m, const vector<u32>& degs_vec, const Fint& modulo);
void print_stats(Fint& mul_counter, std::chrono::microseconds duration);
void print_dlog(const nd_vector<F>& coeff, F* dlog, const F& zero_F, const u32& prime, const string& prefix);

void factorize(Fint n, vector<u32>& factors);
bool getRootOfUnity(const Fint& modulus, const u32& prime, Fint& prim_root);
u64 pow(u64 base, u32 exp);
Fint pow(Fint base, Fint exp);
void compute_M(Fint& capital_M, const vector<u32>& degs_vec, const Fint& modulo, const u64& total_len, const u32& m);
bool is_Fermat_prime(const Fint& prime);
u32 log2(u64 len);
u32 mpz2ul(const Fint& mp);
u64 mpz2ull(const Fint& mp);
void ull2mpz(Fint& mp, const u64& input);

#endif // UTIL_H