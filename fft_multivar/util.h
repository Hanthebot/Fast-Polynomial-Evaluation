#ifndef UTIL_H
#define UTIL_H

#include "fft_multivar.h"

// changed
void init_setup(Fint& modulo, vector<u32>& degs_vec, u64& total_len, u32& m);
void init_evaluation_points(vector<vector<Fint>>& evaluation_points, const u32& m, const Fint& modulo);
void rou_init(vector<Fint>& w, vector<Fint>& dlog, const u32& prime, const Fint& prime_mp, const Fint& rou);
void coeff_init(const nd_vector<Fint>& coeff, const u32& m, const vector<u32>& degs_vec, const Fint& modulo);
void print_stats(Fint& mul_counter, std::chrono::microseconds duration);
void print_dlog(const nd_vector<Fint>& coeff, Fint* dlog, const u32& prime, const string& prefix);

void factorize(Fint n, vector<u32>& factors);
bool getRootOfUnity(const Fint& modulus, const u32& prime, Fint& prim_root);
u64 pow(u64 base, u32 exp);
Fint pow(Fint base, Fint exp);
void compute_M(Fint& capital_M, const vector<u32>& degs_vec, const Fint& modulo, const u64& total_len, const u32& m);
bool is_Fermat_prime(const Fint& prime);
u32 log2(u64 len);
u32 mpz2ul(const Fint& mp);
u64 mpz2ull(const Fint& mp);
void ul2mpz(Fint& mp, const u32& input);
void ull2mpz(Fint& mp, const u64& input);
void enforce_modulus(Fint& mp, const Fint& modulo);
Fint multInverse(const Fint& m, const Fint& prime_mp, const Fint* dlog, const Fint* w);

#endif // UTIL_H