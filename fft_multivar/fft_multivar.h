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
typedef GFelement F;
typedef GaloisField Field;

// changed
void fft(const F* w, F* const& arr, const u32* rev, u32 logn, Fint& mul_counter, F& u, F& v);
void fft_multivar_recur(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);
void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const F* dlog, const u32* rev, u32 logn, Fint& mul_counter, const F& zero_F);
void fft_nd(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);

void print_dlog(const nd_vector<F>& coeff, F* dlog, const F& zero_F, const u32& prime, const string& prefix);

bool getRootOfUnity(const Fint& modulus, const u32& len, const u32& prime, Fint& prim_root);

// changed
Field init_setup(u32& logn, u32& len, size_t& m, u32& prime);

void rou_init(F*& w, F*& dlog, const F& zero_F, u32 prime);

// change
void coeff_init(const nd_vector<F>& coeff, const size_t& m, vector<u32>& degs_vec, u32 prime, Field* const field);

void rev_init(u32*& rev, u32 logn);

void print_stats(Fint& mul_counter, std::chrono::microseconds duration);

template class nd_vector<F>;

#endif // FFT_MULTIVAR_H