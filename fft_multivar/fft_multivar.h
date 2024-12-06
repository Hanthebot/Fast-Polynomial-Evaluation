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
void fft(const vector<F>& w, const span<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter, F& u, F& v);
void fft_multivar_recur(const vector<F>& w, const nd_vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);
void fft_multivar_wrapper(const vector<F>& w, nd_vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter, const F& zero_F);
void fft_nd(const vector<F>& w, const nd_vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);
void compute_manual(const nd_vector<F>& coeff, nd_vector<F>& ans, Fint& mul_counter);

void print_dlog(const nd_vector<F>& coeff, map<F, u32>& dlog, const F& zero_F, const string& prefix);
void check_result(const nd_vector<F>& coeff, const nd_vector<F>& arr, map<F, u32>& dlog, const F& zero_F, Fint& total, Fint& incorrect, const string& prefix = "");

bool getRootOfUnity(const Fint& modulus, const u32& len, Fint& prim_root);

// changed
Field init_setup(const F& zero_F, u32& logn, u32& len, size_t& m);

void rou_init(vector<F>& w, map<F, u32>& dlog, const F& zero_F, u32 logn);

// change
void coeff_init(const nd_vector<F>& coeff, const size_t& m, vector<u32>& degs, const F& zero_F, const vector<u32>& rev, u32 logn);

void rev_init(vector<u32>& rev, u32 logn);

void print_stats(Fint& mul_counter, std::chrono::microseconds duration);

template class nd_vector<F>;

#endif // FFT_MULTIVAR_H