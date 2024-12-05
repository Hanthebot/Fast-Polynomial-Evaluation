#include "galois/gfelement.h"
#include <vector>
#include <map>
#include <chrono>

using namespace shk_galoiscpp;

// Fint: arbitrary size non-negative integer
// F: Field element of finite field of MODULUS

typedef unsigned int u32;
typedef GFelement F;
typedef GaloisField Field;

u32 find_ind(u32 w_exp, u32 modulus);
void fft(const vector<F>& w, vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Field& field, Fint& mul_counter);

bool getRootOfUnity(const Fint& modulus, const u32& len, Fint& prim_root);

Field init_setup(const F& zero_F, u32& logn, u32& len);

void rou_init(vector<F>& w, map<F, u32>& dlog, const F& zero_F, u32 logn);

void coeff_init(vector<F>& coeff, u32& deg, const F& zero_F, const vector<u32>& rev, u32 logn);

void rev_init(vector<u32>& rev, u32 logn);

void print_stats(Fint& mul_counter, std::chrono::microseconds duration);