#include "galois_own/gfelement.h"
#include <vector>
#include <map>
using namespace shk_galoiscpp;

// Fint: arbitrary size non-negative integer
// F: Field element of finite field of MODULUS

typedef unsigned int u32;
typedef GFelement F;
typedef GaloisField Field;

u32 find_ind(u32 w_exp, u32 modulus);
void fft(vector<F>& w, vector<F> &arr, map<F, u32>& dlog, vector<u32>& rev, u32 logn, Field& field, Fint& mul_counter);