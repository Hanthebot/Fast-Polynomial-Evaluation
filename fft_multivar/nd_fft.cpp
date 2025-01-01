#include "fft_multivar.h"

using namespace std;

// replaced by example's rev construction
// u32 find_ind(u32 w_exp, u32 modulus) {
//     u32 ind = 0;
//     for (u32 i = (modulus >> 1); w_exp && i > 0; i >>= 1) {
//         if (w_exp % 2) ind |= i;
//         w_exp >>= 1;
//     }
//     return ind;
// }

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const map<F, u32>& dlog, const u32* rev, u32 logn, Fint& mul_counter, const F& zero_F) {
    u32 len = 1ULL << logn;
    for (const auto& s : arr.getShape()) {
        assert(s == len && "size does not match");
    }
    F* temp_u_vec = new F[arr[0].size()];
    F* temp_v_vec = new F[arr[0].size()];
    span<size_t> shape = arr.getShape().subspan(1);
    span<size_t> units = arr.getUnit().subspan(1);
    nd_vector<F> temp_u {arr.getDim() - 1, shape, units, temp_u_vec, arr[0].size()};
    nd_vector<F> temp_v {arr.getDim() - 1, shape, units, temp_v_vec, arr[0].size()};
    fft_multivar_recur(w, arr, dlog, rev, logn, mul_counter, temp_u, temp_v);
}

void fft_multivar_recur(const F* w, const nd_vector<F>& arr, const map<F, u32>& dlog, 
    const u32* rev, u32 logn, Fint& mul_counter,
    nd_vector<F>& temp_u, nd_vector<F>& temp_v) {
    assert(isEqual(arr.getShape().subspan(1), temp_u.getShape()) && "shape does not match");
    assert(isEqual(arr.getShape().subspan(1), temp_v.getShape()) && "shape does not match");
    if (arr.getDim() == 1) {
        F u = temp_u.get();
        F v = temp_v.get();
        fft(w, arr.ptr(), dlog, rev, logn, mul_counter, u, v);
        return;
    }
    // for temporary storage of (n-1)d array
    nd_vector<F> temp_u_ = temp_u[0];
    nd_vector<F> temp_v_ = temp_v[0];
    for (size_t i = 0; i < arr.getShape()[0]; ++i) {
        fft_multivar_recur(w, arr[i], dlog, rev, logn, mul_counter, temp_u_, temp_v_);
    }
    fft_nd(w, arr, dlog, rev, logn, mul_counter, temp_u, temp_v);
}

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft(const F* w, F* const& arr, const map<F, u32>& dlog, const u32* rev, u32 logn, Fint& mul_counter,
    F& u, F& v) {
    u32 len = 1ULL << logn;
    u32 j_rev;
    for (u32 i = len; i > 1; i >>= 1) { // i: gap between next coefficient 
        for (u32 j = 0; j < len; j += i) { // j: margin between initial point and point of interest
            j_rev = rev[j/i];
            for (u32 k = 0; k < (i >> 1); ++k) { // k: initial starting point
                u = arr[j + k]; // P_e(X^2)
                v = arr[j + k + (i >> 1)] * w[j_rev >> 1]; // x*P_o(X^2)
                arr[j + k] = u + v; // P(X) = P_e(X^2) + x*P_o(X^2)
                arr[j + k + (i >> 1)] = u - v; // P(-X) = P_e(X^2) - x*P_0(X^2)
                mul_counter++;
            }
        }
    }
    for (u32 i = 0; i < len; ++i) {
        if (i < rev[i]) {
            u = arr[i];
            arr[i] = arr[rev[i]];
            arr[rev[i]] = u;
        }
    }
}


/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft_nd(const F* w, const nd_vector<F>& arr, const map<F, u32>& dlog, const u32* rev, u32 logn, Fint& mul_counter,
    nd_vector<F>& temp_u, nd_vector<F>& temp_v) {
    u32 len = 1ULL << logn;
    u32 j_rev;
    for (u32 i = len; i > 1; i >>= 1) {
        for (u32 j = 0; j < len; j += i) {
            j_rev = rev[j/i];
            for (u32 k = 0; k < (i >> 1); ++k) {
                temp_u.set_range(arr[j + k]);
                temp_v.set_mul(arr[j + k + (i >> 1)], w[j_rev >> 1]);
                arr[j + k].set_add(temp_u, temp_v);
                arr[j + k + (i >> 1)].set_sub(temp_u, temp_v);
                mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), arr[0].size());
            }
        }
    }
    for (u32 i = 0; i < len; ++i) {
        if (i < rev[i]) {
            arr[i].swap(arr[rev[i]]);
        }
    }
}