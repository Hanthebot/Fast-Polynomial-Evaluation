#include "fft_multivar.h"

#include <algorithm>
#include <execution>
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
void fft_multivar_wrapper(const vector<F>& w, nd_vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter, const F& zero_F) {
    u32 len = 1ULL << logn;
    for (const auto& s : arr.getShape()) {
        assert(s == len && "size does not match");
    }
    vector<F> temp_u_vec(arr[0].size(), zero_F);
    vector<F> temp_v_vec(arr[0].size(), zero_F);
    nd_vector<F> temp_u {arr.getDim() - 1, arr.getShape().subspan(1), temp_u_vec};  
    nd_vector<F> temp_v {arr.getDim() - 1, arr.getShape().subspan(1), temp_v_vec};
    u32 m = arr.getDim();
    u32 last_dim = arr.getShape()[m - 1];
    vector<u32> parallel_vec(arr.size() / last_dim);
    iota(parallel_vec.begin(), parallel_vec.end(), 0);
    auto arr_begin = arr.span().begin();
    // first dimension
    mutex mul_counter_mutex;
    for_each(execution::par, 
        parallel_vec.begin(), 
        parallel_vec.end(), 
        [&](const u32& i) {
            F u = temp_u.get(i);
            F v = temp_v.get(i);
            span<F> s {arr_begin + (i * last_dim), arr_begin + (i + 1) * last_dim};
            fft(w, s, dlog, rev, logn, mul_counter, u, v, mul_counter_mutex);
        }
    );
    // other dimensions
    span<size_t> shape {arr.getShape()};
    vector<size_t> original_shape_vec;
    original_shape_vec.resize(m, 1);
    span<size_t> original_shape {original_shape_vec};
    copy(shape.begin(), shape.end(), original_shape.begin());
    for (size_t dim = 1; dim < m; ++dim) { // for each dimension, inner work is parallelized
        vector<size_t> temp_shape_vec;
        temp_shape_vec.resize(dim + 2, 1);
        span<size_t> shape {original_shape.subspan(m - dim - 2)};
        copy(shape.begin(), shape.end(), temp_shape_vec.begin());
        span<size_t> temp_shape {temp_shape_vec};
        span<size_t> temp_shape_sub {temp_shape.begin(), temp_shape.begin() + dim + 1};
        u32 fix_size = 1;
        for (size_t i = 0; i < m - dim - 1; ++i) {
            fix_size *= original_shape[i];
        }
        temp_shape[0] = fix_size;
        arr.reshape(temp_shape);
        temp_u.reshape(temp_shape_sub);
        temp_v.reshape(temp_shape_sub);
        for_each(execution::par, 
            parallel_vec.begin(), 
            parallel_vec.begin() + fix_size, 
            [&](const u32& i) {
                nd_vector<F> temp_u_ = temp_u[i];
                nd_vector<F> temp_v_ = temp_v[i];
                fft_nd(w, arr[i], dlog, rev, logn, mul_counter, temp_u_, temp_v_, mul_counter_mutex);
            }
        );
    }
    arr.reshape(shape);
}

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft(const vector<F>& w, const span<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter,
    F& u, F& v, mutex& mul_counter_mutex) {
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
                {
                    lock_guard<mutex> guard(mul_counter_mutex);
                    mul_counter++;
                }
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
void fft_nd(const vector<F>& w, const nd_vector<F> &arr, const map<F, u32>& dlog, const vector<u32>& rev, u32 logn, Fint& mul_counter,
    nd_vector<F>& temp_u, nd_vector<F>& temp_v, mutex& mul_counter_mutex) {
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
                {
                    lock_guard<mutex> guard(mul_counter_mutex);
                    mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), arr[0].size());
                }
            }
        }
    }
    for (u32 i = 0; i < len; ++i) {
        if (i < rev[i]) {
            arr[i].swap(arr[rev[i]]);
        }
    }
}