#include "nd_fft.h"
#include "util.h"

void fft(const F* w, F* const& arr, const u32* rev, u32 logn, Fint& mul_counter, F& u, F& v);
void fft_nd(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);

void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter,
    F*& temp_u_vec, F*& temp_v_vec) {
    u32 len = 1ULL << logn;
    for (const auto& s : arr.getShape()) {
        assert(s == len && "size does not match");
    }
    span<size_t> shape = arr.getShape().subspan(1);
    span<size_t> units = arr.getUnit().subspan(1);
    nd_vector<F> temp_u {arr.getDim() - 1, shape, units, temp_u_vec, arr[0].size()};
    nd_vector<F> temp_v {arr.getDim() - 1, shape, units, temp_v_vec, arr[0].size()};
    vector<size_t> iter_array_shape;
    vector<size_t> iter_temp_shape;
    size_t whole_size = arr.size();
    size_t shape_default = arr.getShape()[0];
    vector<size_t> original_shape(arr.getDim(), shape_default);
    for (size_t dim_i = 1; dim_i <= original_shape.size(); ++dim_i) {
        size_t iter_size = whole_size / pow(shape_default, dim_i);
        iter_array_shape.clear();
        iter_array_shape.resize(dim_i + 1, shape_default);
        iter_temp_shape.resize(dim_i, shape_default);
        iter_temp_shape[0] = iter_array_shape[0] = iter_size;
        arr.reshape(iter_array_shape);
        temp_u.reshape(iter_temp_shape);
        temp_v.reshape(iter_temp_shape);
        if (dim_i == 1) {
            for (size_t arr_i = 0; arr_i < iter_size; ++arr_i) {
                F& u = temp_u.get(arr_i);
                F& v = temp_v.get(arr_i);
                fft(w, arr[arr_i].ptr(), rev, logn, mul_counter, u, v);
            }
        } else {
            for (size_t arr_i = 0; arr_i < iter_size; ++arr_i) {
                nd_vector<F> u = temp_u[arr_i];
                nd_vector<F> v = temp_v[arr_i];
                fft_nd(w, arr[arr_i], rev, logn, mul_counter, u, v);
            }
        }
    }
    arr.reshape(original_shape);
}

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft(const F* w, F* const& arr, const u32* rev, u32 logn, Fint& mul_counter,
    F& u, F& v) {
    u32 len = 1ULL << logn;
    for (u32 i = len; i > 1; i >>= 1) { // i: gap between next coefficient 
        u32 j_count = 0;
        for (u32 j = 0; j < len; j += i) { // j: margin between initial point and point of interest
            // j_rev = rev[j/i];
            u32 j_rev = rev[j_count];
            for (u32 k = 0; k < (i >> 1); ++k) { // k: initial starting point
                u = arr[j + k]; // P_e(X^2)
                v = arr[j + k + (i >> 1)] * w[j_rev >> 1]; // x*P_o(X^2)
                arr[j + k] = u + v; // P(X) = P_e(X^2) + x*P_o(X^2)
                arr[j + k + (i >> 1)] = u - v; // P(-X) = P_e(X^2) - x*P_0(X^2)
            }
            mul_counter += (i >> 1);
            ++j_count;
        }
    }
}


/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft_nd(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter,
    nd_vector<F>& temp_u, nd_vector<F>& temp_v) {
    u32 len = 1ULL << logn;
    for (u32 i = len; i > 1; i >>= 1) {
        u32 j_count = 0;
        for (u32 j = 0; j < len; j += i) {
            u32 j_rev = rev[j_count];
            for (u32 k = 0; k < (i >> 1); ++k) {
                temp_u.set_range(arr[j + k]);
                temp_v.set_mul(arr[j + k + (i >> 1)], w[j_rev >> 1]);
                arr[j + k].set_add(temp_u, temp_v);
                arr[j + k + (i >> 1)].set_sub(temp_u, temp_v);
            }
            mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), arr[0].size() * (i >> 1));
            ++j_count;
        }
    }
}

void fft_NF(const F* w, F* const& arr, const vector<u32>& radix_vec, const vector<u32>& dist, const vector<u32>* rev_rev, const u32& len, Fint& mul_counter, F* const& temp_u, F& temp_v1, F& temp_v2);
void fft_nd_NF(const F* w, const nd_vector<F>& arr, const vector<u32>& radix_vec, const vector<u32>& dist, const vector<u32>* rev_rev, const u32& len, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v1, nd_vector<F>& temp_v2);

void fft_multivar_wrapper_NF(const F* w, nd_vector<F>& arr, const vector<u32>& radix_vec, const vector<u32>& dist, const vector<u32>* rev_rev, 
    Fint& mul_counter, F*& temp_u_vec, F*& temp_v1_vec, F*& temp_v2_vec) {
    nd_vector<F> temp_u {arr.getDim(), arr.getShape(), arr.getUnit(), temp_u_vec, arr.size()};
    // buffer for lower dimension
    span<size_t> shape = arr.getShape().subspan(1);
    span<size_t> units = arr.getUnit().subspan(1);
    nd_vector<F> temp_v1 {arr.getDim() - 1, shape, units, temp_v1_vec, arr[0].size()};
    nd_vector<F> temp_v2 {arr.getDim() - 1, shape, units, temp_v2_vec, arr[0].size()};
    vector<size_t> iter_array_shape;
    vector<size_t> iter_temp_shape;
    size_t whole_size = arr.size();
    size_t shape_default = arr.getShape()[0];
    u32 len = shape_default; // prime - 1 
    vector<size_t> original_shape(arr.getDim(), shape_default);
    for (size_t dim_i = 1; dim_i <= original_shape.size(); ++dim_i) {
        size_t iter_size = whole_size / pow(shape_default, dim_i);
        iter_array_shape.clear();
        // changing shape while nd_vector is holding it: not safe in general
        iter_array_shape.resize(dim_i + 1, shape_default);
        iter_temp_shape.resize(dim_i, shape_default);
        iter_temp_shape[0] = iter_array_shape[0] = iter_size;
        arr.reshape(iter_array_shape);
        temp_u.reshape(iter_array_shape);
        temp_v1.reshape(iter_temp_shape);
        temp_v2.reshape(iter_temp_shape);
        if (dim_i == 1) {
            for (size_t arr_i = 0; arr_i < iter_size; ++arr_i) {
                F* u = temp_u[arr_i].ptr();
                F& v1 = temp_v1.get(arr_i);
                F& v2 = temp_v2.get(arr_i);
                fft_NF(w, arr[arr_i].ptr(), radix_vec, dist, rev_rev, len, mul_counter, u, v1, v2);
            }
        } else {
            for (size_t arr_i = 0; arr_i < iter_size; ++arr_i) {
                nd_vector<F> u = temp_u[arr_i];
                nd_vector<F> v1 = temp_v1[arr_i];
                nd_vector<F> v2 = temp_v2[arr_i];
                fft_nd_NF(w, arr[arr_i], radix_vec, dist, rev_rev, len, mul_counter, u, v1, v2);
            }
        }
    }
    arr.reshape(original_shape);
}

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft_NF(const F* w, F* const& arr, const vector<u32>& radix_vec, const vector<u32>& dist, const vector<u32>* rev_rev, const u32& len, 
    Fint& mul_counter, F* const& temp_u, F& temp_v1, F& temp_v2) {
    for (size_t phase = 1; phase <= radix_vec.size(); ++phase) {
        u32 bi = 0;
        u32 radix = radix_vec[phase - 1];
        u32 phase_num = len / radix;
        for (u32 block = 0; block < len; block += dist[phase - 1]) {
            u32 rev_dist = dist[phase] * rev_rev[phase - 1][bi];
            for (u32 i = 0; i < dist[phase]; i++) {
                u32 j_dist = 0;
                for (u32 j = 0; j < radix; ++j) {
                    temp_u[j] = arr[(block + i + j_dist) % len];
                    j_dist += dist[phase];
                }
                j_dist = 0;
                u32 exp_diff = rev_dist;
                for (u32 j = 0; j < radix; ++j) {
                    temp_v1.setZero();
                    u32 exp = 0;
                    for (u32 k = 0; k < radix; ++k) {
                        temp_v2 = temp_u[k] * w[exp];
                        temp_v1 = temp_v1 + temp_v2;
                        exp += exp_diff;
                        if (exp >= len) {
                            exp %= len;
                        }
                    }
                    mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), radix);
                    arr[(block + i + j_dist) % len] = temp_v1;
                    j_dist += dist[phase];
                    exp_diff += phase_num;
                }
            }
            ++bi;
            mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), block / dist[phase - 1]);
        }
    }
}


/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
void fft_nd_NF(const F* w, const nd_vector<F>& arr, const vector<u32>& radix_vec, const vector<u32>& dist, const vector<u32>* rev_rev, const u32& len, 
    Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v1, nd_vector<F>& temp_v2) {
    for (u32 phase = 1; phase <= radix_vec.size(); ++phase) {
        u32 bi = 0;
        u32 radix = radix_vec[phase - 1];
        u32 phase_num = len / radix;
        for (u32 block = 0; block < len; block += dist[phase - 1]) {
            u32 rev_dist = rev_rev[phase - 1][bi] * dist[phase];
            for (u32 i = 0; i < dist[phase]; i++) {
                u32 j_dist = 0;
                for (u32 j = 0; j < radix; ++j) {
                    temp_u[j].set_range(arr[(block + i + j_dist) % len]);
                    j_dist += dist[phase];
                }
                j_dist = 0;
                u32 exp_diff = rev_dist;
                for (u32 j = 0; j < radix; ++j) {
                    // reset temp to 0
                    for (u64 k = 0; k < temp_v1.size(); ++k) {
                        temp_v1.get(k).setZero();
                    }
                    u32 exp = 0;
                    for (u32 k = 0; k < radix; ++k) {
                        temp_v2.set_mul(temp_u[k], w[exp]);
                        temp_v1.set_add(temp_v1, temp_v2);
                        exp += exp_diff;
                        if (exp >= len) {
                            exp %= len;
                        }
                    }
                    mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), radix * temp_u[0].size());
                    arr[(block + i + j_dist) % len].set_range(temp_v1);
                    j_dist += dist[phase];
                    exp_diff += phase_num;
                }
            }
            ++bi;
            mpz_add_ui(mul_counter.get_mpz_t(), mul_counter.get_mpz_t(), block / dist[phase - 1]);
        }
    }
}