#ifndef ND_FFT_H
#define ND_FFT_H

#include "fft_multivar.h"
#include <set>

// fft_compute.cpp
void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const F* dlog, const u32* rev, u32 logn, Fint& mul_counter, const F& zero_F, F*& temp_u_vec, F*& temp_v_vec);

// nd_fft_helper.cpp
void compute_point_sets(const vector<vector<Fint>>& points_to_evaluate, set<Fint>& val_set);
u32 find_fields(const Fint& modulo, const Fint& capital_M, vector<u32>& fields_used, vector<bool>& isFermat, u32& max_fft_field, set<Fint>& val_set);
void construct_nd_vector_helper(vector<size_t>& shape_vec, span<size_t>& shape, vector<size_t>& unit_vec, span<size_t>& units, const u32& m, const u64& prime);
void init_coeff_field(F*& coeff_data, const u32& m, const u64& buf_size, const F& zero_F);
void coeff_reduce(const nd_vector<Fint>& arr, nd_vector<F>*& coeff, const u32& prime);

#endif // ND_FFT_H