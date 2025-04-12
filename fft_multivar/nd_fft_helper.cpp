#include "nd_fft.h"
#include "util.h"

void compute_point_sets(const vector<vector<Fint>>& points_to_evaluate, set<Fint>& val_set) {
    for (const auto& point : points_to_evaluate) {
        for (const auto& coord : point) {
            val_set.insert(coord);
        }
    }
}

void compute_prime_upper_bound(const Fint& capital_M, vector<u32>& fields_used) {
    // 16 log_2 (capital_M) is a rough upper bound of prime
    Fint temp = 1; // reset temp, now for product of primes
    Fint prime = 3; // finite field with prime 2: with only 1 element, meaningless
    while (temp < capital_M) {
        temp *= prime;
        fields_used.push_back(mpz_get_ui(prime.get_mpz_t()));
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }
}

void verify_fermat(const vector<u32>& fields_used, vector<bool>& isFermat, set<Fint>& val_set, const Fint& modulo) {
    isFermat.clear();
    isFermat.resize(fields_used.size(), false);
    for (size_t i = 0; i < fields_used.size(); ++i) {
        if (is_Fermat_prime(fields_used[i])) {
            Fint prime = fields_used[i];
            Fint temp = prime;
            bool fft_usable = true;
            while (temp < modulo) {
                if (val_set.find(temp) != val_set.end()) {
                    fft_usable = false;
                    break;
                }
                temp *= prime;
            }
            // i.e. safe to use for FFT
            isFermat[i] = fft_usable;
        }
    }
}

u32 find_fields(const Fint& modulo, const Fint& capital_M, vector<u32>& fields_used, vector<bool>& isFermat, 
    u32& max_fft_field, set<Fint>& val_set) {
    fields_used.clear();
    if (is_Fermat_prime(modulo)) {
        fields_used.push_back(modulo.get_ui());
        isFermat.push_back(true);
        max_fft_field = fields_used[0];
        return 0;
    }

    compute_prime_upper_bound(capital_M, fields_used);
    verify_fermat(fields_used, isFermat, val_set, modulo);
    
    cout << "fields used: ";
    for (const auto& field : fields_used) {
        cout << field << " ";
    }
    cout << endl;
    if (fields_used.size() == 0) {
        // impossible case
        cout << "Error: no field found" << endl;
        return 2;
    }
    
    max_fft_field = 1;
    for (size_t field_i = 1; field_i < fields_used.size(); ++field_i) {
        if (fields_used[field_i] > max_fft_field && isFermat[field_i]) {
            max_fft_field = fields_used[field_i];
        }
    }
    return 0;
}

void construct_nd_vector_helper(vector<size_t>& shape_vec, span<size_t>& shape, 
    vector<size_t>& unit_vec, span<size_t>& units, const u32& m, const u64& prime){
    shape_vec.clear(); unit_vec.clear();
    shape_vec.resize(m, prime - 1);
    unit_vec.resize(m, 1);
    for (int i = m - 2; i >= 0; --i) {
        unit_vec[i] = unit_vec[i + 1] * shape_vec[i + 1];
    }
    shape = shape_vec;
    units = unit_vec;
}

void init_coeff_field(F*& coeff_data, const u32& m, const u64& buf_size, const F& zero_F) {
    for (size_t j = 0; j < buf_size; ++j)
        coeff_data[j] = zero_F.getZero();
}

void coeff_reduce_recur(const nd_vector<Fint>& arr, const nd_vector<F>& coeff, const u32& mod) {
    if (arr.getDim() == 0) {
        coeff.ptr()->raw_add(arr.get());
    }
    if (arr.getDim() == 1) {
        for (size_t i = 0; i < arr.getShape()[0]; ++i) {
            (coeff.ptr()+(i % mod))->raw_add(arr.get(i));
        }
        return;
    }
    for (size_t i = 0; i < arr.getShape()[0]; ++i) {
        coeff_reduce_recur(arr[i], coeff[i % mod], mod);
    }
}

void coeff_reduce(const nd_vector<Fint>& arr, nd_vector<F>*& coeff, const u32& prime) {
    // reduce coefficients into field
    coeff_reduce_recur(arr, *coeff, prime - 1); // for modulus
    F* coeff_data = coeff->ptr();
    for (size_t i = 0; i < coeff->size(); ++i) {
        // execute modulus only once per element
        coeff_data[i].enforce_modulus();
    }
}