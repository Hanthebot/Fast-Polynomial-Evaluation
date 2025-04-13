#include <iostream>
#include "nd_fft.h"
#include "util.h"

void compute_point_sets(const vector<vector<Fint>>& points_to_evaluate, set<Fint>& val_set) {
    for (const auto& point : points_to_evaluate) {
        for (const auto& coord : point) {
            val_set.insert(coord);
        }
    }
}

bool verify_compatible(const Fint& prime, const set<Fint>& val_set, const Fint& modulo) {
    Fint temp = prime;
    while (temp < modulo) {
        if (val_set.find(temp) != val_set.end()) {
            return false;
        }
        temp += prime;
    }
    return true;
}

void compute_prime_upper_bound(const Fint& capital_M, vector<u32>& fields_used, 
    vector<bool>& isFermat, const set<Fint>& val_set, const Fint& modulo) {
    // 16 log_2 (capital_M) is a rough upper bound of prime
    Fint temp = 1; // reset temp, now for product of primes
    Fint prime = 3; // finite field with prime 2: with only 1 element, meaningless
    while (temp < capital_M) {
        if (verify_compatible(prime, val_set, modulo)) {
            temp *= prime;
            isFermat.push_back(is_Fermat_prime(prime));
            fields_used.push_back(mpz_get_ui(prime.get_mpz_t()));
        }
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }
}

u32 find_fields(const Fint& modulo, const Fint& capital_M, vector<u32>& fields_used, vector<bool>& isFermat, 
    u32& max_fft_field, const set<Fint>& val_set) {
    isFermat.clear();
    fields_used.clear();
    if (is_Fermat_prime(modulo)) {
        fields_used.push_back(modulo.get_ui());
        isFermat.push_back(true);
        max_fft_field = fields_used[0];
        return 0;
    }

    compute_prime_upper_bound(capital_M, fields_used, isFermat, val_set, modulo);
    
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
        if (fields_used[field_i] > max_fft_field) {
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

void rev_init(u32*& rev, u32 logn) {
    u32 len = 1ULL << logn;
    rev = new u32[len];
    rev[0] = 0;
    for (u32 i = 1; i < len; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (logn - 1);
}

void rev_init_NF(u32*& rev, u32**& rev_rev, const vector<u32>& radix_vec) {
    u32** revs_serial = new u32*[radix_vec.size() + 1];
    revs_serial[0] = new u32[1];
    revs_serial[0][0] = 0;
    rev_rev = new u32*[radix_vec.size() + 1];
    rev_rev[0] = new u32[1];
    rev_rev[0][0] = 0;
    u32 rev_size = 1;
    for (size_t i = 1; i <= radix_vec.size(); ++i) {
        u32 rev_size_last = rev_size;
        rev_size *= radix_vec[i - 1];
        revs_serial[i] = new u32[rev_size];
        rev_rev[i] = new u32[rev_size];
        for (size_t j = 0; j < rev_size; ++j) {
            revs_serial[i][j] = revs_serial[i - 1][j % rev_size_last] * radix_vec[i - 1] + j / rev_size_last;
            rev_rev[i][revs_serial[i][j]] = j;
        }
    }
    rev = revs_serial[radix_vec.size()];
    for (size_t i = 0; i < radix_vec.size(); ++i) {
        delete [] revs_serial[i];
    }
}