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

    compute_prime_upper_bound(capital_M, fields_used, isFermat, val_set, modulo);
    
    if (fields_used.size() == 0) {
        // impossible case
        cout << "Error: no field found" << endl;
        return 2;
    }

    // if the modulo is prime, and field used to compute is larger than prime modulo itself:
    // we can directly apply FFT on finite field of modulo, it will still be denser than the original method
    u32 modulo_ul = mpz2ul(modulo);
    if (fields_used[fields_used.size() - 1] >= modulo_ul && mpz_probab_prime_p(modulo.get_mpz_t(), 10)) {
        fields_used.clear();
        isFermat.clear();
        fields_used.push_back(modulo_ul);
        isFermat.push_back(is_Fermat_prime(modulo));
        // set validation: not done, as coefficients will be in [1, modulo) (from input validation)
    }
    
    max_fft_field = 1;
    cout << "fields used: ";
    for (const auto& field : fields_used) {
        cout << field << " ";
        if (field > max_fft_field) {
            max_fft_field = field;
        }
    }
    cout << endl;
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

void init_coeff_field(Fint*& coeff_data, const u32& m, const u64& buf_size) {
    for (size_t j = 0; j < buf_size; ++j)
        coeff_data[j] = 0;
}

void coeff_reduce_recur(const nd_vector<Fint>& arr, const nd_vector<Fint>& coeff, const u32& mod) {
    if (arr.getDim() == 0) {
        *(coeff.ptr()) += arr.get();
    }
    if (arr.getDim() == 1) {
        for (size_t i = 0; i < arr.getShape()[0]; ++i) {
            *(coeff.ptr()+(i % mod)) += arr.get(i);
        }
        return;
    }
    for (size_t i = 0; i < arr.getShape()[0]; ++i) {
        coeff_reduce_recur(arr[i], coeff[i % mod], mod);
    }
}

void coeff_reduce(const nd_vector<Fint>& arr, nd_vector<Fint>*& coeff, const u32& prime) {
    // reduce coefficients into field
    Fint prime_mp;
    ul2mpz(prime_mp, prime);
    coeff_reduce_recur(arr, *coeff, prime - 1); // for modulus
    Fint* coeff_data = coeff->ptr();
    Fint temp;
    for (size_t i = 0; i < coeff->size(); ++i) {
        // execute modulus only once per element
        if (coeff_data[i] >= prime_mp) {
            coeff_data[i] %= prime_mp;
        } else if (coeff_data[i] < 0) {
            while (coeff_data[i] < 0) {
                coeff_data[i] += prime_mp;
            }
        }
        // enforce_modulus(coeff_data[i], prime_mp);
    }
}

void rev_init(vector<u32>& rev, u32 logn) {
    u32 len = 1ULL << logn;
    rev.resize(len);
    rev[0] = 0;
    for (u32 i = 1; i < len; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (logn - 1);
}

void rev_init_NF(vector<u32>& rev, vector<vector<u32>>& rev_rev, const vector<u32>& radix_vec) {
    vector<u32>* revs_serial = new vector<u32>[radix_vec.size() + 1];
    revs_serial[0].resize(1);
    revs_serial[0][0] = 0;
    rev_rev.resize(radix_vec.size() + 1);
    rev_rev[0].resize(1);
    rev_rev[0][0] = 0;
    u32 rev_size = 1;
    for (size_t i = 1; i <= radix_vec.size(); ++i) {
        u32 rev_size_last = rev_size;
        rev_size *= radix_vec[i - 1];
        revs_serial[i].resize(rev_size);
        rev_rev[i].resize(rev_size);
        for (size_t j = 0; j < rev_size; ++j) {
            revs_serial[i][j] = revs_serial[i - 1][j % rev_size_last] * radix_vec[i - 1] + j / rev_size_last;
            rev_rev[i][revs_serial[i][j]] = j;
        }
    }
    rev = revs_serial[radix_vec.size()];
}