#include "fft_multivar.h"
#include <set>

using namespace std;

void fft(const F* w, F* const& arr, const u32* rev, u32 logn, Fint& mul_counter, F& u, F& v);
void fft_multivar_recur(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);
void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const F* dlog, const u32* rev, u32 logn, Fint& mul_counter, const F& zero_F, F*& temp_u_vec, F*& temp_v_vec);
void fft_nd(const F* w, const nd_vector<F>& arr, const u32* rev, u32 logn, Fint& mul_counter, nd_vector<F>& temp_u, nd_vector<F>& temp_v);

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
bool is_Fermat_prime(const Fint& modulo) {
    // fermat prime: 2 ^ {2^n} + 1
    // n = 0, 1, 2, 3, 4
    return modulo == 3 || modulo == 5 || modulo == 17 || modulo == 257 || modulo == 65537;
}

u32 log2(u64 len) {
    // 2^logn + 1 = modulo
    // logn = log2(modulo - 1)
    // len = 2^logn
    u32 logn = 0;
    while (len > 1) {
        len >>= 1;
        ++logn;
    }
    return logn;
}

void compute_point_sets(const vector<vector<Fint>>& points_to_evaluate, set<Fint>& val_set) {
    for (const auto& point : points_to_evaluate) {
        for (const auto& coord : point) {
            val_set.insert(coord);
        }
    }
}

void compute_prime_upper_bound(const Fint& capital_M, vector<u32>& fields_used) {
    // we perform rough estimation of log_2 (capital_M)
    Fint temp = capital_M;
    // Fint upper_bound = 1;
    // while (temp > 1) {
    //     temp >>= 1;
    //     ++upper_bound;
    // }
    // upper_bound *= 16; // 16 log_2 (capital_M) is a rough upper bound of prime
    temp = 1; // reset temp, now for product of primes
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

int evaluate_all_point(const nd_vector<Fint>& arr, const Fint& modulo, const Fint& capital_M, Fint& mul_counter, const vector<vector<Fint>>& points_to_evaluate, vector<Fint>& results) {
    // figure out which fields to use
    u32 m = arr.getDim(), max_fft_field, field_signal;
    vector<bool> isFermat;
    vector<u32> fields_used;
    set<Fint> val_set;
    compute_point_sets(points_to_evaluate, val_set);
    
    if ((field_signal = find_fields(modulo, capital_M, fields_used, isFermat, max_fft_field, val_set))) {
        cout << "Error: from finding fields " << field_signal << endl;
        return 2;
    }
    
    Field** fields = new Field*[fields_used.size()];
    F** w_s = new F*[fields_used.size()];
    F** dlog_s = new F*[fields_used.size()];
    F* zero_F_s = new F[fields_used.size()];
    F* coeff_data = new F[pow(max_fft_field - 1, m)]; // for reduced coefficients
    Fint** point_evaluation = new Fint*[fields_used.size()];
    nd_vector<F>** coeff_s = new nd_vector<F>*[fields_used.size()];
    u32** rev_s = new u32*[fields_used.size()];
    vector<size_t>* shape_vec_s = new vector<size_t>[fields_used.size()];
    vector<size_t>* unit_vec_s = new vector<size_t>[fields_used.size()];
    
    F* temp_u_vec = new F[pow(max_fft_field - 1, m - 1)];
    F* temp_v_vec = new F[pow(max_fft_field - 1, m - 1)];
    
    for (size_t i = 0; i < fields_used.size(); ++i) {
        span<size_t> shape, units;
        u64 buf_size = pow(fields_used[i] - 1, m);
        
        fields[i] = field_setup(m, fields_used[i]);
        zero_F_s[i].setField(fields[i]);
        rou_init(w_s[i], dlog_s[i], zero_F_s[i], fields_used[i]);
        fields[i]->setDlogW(dlog_s[i], w_s[i]);
        point_evaluation[i] = new Fint[points_to_evaluate.size()];
        if (isFermat[i]) {
            // receiving coefficients
            construct_nd_vector_helper(shape_vec_s[i], shape, unit_vec_s[i], units, m, fields_used[i]);
            
            init_coeff_field(coeff_data, arr.getDim(), buf_size, zero_F_s[i]);
    
            coeff_s[i] = new nd_vector<F>(m, shape, units, coeff_data, buf_size);
            coeff_reduce(arr, coeff_s[i], fields_used[i]);

            u32 local_logn = log2(fields_used[i] - 1);
            rev_init(rev_s[i], local_logn);
            fft_multivar_wrapper(w_s[i], *coeff_s[i], dlog_s[i], rev_s[i], local_logn, 
                mul_counter, zero_F_s[i], temp_u_vec, temp_v_vec);
            for (size_t j = 0; j < points_to_evaluate.size(); ++j) {
                nd_vector<F> elem = *coeff_s[i];
                for (const Fint& coord : points_to_evaluate[j]) {
                    Fint coord_mod = coord % fields_used[i];
                    u64 coord_int = mpz_get_ui(coord_mod.get_mpz_t());
                    elem = elem[(dlog_s[i][coord_int]).toSize()];
                }
                point_evaluation[i][j] = elem.get().getX();
            }
        } else {
            for (size_t j = 0; j < points_to_evaluate.size(); ++j) {
                point_evaluation[i][j] = evaluate_brutal(arr, points_to_evaluate[j], mul_counter);
                point_evaluation[i][j] %= fields_used[i];
            }
        }
    }

    // reconstruct the solution using Chinese remainder theorem
    vector<Fint> interpolate(fields_used.size(), Fint(0));
    Fint common_prod = 1;
    for (size_t i = 0; i < fields_used.size(); ++i) {
        common_prod *= fields[i]->getModulus();
    }
    
    for (size_t i = 0; i < fields_used.size(); ++i) {
        Fint temp = common_prod / fields[i]->getModulus();
        Fint inv = multInverse(temp % fields[i]->getModulus(), fields[i]);
        interpolate[i] = inv * temp;
    }

    results.reserve(points_to_evaluate.size());
    for (size_t j = 0; j < points_to_evaluate.size(); ++j) {
        if (DEBUG) {
            cout << "(" << points_to_evaluate[j][0];
            for (size_t i = 1; i < points_to_evaluate[j].size(); ++i) {
                cout << ", " << points_to_evaluate[j][i];
            }
            cout << ") :\n";
        }
        Fint temp = 0;
        for (size_t i = 0; i < fields_used.size(); ++i) {
            temp += interpolate[i] * point_evaluation[i][j];
            if (DEBUG) {
                cout << "\t" << fields_used[i] << ": " << point_evaluation[i][j] << endl;
            }
        }
        if (DEBUG) {
            cout << "\toverall: " << temp << " w/ mod " << common_prod;
        }
        temp %= common_prod;
        results.push_back(temp % modulo);
        if (DEBUG) {
            cout << " with mod " << modulo << " => " << results.back() << endl;
        }
    }

    return (fields_used.size() == 1 && isFermat[0]) ? 1 : 0; // success
}

void fft_multivar_wrapper(const F* w, nd_vector<F>& arr, const F* dlog, const u32* rev, u32 logn, Fint& mul_counter, const F& zero_F,
    F*& temp_u_vec, F*& temp_v_vec) {
    u32 len = 1ULL << logn;
    for (const auto& s : arr.getShape()) {
        assert(s == len && "size does not match");
    }
    span<size_t> shape = arr.getShape().subspan(1);
    span<size_t> units = arr.getUnit().subspan(1);
    nd_vector<F> temp_u {arr.getDim() - 1, shape, units, temp_u_vec, arr[0].size()};
    nd_vector<F> temp_v {arr.getDim() - 1, shape, units, temp_v_vec, arr[0].size()};
    fft_multivar_recur(w, arr, rev, logn, mul_counter, temp_u, temp_v);
}

void fft_multivar_recur(const F* w, const nd_vector<F>& arr, 
    const u32* rev, u32 logn, Fint& mul_counter,
    nd_vector<F>& temp_u, nd_vector<F>& temp_v) {
    assert(isEqual(arr.getShape().subspan(1), temp_u.getShape()) && "shape does not match");
    assert(isEqual(arr.getShape().subspan(1), temp_v.getShape()) && "shape does not match");
    if (arr.getDim() == 1) {
        F u = temp_u.get();
        F v = temp_v.get();
        fft(w, arr.ptr(), rev, logn, mul_counter, u, v);
        return;
    }
    // for temporary storage of (n-1)d array
    nd_vector<F> temp_u_ = temp_u[0];
    nd_vector<F> temp_v_ = temp_v[0];
    for (size_t i = 0; i < arr.getShape()[0]; ++i) {
        fft_multivar_recur(w, arr[i], rev, logn, mul_counter, temp_u_, temp_v_);
    }
    fft_nd(w, arr, rev, logn, mul_counter, temp_u, temp_v);
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
    for (u32 i = 0; i < len; ++i) {
        if (i < rev[i]) {
            arr[i].swap(arr[rev[i]]);
        }
    }
}

Fint evaluate_nd(const nd_vector<Fint>& arr, const span<const Fint>& points, Fint& mul_counter) {
    if (arr.getDim() == 0) {
        return arr.get();
    }
    Fint result = 0;
    for (size_t i = arr.getShape()[0] - 1; ; --i) {
        result *= points[0];
        result += evaluate_nd(arr[i], points.subspan(1), mul_counter);
        if (i == 0) break;
    }
    return result;
}

Fint evaluate_brutal(const nd_vector<Fint>& arr, const vector<Fint>& points, Fint& mul_counter) {
    span<const Fint> points_span = points;
    Fint result = evaluate_nd(arr, points_span, mul_counter);
    Fint acc_count = 1, temp;
    for (size_t i = 0; i < arr.getDim(); ++i) {
        mpz_mul_ui(acc_count.get_mpz_t(), acc_count.get_mpz_t(), arr.getShape()[i]);
        mul_counter += acc_count;
    }
    
    return result;
}