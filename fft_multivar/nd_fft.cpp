#include "fft_multivar.h"
#include "nd_fft.h"
#include "util.h"
#include "field_data.h"
#include <omp.h>

using namespace std;

int evaluate_all_point(const nd_vector<Fint>& arr, const Fint& modulo, const Fint& capital_M, Fint& mul_counter, const vector<vector<Fint>>& points_to_evaluate, vector<Fint>& results, EvalIO& meta) {
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
    
    FieldData* FD = new FieldData[fields_used.size()]; // field data
    
    Fint** point_evaluation = new Fint*[fields_used.size()];
    Fint** coeff_data_s = new Fint*[fields_used.size()]; // for reduced coefficients
    nd_vector<Fint>** coeff_s = new nd_vector<Fint>*[fields_used.size()];
    vector<size_t>* shape_vec_s = new vector<size_t>[fields_used.size()];
    vector<size_t>* unit_vec_s = new vector<size_t>[fields_used.size()];
    
    Fint* temp_u_vec = new Fint[pow(max_fft_field - 1, m)]; // temporary storage for FFT
    Fint* temp_v1_vec = new Fint[pow(max_fft_field - 1, m - 1)]; // temporary storage for value storage
    Fint* temp_v2_vec = new Fint[pow(max_fft_field - 1, m - 1)]; // temporary storage for value storage
    
    chrono::high_resolution_clock::time_point start, end;
    long long int duration = 0;
    for (size_t i = 0; i < fields_used.size(); ++i) {
        start = chrono::high_resolution_clock::now();
        span<size_t> shape, units;
        u64 buf_size = pow(fields_used[i] - 1, m);
        
        FD[i].prime = fields_used[i];
        ul2mpz(FD[i].prime_mp, fields_used[i]);
        factorize(FD[i].prime - 1, FD[i].radix_vec);
        if (!getRootOfUnity(FD[i].prime_mp, FD[i].prime, FD[i].rou)) {
            cout << "Error: root of unity not found for field " << fields_used[i] << endl;
            return 2;
        }
        rou_init(FD[i].w, FD[i].dlog, FD[i].prime, FD[i].prime_mp, FD[i].rou);
        
        if (isFermat[i]) {
            rev_init(FD[i].rev, log2(FD[i].prime - 1));
        } else {
            rev_init_NF(FD[i].rev, FD[i].rev_rev, FD[i].radix_vec);
            FD[i].dist.resize(FD[i].radix_vec.size() + 1, 0);
            FD[i].dist[0] = FD[i].prime - 1;
            for (size_t dist_i = 1; dist_i <= FD[i].radix_vec.size(); ++dist_i) {
                FD[i].dist[dist_i] = FD[i].dist[dist_i - 1] / FD[i].radix_vec[dist_i - 1];
            }
        }
        accumulate_time(start, end, meta.times.precomp);
        
        // for ease of use
        Fint* w = FD[i].w.data();
        point_evaluation[i] = new Fint[points_to_evaluate.size()];
        construct_nd_vector_helper(shape_vec_s[i], shape, unit_vec_s[i], units, m, FD[i].prime);
        
        coeff_data_s[i] = new Fint[pow(max_fft_field - 1, m)];
        Fint* coeff_data = coeff_data_s[i];
        init_coeff_field(coeff_data, arr.getDim(), buf_size);
        
        coeff_s[i] = new nd_vector<Fint>(m, shape, units, coeff_data, buf_size);
        coeff_reduce(arr, coeff_s[i], FD[i].prime);
        accumulate_time(start, end, meta.times.mem_clean);
        if (isFermat[i]) {
            u32 local_logn = log2(FD[i].prime - 1);
            fft_multivar_wrapper(w, *coeff_s[i], FD[i].rev.data(), local_logn, 
                mul_counter, temp_v1_vec, temp_v2_vec);
        } else {
            fft_multivar_wrapper_NF(w, *coeff_s[i], FD[i].radix_vec, FD[i].dist, FD[i].rev_rev.data(),
            mul_counter, temp_u_vec, temp_v1_vec, temp_v2_vec);
        }
        duration = accumulate_time(start, end, meta.times.fft);
        if (meta.print_time) {
            cout << "Field " << fields_used[i] << " evaluated in " << duration << " μs" << endl;
        }
    }
    for (size_t i = 0; i < fields_used.size(); ++i) {
        start = chrono::high_resolution_clock::now();
        Fint* dlog = FD[i].dlog.data();
        for (size_t j = 0; j < points_to_evaluate.size(); ++j) {
            nd_vector<Fint> elem = *coeff_s[i];
            for (const Fint& coord : points_to_evaluate[j]) {
                Fint coord_mod = coord % FD[i].prime;
                u64 coord_int = mpz2ull(coord_mod);
                elem = elem[FD[i].rev[mpz2ul(dlog[coord_int])]];
            }
            enforce_modulus(elem.get(), FD[i].prime);
            point_evaluation[i][j] = elem.get();
        }
        accumulate_time(start, end, meta.times.retrieval);
    }
    
    // reconstruct the solution using Chinese remainder theorem
    start = chrono::high_resolution_clock::now();
    vector<Fint> interpolate(fields_used.size(), Fint(0));
    Fint common_prod = 1;
    for (size_t i = 0; i < fields_used.size(); ++i) {
        common_prod *= FD[i].prime_mp;
    }
    
    for (size_t i = 0; i < fields_used.size(); ++i) {
        Fint temp = common_prod / FD[i].prime_mp;
        Fint inv = multInverse(temp % FD[i].prime_mp, FD[i].prime_mp, FD[i].dlog.data(), FD[i].w.data());
        interpolate[i] = inv * temp;
    }

    results.reserve(points_to_evaluate.size());
    for (size_t j = 0; j < points_to_evaluate.size(); ++j) {
        if (meta.debug) {
            cout << "(" << points_to_evaluate[j][0];
            for (size_t i = 1; i < points_to_evaluate[j].size(); ++i) {
                cout << ", " << points_to_evaluate[j][i];
            }
            cout << ") :\n";
        }
        Fint temp = 0;
        for (size_t i = 0; i < fields_used.size(); ++i) {
            temp += interpolate[i] * point_evaluation[i][j];
            if (meta.debug) {
                cout << "\t" << fields_used[i] << ": " << point_evaluation[i][j];
                Fint t = 0;
                Fint brute = evaluate_brutal(arr, points_to_evaluate[j], t) % fields_used[i];
                cout << ((point_evaluation[i][j] == brute) ? " == " : " != ") << brute << endl;
            }
        }
        if (meta.debug) {
            cout << "\toverall: " << temp << " w/ mod " << common_prod;
        }
        temp %= common_prod;
        results.push_back(temp % modulo);
        if (meta.debug) {
            cout << " with mod " << modulo << " => " << results.back() << endl;
        }
    }
    accumulate_time(start, end, meta.times.crt);

    return (fields_used.size() == 1 && isFermat[0]) ? 1 : 0; // success
}
