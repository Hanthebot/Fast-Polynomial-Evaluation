#include "fft_multivar.h"
#include "nd_fft.h"
#include "util.h"

using namespace std;

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
