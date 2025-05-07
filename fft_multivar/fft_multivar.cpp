#include <iostream>
#include <cmath>
#include "fft_multivar.h"
#include "util.h"

using namespace std;

int main(int argc, char* argv[]) {
    vector<u32> degs_vec;
    vector<vector<Fint>> evaluation_points;
    vector<Fint> results;
    Fint mul_counter = 0, modulo = 0, capital_M;
    u64 total_len;
    u32 m;
    
    bool print_coeff, print_rst, print_verify_rst;
    cout << "Print coefficients? (0/1): ";
    cin >> print_coeff;
    cout << "Print result? (0/1): ";
    cin >> print_rst;
    cout << "Verify result? (0/1): ";
    cin >> print_verify_rst;

    EvalIO meta;
    if (argc > 1) {
        meta.print_time = atoi(argv[1]);
    }
    if (argc > 2) {
        meta.debug = atoi(argv[2]);
    }
    // initializing field and consts
    init_setup(modulo, degs_vec, total_len, m);

    compute_M(capital_M, degs_vec, modulo, total_len, m);
    // Field my_field = init_setup(logn, len, m, prime);
    // zero_F.setField(&my_field);
    init_evaluation_points(evaluation_points, m, modulo);

    // receiving coefficients
    vector<size_t> shape_vec(m, 0);
    for (size_t i = 0; i < m; ++i) {
        shape_vec[i] = degs_vec[i] + 1;
    }
    span<size_t> shape{shape_vec};
    vector<size_t> units_vec(m, 1);
    units_vec[m - 1] = 1;
    for (int i = m - 2; i >= 0; --i) {
        units_vec[i] = units_vec[i + 1] * shape_vec[i + 1];
    }
    span<size_t> units{units_vec};
    Fint* vec_data = new Fint[total_len];
    nd_vector<Fint> coeff {m, shape, units, vec_data, total_len};

    coeff_init(coeff, m, degs_vec, modulo);
    
    if (print_coeff) {
        cout << "===coeff===\n" << coeff << endl;
    }
    auto start = chrono::high_resolution_clock::now();
    {
        int result = evaluate_all_point(coeff, modulo, capital_M, mul_counter, evaluation_points, results, meta);
        switch (result) {
            case 0:
                cout << "Success!" << endl;
                break;
            case 1:
                cout << "Success with Fermat prime!" << endl;
                break;
            case 2:
                cout << "Failure: non-32-bit prime not supported yet!" << endl;
                break;
        }
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    print_stats(mul_counter, duration);
    if (meta.print_time) {
        cout << "Time taken for precomputable: " << from_micro(meta.times.precomp) << " s" << endl
        << "Time taken for memory cleaning: " << from_micro(meta.times.mem_clean) << " s" << endl
        << "Time taken for FFT: " << from_micro(meta.times.fft) << " s" << endl
        << "Time taken for retrieval: " << from_micro(meta.times.retrieval) << " s" << endl
        << "Time taken for CRT: " << from_micro(meta.times.crt) << " s" << endl;
    }

    if (print_rst) {
        long long int naive_time = 0;
        cout << "\n===computation===" << endl;
        for (size_t i = 0; i < results.size(); ++i) {
            // mul_counter: doesn't matter anymore
            cout << "f(" << evaluation_points[i][0];
            for (size_t j = 1; j < m; ++j) {
                cout << ", " << evaluation_points[i][j];
            }
            cout << ") = " << results[i] << " ";
            if (print_verify_rst) {
                start = chrono::high_resolution_clock::now();
                Fint brute_eval = evaluate_brutal(coeff, evaluation_points[i], mul_counter) % modulo;
                accumulate_time(start, stop, naive_time);
                if (results[i] != brute_eval) {
                    cout << "!= " << brute_eval;
                } else {
                    cout << "== " << brute_eval;
                }
            }
            cout << endl;
        }
        mpf_class time_spent = from_micro(naive_time);
        double res_size = (double) results.size();
        if (res_size > 0) {
            time_spent /= res_size;
        }
        cout << "Time taken for naive: " << time_spent << " s" << endl;
    }

    return 0;
}