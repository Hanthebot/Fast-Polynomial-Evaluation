#include <iostream>
#include <cmath>
#include "fft_multivar.h"

using namespace std;
using namespace shk_galoiscpp;

int main() {
    map<F, u32> dlog;
    vector<u32> rev, degs_vec;
    vector<F> w;
    F zero_F(nullptr, 0);
    Fint mul_counter = 0;
    u32 logn, len;
    size_t m;
    
    bool print_coeff, print_rst;
    cout << "Print coefficients? (0/1): ";
    cin >> print_coeff;
    cout << "Print result? (0/1): ";
    cin >> print_rst;
    // initializing field and consts
    Field my_field = init_setup(zero_F, logn, len, m);
    zero_F.setField(&my_field);
    
    // creating a vector of reversing bits
    rev_init(rev, logn);

    // receiving coefficients
    vector<size_t> shape_vec(m, len);
    span<size_t> shape{shape_vec};
    size_t total_len = static_cast<size_t>(pow(len, m));
    vector<F> vec_data (total_len, zero_F);
    nd_vector<F> coeff {m, shape, vec_data};

    coeff_init(coeff, m, degs_vec, zero_F, rev, logn);
    
    // initializing rou and dlog map
    rou_init(w, dlog, zero_F, logn);

    if (print_coeff) {
        cout << "===coeff===\n" << coeff << endl;
    }

    // manual computation
    vector<F> vec_verify_data (total_len, zero_F);
    nd_vector<F> arr {m, shape, vec_verify_data};

    auto start = chrono::high_resolution_clock::now();
    Fint mul_counter_v = 0;
    compute_manual(coeff, arr, mul_counter_v);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    print_stats(mul_counter_v, duration);

    start = chrono::high_resolution_clock::now();
    fft_multivar_wrapper(w, coeff, dlog, rev, logn, mul_counter, zero_F);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    print_stats(mul_counter, duration);

    if (print_rst) {
        cout << "\n===computation===" << endl;
        print_dlog(coeff, dlog, zero_F, "f(");
    }

    Fint total = 0, incorrect = 0;
    cout << "\n===Verification===" << endl;
    check_result(coeff, arr, dlog, zero_F, total, incorrect);
    cout << "Total: " << total << ", Incorrect: " << incorrect << endl;

    return 0;
}