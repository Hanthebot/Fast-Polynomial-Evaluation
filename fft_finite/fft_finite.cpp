#include <iostream>
#include "fft_finite.h"

using namespace std;
using namespace shk_galoiscpp;

int main() {
    map<F, u32> dlog;
    vector<u32> rev;
    vector<F> coeff, w;
    F zero_F(nullptr, 0);
    Fint mul_counter = 0;
    u32 logn, len, deg;
    
    bool print_coeff, print_rst;
    cout << "Print coefficients? (0/1): ";
    cin >> print_coeff;
    cout << "Print result? (0/1): ";
    cin >> print_rst;
    // initializing field and consts
    Field my_field = init_setup(zero_F, logn, len);
    zero_F.setField(&my_field);
    
    // creating a vector of reversing bits
    rev_init(rev, logn);

    // receiving coefficients
    coeff_init(coeff, deg, zero_F, rev, logn);

    // initializing rou and dlog map
    rou_init(w, dlog, zero_F, logn);

    if (print_coeff) {
        cout << "===coeff===\n";
        cout << "f(x) = ";
        for (u32 i = 0; i < len; ++i) {
            cout <<  coeff[i] << "x^{" << i << "} ";
        }
        cout << endl;
    }

    auto start = chrono::high_resolution_clock::now();
    fft(w, coeff, dlog, rev, logn, my_field, mul_counter);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    print_stats(mul_counter, duration);

    if (print_rst) {
        u32 counter = 0;
        for (F i = zero_F.getOne(); counter < len; ++i) {
            assert(dlog.find(i) != dlog.end() && "i not found");
            cout << "f(" << i << "=w^" << dlog[i] << ") : " << coeff[dlog[i]] << endl;
            ++counter;
        }
    }

    return 0;
}