#include <iostream>
#include <vector>
#include <map>
#include <chrono>
#include "fft_finite.h"

using namespace std;
using namespace shk_galoiscpp;
using namespace std::chrono;

const u32 LOG_N = 16;
const u32 LEN = 1 << LOG_N;
const Fint MODULUS = LEN + 1;
const Fint ROU = 3;

int main() {
    
    Field my_field(MODULUS, 1, ROU);
    Fint mul_counter = 0;
    map<F, u32> dlog;
    
    u32 deg;
    cout << "Enter degree of polynomial: " << endl;
    cin >> deg;
    assert(deg < LEN && "degree larger than length");
    vector<F> coeff;
    vector<F> w; // w[i]: rou ^ i

    F zero_F(&my_field, 0);
    coeff.resize(LEN, zero_F);
    Fint temp;
    cout << "Enter coefficients for " << deg << " degree polynomial, from degree 0: " << endl;
    for (u32 i = 0; i <= deg; ++i) {
        cin >> temp;
        coeff[i].setX(temp);
    }

    // cout << "received: ";
    // for (u32 i = deg; i > 0; --i) {
    //     cout << coeff[i] << "x^{" << i << "} + ";
    // }
    // cout << coeff[0] << endl;

    // creating a vector of reversing bits
    vector<u32> rev;
    rev.resize(LEN, 0);
    rev[0] = 0;
    for (u32 i = 1; i < LEN; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (LOG_N - 1);

    
    w.clear();
    w.resize(LEN, coeff[0].getOne());
    w[0] = coeff[0].getOne();
    F rou = coeff[0].getRootOfUnity();
    dlog[w[0]] = 0;
    w[1] = rou;
    dlog[w[1]] = 1;
    for (u32 i = 2; i < LEN; ++i) {
        w[i] = w[i - 1] * w[1];
        dlog[w[i]] = i;
    }
    mul_counter += LEN-2;

    auto start = high_resolution_clock::now();
    fft(w, coeff, dlog, rev, LOG_N, my_field, mul_counter);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    // Fint dur = duration.count();

    cout << "It took " << duration.count() << " s" << endl;
    cout << "Multiplication counter: " << mul_counter << endl;
    // auto speed = dur / mul_counter / 1000000;
    // cout << "Speed " << speed << " s/M mult" << endl;

    // uncomment to print
    // for (F i = zero_F.getOne(); i <= LEN; ++i) {
    //     if (dlog.find(i) == dlog.end() && "i not found") break;
    //     cout << "f(" << i << "=w^" << dlog[i] << ") : " << coeff[rev[dlog[i]]] << endl;
    // }

    return 0;
}