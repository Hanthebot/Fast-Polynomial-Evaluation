#include "fft_finite.h"
#include <iostream>
#include <random>

using namespace std;
using namespace shk_galoiscpp;
using namespace std::chrono;

bool getRootOfUnity(const Fint& modulus, const u32& len, Fint& prim_root) {
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist(1, len);
    Fint tester = len >> 1, temp;
    u32 trial;
    for (trial = 0; trial < 1000; ++trial) {
        prim_root = dist(rng);
        // compute prim_root ^ ((modulus - 1) / 2) = 1 (mod modulus)
        // if so: it's not a primitive root
        // idea source: https://math.stackexchange.com/questions/124408/finding-a-primitive-root-of-a-prime-number
        mpz_powm(temp.get_mpz_t(), prim_root.get_mpz_t(), tester.get_mpz_t(), modulus.get_mpz_t());
        if (temp != 1) // succeed in finding rou
            break;
    }
    return (trial <= 999);
}

Field init_setup(const F& zero_F, u32& logn, u32& len) {
    cout << "Enter logn: " << endl;
    cin >> logn;
    len = 1 << logn;
    Fint modulus = len + 1;
    Fint rou;
    // alternatively: simply use rou = 3;
    assert(getRootOfUnity(modulus, len, rou) && "rou not found");
    Field my_field(modulus, 1, rou);
    cout << my_field << endl;

    // init_setup(my_field, dlog);
    return my_field;
}

void rou_init(vector<F>& w, map<F, u32>& dlog, const F& zero_F, u32 logn) {
    u32 len = 1ULL << logn;
    w.resize(len, zero_F.getOne());

    // w[0] = zero_F.getOne();
    dlog[w[0]] = 0;
    w[1] = zero_F.getRootOfUnity();
    dlog[w[1]] = 1;
    for (u32 i = 2; i < len; ++i) {
        w[i] = w[i - 1] * w[1];
        dlog[w[i]] = i;
    }
    // mul_counter += len-2;
}

void coeff_init(vector<F>& coeff, u32& deg, const F& zero_F, const vector<u32>& rev, u32 logn) {
    u32 len = 1ULL << logn;

    cout << "Enter degree of polynomial: " << endl;
    cin >> deg;
    assert(deg < len && "degree larger than length");
    
    coeff.resize(len, zero_F);
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
    // cout << coeff[0] << endl;s
}

void rev_init(vector<u32>& rev, u32 logn) {
    u32 len = 1ULL << logn;
    rev.clear();
    rev.resize(len, 0);
    rev[0] = 0;
    for (u32 i = 1; i < len; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (logn - 1);
}

void print_stats(Fint& mul_counter, microseconds duration) {
    cout << "\n===Result===\n";
    mpf_class speed, counter;

    mpf_set_si(speed.get_mpf_t(), duration.count());
    mpf_set_si(counter.get_mpf_t(), mul_counter.get_si());
    speed /= counter;

    cout << "It took " << duration.count() << " μs" << endl;
    cout << "Multiplication counter: " << mul_counter << endl;
    cout << "Speed: " << speed << " s/M mult" << endl;
}
