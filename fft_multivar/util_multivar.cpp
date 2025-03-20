#include "fft_multivar.h"
#include <iostream>
#include <random>

using namespace std;
using namespace shk_galoiscpp;
using namespace std::chrono;

void factorize(Fint n, vector<Fint>& factors) {
    Fint factor = 2;
    while (factor * factor <= n) {
        if (mpz_divisible_p(n.get_mpz_t(), factor.get_mpz_t())) {
            factors.push_back(factor);
            while (mpz_divisible_p(n.get_mpz_t(), factor.get_mpz_t())) {
                mpz_tdiv_q(n.get_mpz_t(), n.get_mpz_t(), factor.get_mpz_t());
            }
        }
        mpz_nextprime(factor.get_mpz_t(), factor.get_mpz_t());
    }
    if (n > 1) {
        factors.push_back(n);
    }
}

bool getRootOfUnity(const Fint& modulus, const u32& len, const u32& prime, Fint& prim_root) {
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist(1, prime - 1);
    // finding primitive root of prime using Euler's Totient Theorem
    // reference: https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
    Fint prime_mpz = prime;
    assert(mpz_probab_prime_p(prime_mpz.get_mpz_t(), 10) && "prime not prime");
    prime_mpz -= 1; // totient
    vector<Fint> factors;
    factorize(prime_mpz, factors);
    Fint tester, temp;
    u32 trial;
    for (trial = 0; trial < 1000; ++trial) {
        prim_root = dist(rng);
        for (const auto& factor : factors) {
            mpz_tdiv_q(tester.get_mpz_t(), prime_mpz.get_mpz_t(), factor.get_mpz_t());
            mpz_powm(temp.get_mpz_t(), prim_root.get_mpz_t(), tester.get_mpz_t(), modulus.get_mpz_t());
            if (temp == 1) // fail in finding rou
                break;
        }
        if (temp != 1) // success in finding rou
            break;
    }
    return (trial <= 999);
}

Field init_setup(u32& logn, u32& len, size_t& m, u32& prime) {
    cout << "Enter m: " << endl;
    cin >> m;
    cout << "Enter prime, upper bound of deg is prime - 2: " << endl;
    cin >> prime;
    logn = 1;
    while ((1ULL << logn) < prime - 1)
        ++logn;
    len = 1 << logn;
    Fint modulus = prime;
    Fint rou;
    // alternatively: simply use rou = 3;
    assert(getRootOfUnity(modulus, len, prime, rou) && "rou not found");
    Field my_field(modulus, 1, rou);
    cout << my_field << endl;

    // init_setup(my_field, dlog);
    return my_field;
}

void rou_init(F*& w, F*& dlog, const F& zero_F, u32 prime) {
    w = new F[prime - 1];
    dlog = new F[prime];
    w[0] = zero_F.getOne();
    dlog[w[0].toSize()] = zero_F;
    w[1] = zero_F.getRootOfUnity();
    dlog[w[1].toSize()] = zero_F.getOne();
    Fint count = 2;
    for (u32 i = 2; i < prime - 1; ++i) {
        w[i] = w[i - 1] * w[1];
        dlog[w[i].toSize()] = count;
        count++;
    }
    // mul_counter += len-2;
}

void recur_input(const nd_vector<F>& coeff, const span<u32>& degs, Field* const field) {
    if (degs.size() == 0) return;
    if (degs.size() == 1) {
        Fint temp;
        for (u32 i = 0; i <= degs[0]; ++i) {
            // in case degree < len, not using operator>> on nd_vector directly, but get()
            cin >> temp;
            coeff[i].get() = {field, temp};
        }
        for (u32 i = degs[0] + 1; i < coeff.getShape()[0]; ++i) {
            coeff[i].get() = {field, 0};
        }
        return;
    }
    for (u32 i = 0; i <= degs[0]; ++i) {
        recur_input(coeff[i], degs.subspan(1), field);
    }
}

void coeff_init(const nd_vector<F>& coeff, const size_t& m, vector<u32>& degs_vec, u32 prime, Field* const field) {
    degs_vec.resize(m, 0);
    cout << "Enter degree of " << m << "-variate polynomial: " << endl;
    for (size_t i = 0; i < m; ++i) {
        cin >> degs_vec[i];
        assert(degs_vec[i] < prime - 1 && "degree larger than length");
    }
    
    cout << "Enter coefficients for <" << degs_vec[0];
    for (u32 i = 1; i < m; ++i)
        cout << ", " << degs_vec[i];
    cout << "> degree polynomial, from degree 0: " << endl;
    span<u32> degs{degs_vec};
    recur_input(coeff, degs, field);
}

void rev_init(u32*& rev, u32 logn) {
    u32 len = 1ULL << logn;
    rev = new u32[len];
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

void print_dlog(const nd_vector<F>& coeff, F* dlog, const F& zero_F, const u32& prime, const string& prefix) {
    if (coeff.getDim() == 0) {
        cout << prefix << coeff.get() << "\n";
        return;
    }
    if (coeff.getDim() == 1) {
        u32 counter = 1;
        for (F i = zero_F.getOne(); counter < prime; ++i) {
            // assert(dlog.find(i) != dlog.end() && "i not found");
            cout << prefix << i << ") : " << coeff[dlog[counter].toSize()].get() << endl;
            ++counter;
        }
        cout << endl;
        return;
    } else {
        u32 counter = 1;
        for (F i = zero_F.getOne(); counter < prime; ++i) {
            print_dlog(coeff[dlog[counter].toSize()], dlog, zero_F, prime, prefix + i.getX().get_str() + ", ");
            ++counter;
        }
        cout << endl;
        return;
    }
}
