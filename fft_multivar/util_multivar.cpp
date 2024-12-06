#include "fft_multivar.h"
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

Field init_setup(const F& zero_F, u32& logn, u32& len, size_t& m) {
    cout << "Enter m: " << endl;
    cin >> m;
    cout << "Enter logn, upper bound of deg is 2^logn: " << endl;
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

void recur_input(const nd_vector<F>& coeff, const span<u32>& degs) {
    if (degs.size() == 0) return;
    if (degs.size() == 1) {
        Fint temp;
        for (u32 i = 0; i <= degs[0]; ++i) {
            // in case degree < len, not using operator>> on nd_vector directly, but get()
            cin >> temp;
            coeff[i].get() = temp;
        }
        return;
    }
    for (u32 i = 0; i <= degs[0]; ++i) {
        recur_input(coeff[i], degs.subspan(1));
    }
}

void coeff_init(const nd_vector<F>& coeff, const size_t& m, vector<u32>& degs_vec, const F& zero_F, const vector<u32>& rev, u32 logn) {
    u32 len = 1ULL << logn;

    degs_vec.resize(m, 0);
    cout << "Enter degree of " << m << "-variate polynomial: " << endl;
    for (size_t i = 0; i < m; ++i) {
        cin >> degs_vec[i];
        assert(degs_vec[i] < len && "degree larger than length");
    }
    
    Fint temp;
    cout << "Enter coefficients for <" << degs_vec[0];
    for (u32 i = 1; i < m; ++i)
        cout << ", " << degs_vec[i];
    cout << "> degree polynomial, from degree 0: " << endl;
    span<u32> degs{degs_vec};
    recur_input(coeff, degs);
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

void print_dlog(const nd_vector<F>& coeff, map<F, u32>& dlog, const F& zero_F, const string& prefix) {
    if (coeff.getDim() == 0) {
        cout << prefix << coeff.get() << "\n";
        return;
    }
    if (coeff.getDim() == 1) {
        u32 counter = 0;
        for (F i = zero_F.getOne(); counter < coeff.getShape()[0]; ++i) {
            assert(dlog.find(i) != dlog.end() && "i not found");
            cout << prefix << i << ") : " << coeff[dlog[i]].get() << endl;
            ++counter;
        }
        cout << endl;
        return;
    } else {
        u32 counter = 0;
        for (F i = zero_F.getOne(); counter < coeff.getShape()[0]; ++i) {
            assert(dlog.find(i) != dlog.end() && "i not found");
            print_dlog(coeff[dlog[i]], dlog, zero_F, prefix + i.getX().get_str() + ", ");
            ++counter;
        }
        cout << endl;
        return;
    }
}

// void compute_manual(nd_vector<F>& ans, const nd_vector<F>& coeff, const F& carry) {
//     if (coeff.getDim() == 0) {
//         return;
//     } else if (coeff.getDim() == 1) {
//         F temp = carry.getOne();
//         for (size_t i = i; i <= ans.getShape()[0]; ++i) {
//             ans.get(i) += carry * temp;
//             temp += 1;
//         }
//         return;
//     }
//     if (coeff.getDim() == 1) {
//         for (size_t i = 0; i < coeff.getShape()[0]; ++i) {
//             ans[i].get() = coeff[i].get();
//         }
//         return;
//     } else {
//         for (size_t i = 0; i < coeff.getShape()[0]; ++i) {
//             compute_manual(ans[i] *= coeff[i]);
//         }
//         return;
//     }
// }
