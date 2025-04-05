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

bool getRootOfUnity(const Fint& modulus, const u32& prime, Fint& prim_root) {
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
    u64 trial;
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

/**
 * @brief Computes the power of a base raised to an exponent.
 * assumes base and exponent are non-negative integers.
 */
u64 pow(u64 base, u32 exp) {
    if (exp == 0) return 1;
    u64 result = 1;
    while (exp > 0) {
        if ((exp & 1) != 0) {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    return result;
}

/**
 * @brief Computes the power of a base raised to an exponent.
 * assumes base and exponent are non-negative integers.
 */
Fint pow(Fint base, Fint exp) {
    if (exp == 0) return 1;
    Fint result = 1;
    while (exp > 0) {
        if ((exp & 1) != 0) {
            result *= base;
        }
        base *= base;
        exp >>= 1;
    }
    return result;
}

void init_setup(Fint& modulo, vector<u32>& degs_vec, u64& total_len, u32& m) {
    cout << "Enter m: " << endl;
    cin >> m;
    modulo = 1;
    Fint MOD_UPPER_BOUND; 
    mpz_set_ui(MOD_UPPER_BOUND.get_mpz_t(), (1ULL << 32) - 1);
    while (modulo <= 1 || modulo > MOD_UPPER_BOUND) {
        cout << "Enter modulo: " << endl;
        cin >> modulo;
    }
    degs_vec.resize(m, 0);
    cout << "Enter degree of " << m << "-variate polynomial: " << endl;
    total_len = 1;
    for (size_t i = 0; i < m; ++i) {
        cin >> degs_vec[i];
        total_len *= degs_vec[i] + 1;
    }
}

void init_evaluation_points(vector<vector<Fint>>& evaluation_points, const u32& m, const Fint& modulo) {
    cout << "Enter number of evaluation points: " << endl;
    int num_points;
    cin >> num_points;
    evaluation_points.resize(num_points, vector<Fint>(m, 0));
    cout << "Enter evaluation points: " << endl;
    for (int i = 0; i < num_points; ++i) {
        for (size_t j = 0; j < m; ++j) {
            while (evaluation_points[i][j] <= 0 || evaluation_points[i][j] >= modulo) {
                cin >> evaluation_points[i][j];
            }
        }
    }
}

Field* field_setup(const u32& m, const u32& prime) {
    Fint modulus = prime;
    Fint rou;
    assert(getRootOfUnity(modulus, prime, rou) && "rou not found");
    Field* field = new Field(modulus, 1, rou);
    // cout << *field;

    return field;
}

void rou_init(F*& w, F*& dlog, const F& zero_F, const u32& prime) {
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

void recur_input(const nd_vector<Fint>& coeff, const span<const u32>& degs, const Fint& modulo) {
    if (degs.size() == 0) return;
    if (degs.size() == 1) {
        Fint temp;
        for (u32 i = 0; i <= degs[0]; ++i) {
            // in case degree < len, not using operator>> on nd_vector directly, but get()
            cin >> temp;
            *(coeff.ptr() + i) = temp % modulo;
        }
        return;
    }
    for (u32 i = 0; i <= degs[0]; ++i) {
        recur_input(coeff[i], degs.subspan(1), modulo);
    }
}

void coeff_init(const nd_vector<Fint>& coeff, const u32& m, const vector<u32>& degs_vec, const Fint& modulo) {
    cout << "Enter coefficients for <" << degs_vec[0];
    for (u32 i = 1; i < degs_vec.size(); ++i)
        cout << ", " << degs_vec[i];
    cout << "> degree polynomial, from degree 0: " << endl;
    span<const u32> degs = degs_vec; // to preserve const-ness
    cout << "modulo: " << modulo << endl;
    recur_input(coeff, degs, modulo);
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

void compute_M(Fint& capital_M, const vector<u32>& degs_vec, const Fint& modulo, 
    const u64& total_len, const u32& m) {
    /*
        we modify the original article's structure 
        to find a lower upper bound of M
        M = {prod_{i=0..m}(d_i + 1)} * (modulo - 1) ^ (sum(d_i) + 1)
        - first term represents the number of terms in the polynomial, precomputed in total_len
        - second term is the max value each term can hold
          where d_i is the degree of i-th variable
    */
    Fint sum_degs = 0;
    Fint modulo_max = modulo - 1;
    u64 fint_temp = 1ULL << 32;
    Fint UL_BOUND, mul_temp;
    mpz_set_ui(UL_BOUND.get_mpz_t(), fint_temp);
    for (size_t i = 0; i < m; ++i) {
        sum_degs += degs_vec[i];
    }
    sum_degs += 1;
    if (sum_degs < UL_BOUND) {
        mpz_pow_ui(capital_M.get_mpz_t(), modulo_max.get_mpz_t(), sum_degs.get_ui());
    } else {
        capital_M = pow(modulo_max, sum_degs);
    }
    mpz_set_ui(mul_temp.get_mpz_t(), total_len);
    capital_M *= mul_temp;
}

u64 mpz2ull(Fint& mp) {
    u64 result = 0;
    mpz_export(&result, 0, -1, sizeof(result), 0, 0, mp.get_mpz_t());
    return result;
}

void ull2mpz(Fint& mp, const u64& input) {
    mpz_import(mp.get_mpz_t(), 1, -1, sizeof(u64), 0, 0, &input);
}