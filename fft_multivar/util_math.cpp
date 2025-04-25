#include "util.h"
#include <iostream>
#include <random>

using namespace std;
using namespace std::chrono;

void get_all_factors(Fint n, vector<Fint>& factors) {
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

void factorize(Fint n, vector<Fint>& factors) {
    Fint factor = 2;
    while (factor * factor <= n) {
        if (mpz_divisible_p(n.get_mpz_t(), factor.get_mpz_t())) {
            while (mpz_divisible_p(n.get_mpz_t(), factor.get_mpz_t())) {
                mpz_tdiv_q(n.get_mpz_t(), n.get_mpz_t(), factor.get_mpz_t());
                factors.push_back(factor);
            }
        }
        mpz_nextprime(factor.get_mpz_t(), factor.get_mpz_t());
    }
    if (n > 1) {
        factors.push_back(n);
    }
}

void factorize(Fint n, vector<u32>& factors) {
    vector<Fint> factors_mpz;
    factorize(n, factors_mpz);
    factors.resize(factors_mpz.size(), 0);
    for (size_t i = 0; i < factors_mpz.size(); ++i) {
        factors[i] = mpz2ul(factors_mpz[i]);
    }
}

bool getRootOfUnity(const Fint& modulus, const u32& prime, Fint& prim_root) {
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist(1, prime - 1);
    // finding primitive root of prime using Euler's Totient Theorem
    // reference: https://www.geeksforgeeks.org/primitive-root-of-a-prime-number-n-modulo-n/
    Fint prime_mpz = modulus;
    assert(mpz_probab_prime_p(prime_mpz.get_mpz_t(), 10) && "prime not prime");
    prime_mpz -= 1; // totient
    vector<Fint> factors;
    get_all_factors(prime_mpz, factors);
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

/**
    @param arr a polynomial in coefficient form, returned as map of RoU
    @param logn value where 2 ^ logn + 1 = modulo 
 */
bool is_Fermat_prime(const Fint& modulo) {
    // fermat prime: 2 ^ {2^n} + 1
    // n = 0, 1, 2, 3, 4
    return modulo == 3 || modulo == 5 || modulo == 17 || modulo == 257 || modulo == 65537;
}

u32 log2(u64 len) {
    // 2^logn + 1 = modulo
    // logn = log2(modulo - 1)
    // len = 2^logn
    u32 logn = 0;
    while (len > 1) {
        len >>= 1;
        ++logn;
    }
    return logn;
}

u32 mpz2ul(const Fint& mp) {
    return mpz_get_ui(mp.get_mpz_t());
}

u64 mpz2ull(const Fint& mp) {
    u64 result = 0;
    mpz_export(&result, 0, -1, sizeof(result), 0, 0, mp.get_mpz_t());
    return result;
}

void ul2mpz(Fint& mp, const u32& input) {
    mpz_set_ui(mp.get_mpz_t(), input);
}

void ull2mpz(Fint& mp, const u64& input) {
    mpz_import(mp.get_mpz_t(), 1, -1, sizeof(u64), 0, 0, &input);
}

void enforce_modulus(Fint& mp, const Fint& modulo) {
    if (mp >= modulo)
       mp %= modulo;
    else if (mp < 0) {
       while (mp < 0)
          mp += modulo;
    }
}

inline Fint inverseModular(Fint a, Fint mod)
{
   Fint result;
   int flag = mpz_invert(result.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());

   assert(flag != 0 && "ErrorNoInverse: No inverse found");

   return result;
}

Fint multInverse(const Fint& m, const Fint& prime_mp, const Fint* dlog = nullptr, const Fint* w = nullptr)
{
   if (m == 0) assert(0 && "ErrorNoInverse: No inverse found");
   if (m == 1) {
      return 1;
   }
   if (dlog && w) {
      Fint temp = prime_mp - 1 - dlog[(ptrdiff_t) mpz_get_ui(m.get_mpz_t())];
      return w[(ptrdiff_t)
         mpz_get_ui(
            temp.get_mpz_t()
         )];
   }

   return inverseModular(m, prime_mp);
}