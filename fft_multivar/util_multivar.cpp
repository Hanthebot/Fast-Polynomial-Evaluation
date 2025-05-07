#include "util.h"
#include <iostream>
#include <random>

using namespace std;
using namespace std::chrono;

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

void rou_init(vector<Fint>& w, vector<Fint>& dlog, const u32& prime, const Fint& prime_mp, const Fint& rou) {
    w.resize(prime - 1);
    dlog.resize(prime);
    w[0] = 1;
    dlog[mpz2ul(w[0])] = 0;
    w[1] = rou;
    dlog[mpz2ul(w[1])] = 1;
    Fint count = 2;
    for (u32 i = 2; i < prime - 1; ++i) {
        w[i] = w[i - 1] * w[1];
        enforce_modulus(w[i], prime_mp);
        dlog[mpz2ul(w[i])] = count;
        count++;
    }
    // mul_counter += len-2;
}

void recur_input(const nd_vector<Fint>& coeff, const span<const u32>& degs, const Fint& modulo) {
    if (degs.size() == 0) return;
    if (degs.size() == 1) {
        Fint temp;
        for (u32 i = 0; i <= degs[0]; ++i) {
            // in case degree < len, not using operator>> on nd_vector directly, but get()
            cin >> temp;
            enforce_modulus(temp, modulo);
            *(coeff.ptr() + i) = temp;
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

void print_stats(Fint& mul_counter, microseconds duration) {
    cout << "\n===Result===\n";
    mpf_class speed = ill2mpf(duration.count());
    mpf_class counter = ill2mpf(mul_counter.get_ui());
    if (counter == 0) counter = 1;
    speed /= counter;
    
    cout << "It took " << from_micro(duration.count()) << " s" << endl;
    cout << "Multiplication counter: " << mul_counter << endl;
    cout << "Speed: " << speed << " s/M mult" << endl;
}

void print_dlog(const nd_vector<Fint>& coeff, Fint* dlog, const u32& prime, const string& prefix) {
    if (coeff.getDim() == 0) {
        cout << prefix << coeff.get() << "\n";
        return;
    }
    if (coeff.getDim() == 1) {
        for (u32 counter = 1; counter < prime; ++counter) {
            // assert(dlog.find(i) != dlog.end() && "i not found");
            cout << prefix << counter << ") : " << coeff[mpz2ul(dlog[counter])].get() << endl;
        }
        cout << endl;
        return;
    } else {
        for (u32 counter = 1; counter < prime; ++counter) {
            print_dlog(coeff[mpz2ul(dlog[counter])], dlog, prime, prefix + to_string(counter) + ", ");
        }
        cout << endl;
        return;
    }
}