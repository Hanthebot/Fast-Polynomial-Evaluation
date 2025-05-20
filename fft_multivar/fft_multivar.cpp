#include <iostream>
#include <fstream>
#include <cmath>
#include "fft_multivar.h"
#include "util.h"
#include "crt_data.h"

using namespace std;

int main(int argc, char* argv[]) {
    vector<u32> degs_vec;
    vector<vector<Fint>> evaluation_points;
    vector<Fint> results;
    Fint mul_counter = 0, modulo = 0, capital_M;

    // initializing field and consts
    modulo = 4093;
    std::ostream& my_cout = std::cout;
    for (u32 m = 3; m <= 8; ++m) {
        u32 add_ind = 0;
        u32 eg_count = 1;
        degs_vec.clear();
        degs_vec = vector<u32>(m, 1);
        u64 total_len = 1;
        for (size_t i = 0; i < m; ++i) {
            total_len *= degs_vec[i] + 1;
        }
        my_cout << "===m: " << m << "===" << endl;
        while (total_len < 5000000) {
            if (total_len >= 10000) {
                eg_count += 1;
                // Forward output after this line to a file stream
                string output_name = "multivar_sample/output_";
                for (size_t i = 0; i < m; ++i) {
                    output_name += to_string(degs_vec[i]);
                    if (i != m - 1) {
                        output_name += "_";
                    }
                }
                output_name += ".txt";
                ofstream outfile(output_name);
                compute_M(capital_M, degs_vec, modulo);
                init_evaluation_points(evaluation_points, m, modulo, outfile);
    
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
                nd_vector<Fint> coeff {m, shape, units, nullptr, total_len};
    
                coeff_init(coeff, m, degs_vec, modulo, outfile);
                
                outfile << "Terms: " << total_len << endl;
                CRTData CRT(coeff, mpz2ul(modulo));
                outfile << "Table count: " << CRT.getTableCount() << endl
                << "Table size: " << CRT.getTableSize() << endl;

                Fint mul_count = 0, access_count = 0;
                CRT.evaluate_stats(mul_count, access_count);
                outfile << "Multiplication count: " << mul_count << endl;
                outfile << "Access count: " << access_count << endl;
                
                auto start = chrono::high_resolution_clock::now();
                for (const auto& evaluation_point : evaluation_points) {
                    CRT.evaluate(evaluation_point);
                }
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
                mpf_class time_spent_fft = from_micro(duration.count());
                outfile << "Time taken for FFT: " << (time_spent_fft / (double) evaluation_points.size()) << " s" << endl;
    
                start = chrono::high_resolution_clock::now();
                for (size_t i = 0; i < evaluation_points.size(); ++i) {
                    evaluate_brutal(coeff, evaluation_points[i], mul_counter) % modulo;
                }
                stop = chrono::high_resolution_clock::now();
                long long int naive_time = chrono::duration_cast<chrono::microseconds>(stop - start).count();
                mpf_class time_spent = from_micro(naive_time);
                if (evaluation_points.size() > 0) {
                    time_spent /= (double) evaluation_points.size();
                }
                outfile << "Time taken for naive: " << time_spent << " s" << endl;
                outfile.close();
                my_cout << "eg_count: " << eg_count << "\r";
            }
            ++degs_vec[add_ind];
            add_ind = (add_ind + 1) % m;
            total_len = 1;
            for (size_t i = 0; i < m; ++i) {
                total_len *= degs_vec[i] + 1;
            }
        }
        my_cout << endl;
    }
    return 0;
}