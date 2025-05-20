#include "fft_multivar.h"

Fint evaluate_nd(const nd_vector<Fint>& arr, const span<const Fint>& points, Fint& mul_counter) {
    if (arr.getDim() == 0) {
        return 1; // arr.get();
    }
    Fint result = 0;
    for (size_t i = arr.getShape()[0] - 1; ; --i) {
        result *= points[0];
        result += evaluate_nd(arr[i], points.subspan(1), mul_counter);
        if (i == 0) break;
    }
    return result;
}

Fint evaluate_brutal(const nd_vector<Fint>& arr, const vector<Fint>& points, Fint& mul_counter) {
    span<const Fint> points_span = points;
    Fint result = evaluate_nd(arr, points_span, mul_counter);
    Fint acc_count = 1, temp;
    for (size_t i = 0; i < arr.getDim(); ++i) {
        mpz_mul_ui(acc_count.get_mpz_t(), acc_count.get_mpz_t(), arr.getShape()[i]);
        mul_counter += acc_count;
    }
    
    return result;
}