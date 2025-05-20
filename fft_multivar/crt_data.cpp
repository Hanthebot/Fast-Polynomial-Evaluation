#include "crt_data.h"
#include <set>
using namespace std;

CRTData::CRTData(nd_vector<Fint>& arr, u32 prime) {
    modulo = prime;
    prime_mp = Fint(prime);
    isPrime = mpz_probab_prime_p(prime_mp.get_mpz_t(), 10);
    coeff = &arr;
    Fint dim_mp;
    ul2mpz(dim_mp, arr.getDim());
    buffer_size = pow(prime_mp, dim_mp);
    findFieldsCRT();
    preprocess();
}

// Recursively count tables
u32 CRTData::getTableCount() const {
    if (!subField) return 1;
    u32 count = 0;
    for (u32 i = 0; i < fields_used.size(); ++i) {
        count += subField[i]->getTableCount();
    }
    return count;
}

// Recursively count table size
Fint CRTData::getTableSize() const {
    if (!subField) return buffer_size * ((log2(modulo) + 7) >> 3);
    Fint table_size = 0;
    for (u32 i = 0; i < fields_used.size(); ++i) {
        table_size += subField[i]->getTableSize();
    }
    return table_size;
}

void CRTData::printTableSize(string st) const {
    cout << st << "Table size: " << getTableSize() << endl;
    cout << st << "Table count: " << getTableCount() << endl;
    for (u32 i = 0; i < fields_used.size(); ++i) {
        subField[i]->printTableSize(st + "\t");
    }
}

// Direct access for leaf node
Fint CRTData::access(const vector<Fint>& eval_point) {
    // Evaluate using FFT tables for this field
    nd_vector<Fint> elem = *preprocessed;
    for (u32 d = 0; d < eval_point.size(); ++d) {
        Fint coord_mod = eval_point[d] % modulo;
        u64 coord_int = mpz2ull(coord_mod);
        elem = elem[coord_int]; // assume dlog, rev permutation done
    }
    // enforce_modulus(elem.get(), modulo);
    return 0; //elem.get();
}

// Recursive evaluation
Fint CRTData::evaluate(const vector<Fint>& eval_point) {
    if (!subField && !preprocessed) {
        return 0; // preprocessing not done yet
    }
    if (!subField) {
        return access(eval_point);
    } else {
        // CRT reconstruction
        Fint result = 0;
        for (u32 i = 0; i < fields_used.size(); ++i) {
            Fint val = subField[i]->evaluate(eval_point);
            result += interpolate[i] * val;
        }
        result %= common_prod;
        return result % modulo;
    }
}

void CRTData::evaluate_stats(Fint& mul_count, Fint& access_count) {
    if (!preprocessed) {
        return; // preprocessing not done yet
    }
    if (!subField) {
        access_count += 1;
        return;
    } else {
        // CRT reconstruction
        for (u32 i = 0; i < fields_used.size(); ++i) {
            subField[i]->evaluate_stats(mul_count, access_count);
        }
        Fint temp;
        ul2mpz(temp, fields_used.size());
        mul_count += temp;
        return;
    }
}

void CRTData::checkSpaceOptimum(const span<size_t>& shape) {
    // update the fields_used
    if (fields_used.size() <= 1) return;
    Fint prime_pow_acc = 0;
    u32 m = coeff->getDim();
    for (const u32& field : fields_used) {
        prime_pow_acc += pow((Fint) field, (Fint) m) * ((log2(field) + 7) >> 3);
    }
    if (prime_pow_acc > buffer_size * ((log2(modulo) + 7) >> 3)) {
        // inefficient in space
        fields_used.clear();
    }
}

void CRTData::findFieldsCRT() {
    // find fields to be used
    fields_used.clear();
    set<Fint> val_set;
    vector<bool> isFermat;
    u32 max_fft_field;
    vector<u32> deg_vec(coeff->getDim(), 0);
    for (u32 i = 0; i < coeff->getDim(); ++i)
        deg_vec[i] = coeff->getShape()[i] - 1;

    compute_M(M, deg_vec, modulo);
    // largest field used: guaranteed to be less than modulo, if modulo prime
    u32 signal = find_fields(modulo, M, fields_used, isFermat, max_fft_field, val_set);
    assert(signal == 0 && "Error: from finding fields");
    checkSpaceOptimum(coeff->getShape());
    // check if fields are being efficient in terms of size
    if (isPrime && fields_used.size() == 1) {
        // direct evaluation
        fields_used.clear();
    }
}

void CRTData::preprocess() {
    // create a new nd_vector
    size_t m = coeff->getDim();
    shape_vec.clear();
    shape_vec.resize(m, modulo);
    unit_vec.clear();
    unit_vec.resize(m, 1);
    for (int i = m - 2; i >= 0; --i) {
        unit_vec[i] = unit_vec[i + 1] * shape_vec[i + 1];
    }
    span<size_t> shape = shape_vec;
    span<size_t> units = unit_vec;

    preprocessed = new nd_vector<Fint>(m, shape, units, nullptr, (size_t) mpz2ull(buffer_size));

    // i.e. direct computation
    if (fields_used.empty()) return;
    FD = new FieldData[fields_used.size()]; // field data
    // set up FD
    for (size_t i = 0; i < fields_used.size(); ++i) {
        FD[i].prime = fields_used[i];
        ul2mpz(FD[i].prime_mp, fields_used[i]);
    }
    // prepare for interpolation, if CRT computation
    common_prod = 1;
    for (u32 i = 0; i < fields_used.size(); ++i) {
        common_prod *= FD[i].prime_mp;
    }
    interpolate.resize(fields_used.size());
    for (u32 i = 0; i < fields_used.size(); ++i) {
        Fint temp = common_prod / FD[i].prime_mp;
        Fint inv = 1; // multInverse(temp % subField[i].prime_mp, subField[i].prime_mp, subField[i].dlog.data(), subField[i].w.data());
        interpolate[i] = mpz2ul(inv * temp);
    }
    
    subField = new CRTData*[fields_used.size()];
    for (size_t i = 0; i < fields_used.size(); ++i) {
        subField[i] = new CRTData(*preprocessed, fields_used[i]);
    }
}