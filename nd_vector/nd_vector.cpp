#include "nd_vector.h"
#include <algorithm>
#include <vector>

bool isEqual(const std::span<size_t>& a, const std::span<size_t>& b) {
    bool rst = a.size() == b.size();
    for (size_t i = 0; i < a.size(); ++i)
        rst &= (a[i] == b[i]);
    return rst;
}

template<typename T>
void nd_vector<T>::set_val(const T& val, size_t i) const {
    assert(i < data.size() && "index out of bounds");
    data[i] = val;
}

template<typename T>
nd_vector<T> nd_vector<T>::dim_prepend(size_t v) {
    std::vector<size_t>* const shape_arr = new std::vector<size_t>{shape.size() + v}; // potential memory leak 
    std::span<size_t> new_shape{*shape_arr};
    std::vector<size_t>* const unit_arr = new std::vector<size_t>{unit.size() + v}; // potential memory leak 
    std::span<size_t> new_unit{*unit_arr};
    std::ranges::copy(shape.begin(), shape.end(), new_shape.begin() + v);
    std::ranges::copy(unit.begin(), unit.end(), new_unit.begin() + v);
    new_unit[v - 1] = unit[0] * shape[0];
    for (size_t i = 0; i < v; ++i) {
        new_shape[i] = 1;
        new_unit[i] = new_unit[v - 1];
    }
    // for (size_t i = 0; i < shape.size(); ++i)
    //     new_shape[i + v] = shape[i];
    return {dim + v, new_shape, new_unit, data};
}

template<typename T>
nd_vector<T> nd_vector<T>::dim_append(size_t v) {
    std::vector<size_t>*  shape_arr = new std::vector<size_t>{shape.size() + v}; // potential memory leak 
    std::span<size_t> new_shape{*shape_arr};
    std::ranges::copy(shape.begin(), shape.end(), new_shape.begin());
    std::vector<size_t>*  unit_arr = new std::vector<size_t>{unit.size() + v}; // potential memory leak 
    std::span<size_t> new_unit{*unit_arr};
    std::ranges::copy(unit.begin(), unit.end(), new_unit.begin());
    size_t shape_size = shape.size();
    for (size_t i = 0; i < v; ++i) {
        new_shape[shape_size + i] = 1;
        new_unit[shape_size + i] = 1;
    }
    return {dim + v, new_shape, new_unit, data};
}

template class nd_vector<int>;
