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
    for (size_t i = 0; i < v; ++i)
        new_shape[i] = 1;
    std::ranges::copy(shape.begin(), shape.end(), new_shape.begin() + v);
    // for (size_t i = 0; i < shape.size(); ++i)
    //     new_shape[i + v] = shape[i];
    return {dim + v, new_shape, data};
}

template<typename T>
nd_vector<T> nd_vector<T>::dim_append(size_t v) {
    std::vector<size_t>*  shape_arr = new std::vector<size_t>{shape.size() + v}; // potential memory leak 
    std::span<size_t> new_shape{*shape_arr};
    std::ranges::copy(shape.begin(), shape.end(), new_shape.begin());
    size_t shape_size = shape.size();
    for (size_t i = 0; i < v; ++i)
        new_shape[shape_size + i] = 1;
    return {dim + v, new_shape, data};
}

template class nd_vector<int>;
