#include "nd_vector.h"
#include <algorithm>
#include <vector>

template<typename T>
nd_vector<T>::nd_vector(const size_t& dim, const std::span<size_t>& shape, const std::span<T>& data) : dim{dim}, shape{shape}, data{data} {
    assert(shape.size() == dim && "dimension and shape does not match");
    u32 length = 1;
    for (const int& v : shape)
        length *= v;
    assert(length == data.size() && "size does not match");
}

template<typename T>
nd_vector<T> nd_vector<T>::operator[](size_t i) {
    assert(dim >= 1 && "no inner dimension");
    assert(i < shape[0] && "index out of bounds");
    size_t unit = data.size() / shape[0];
    auto beg_iter = data.begin() + i * unit;
    return {dim - 1, 
        {next(shape.begin()), shape.end()}, 
        {beg_iter, beg_iter + unit}
    };
}

template<typename T>
void nd_vector<T>::set_range(const nd_vector<T>& vec) {
    // set range if latter dimension are the same, or vec[dim] being 1
    if (this->dim == vec.dim) {
        if (this == &vec) return;
        assert(isEqual(this->shape, vec.shape) && "shape does not match");
        copy(vec.data.begin(), vec.data.end(), this->data.begin());
        return;
    }
    assert(this->dim > vec.dim && "dimension does not match");
    for (int i = 0; i < this->shape[0]; ++i) {
        (*this)[i].set_range(vec);
    }
}

template<typename T>
void nd_vector<T>::set_val(const T& val, size_t i) {
    assert(i < data.size() && "index out of bounds");
    data[i] = val;
}

template<typename T>
nd_vector<T> nd_vector<T>::dim_prepend(size_t v) {
    std::vector<size_t>* const shape_arr = new std::vector<size_t>{shape.size() + v}; // potential memory leak 
    std::span<size_t> new_shape{*shape_arr};
    for (int i = 0; i < v; ++i)
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
    for (int i = 0; i < v; ++i)
        new_shape[shape_size + i] = 1;
    return {dim + v, new_shape, data};
}

template<typename T>
bool isEqual(const std::span<T>& a, const std::span<T>& b) {
    bool rst = a.size() == b.size();
    for (int i = 0; i < a.size(); ++i)
        rst &= a[i] == b[i];
    return rst;
}

template class nd_vector<int>;
