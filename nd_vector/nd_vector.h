#ifndef ND_VECTOR_H
#define ND_VECTOR_H

#include <iostream>
#include <span>
#include <vector>
#include <cassert>

// using namespace std;

typedef unsigned int u32;

bool isEqual(const std::span<size_t>& a, const std::span<size_t>& b);

template<typename T>
class nd_vector {
    private:
        size_t dim;
        std::span<size_t> shape;
        std::span<size_t> unit;
        T* data;
        size_t size_;
    public:
        /**
        * @brief Constructs an nd_vector object.
        * 
        * @param dim The dimension of the vector.
        * @param shape The shape of the vector.
        * @param data The data of the vector.
        */
        nd_vector(const size_t& dim, const std::span<size_t>& shape, const std::span<size_t>& unit, T* const data, const size_t& size_);

        /**
        * @brief Accesses the sub-vector at the given index.
        * 
        * @param i The index to access.
        * @return nd_vector<T> The sub-vector at the given index.
        */
        nd_vector<T> operator[](size_t i) const;

        /**
        * @brief Sets the value at the given index.
        * 
        * @param val The value to set.
        * @param i The index at which to set the value.
        */
        void set_val(const T& val, size_t i) const;

        /**
        * @brief Sets a range of values starting at the given index.
        * 
        * @param vec The vector containing the values to set.
        */
        void set_range(const nd_vector<T>& vec) const;

        /**
        * @brief Gets the value at the given index.
        * 
        * @param i The index to access.
        * @return const T& The value at the given index.
        */
        const T& get(size_t i = 0) const;
        T& get(size_t i = 0);

        /**
        * @brief Returns a span of the data.
        * 
        * @return const std::span<T> A span of the data.
        */
        T* const& ptr() const;
        T*& ptr();

        /**
        * @brief Returns a span of the shape.
        * 
        * @return const std::span<size_t> A span of the shape.
        */
        const std::span<size_t> getShape() const;
        const std::span<size_t> getUnit() const;

        /**
        * @brief Returns the size of the vector.
        * 
        * @return size_t The size of the vector.
        */
        size_t size() const;

        /**
        * @brief Returns the dimension of the vector.
        * 
        * @return size_t The dimension of the vector.
        */
        size_t getDim() const;

        /**
        * @brief Prepends a dimension to the vector. For potential future use.
        * 
        * @param v The dimension to prepend.
        * @return nd_vector<T> The resulting vector with the prepended dimension.
        */
        nd_vector<T> dim_prepend(size_t v);

        /**
        * @brief Appends a dimension to the vector. For potential future use.
        * 
        * @param v The dimension to append.
        * @return nd_vector<T> The resulting vector with the appended dimension.
        */
        nd_vector<T> dim_append(size_t v);

        void reshape(const std::span<size_t>& new_shape) {
            u32 new_size = 1;
            for (const auto& v : new_shape)
                new_size *= v;
            assert((size_ == new_size) && "shape does not match");
            shape = new_shape;
            dim = new_shape.size();
            std::vector<size_t>* new_unit = new std::vector<size_t>(dim, 1);
            unit = std::span<size_t>(*new_unit);
            unit[0] = size_ / new_shape[0];
            for (size_t i = 1; i < dim; ++i) {
                unit[i] = unit[i - 1] / new_shape[i];
            }
        }
        
        /**
        * @brief Sets the values of the vector to the product of the values of the given vectors.
        *
        * @param a The first vector.
        * @param b The second vector. 
        */
        void set_mul(const nd_vector<T>& a, const nd_vector<T>& b) const {
            assert(a.size() == b.size() && "size does not match");
            assert(size_ == a.size() && "size does not match");
            for (size_t i = 0; i < a.size(); ++i)
                data[i] = a.data[i] * b.data[i];
        }
        
        /**
        * @brief Sets the values of the vector to the product of the values of the given vector and scalar.
        *
        * @param a The first vector.
        * @param b The scalar. 
        */
        void set_mul(const nd_vector<T>& a, const T& b) const {
            assert(size_ == a.size() && "size does not match");
            for (size_t i = 0; i < size_; ++i)
                data[i] = a.data[i] * b;
        }

        /**
        * @brief Sets the values of the vector to the sum of the values of the given vectors.
        *
        * @param a The first vector.
        * @param b The second vector. 
        */
        void set_add(const nd_vector<T>& a, const nd_vector<T>& b) const {
            assert(a.size() == b.size() && "size does not match");
            assert(size_ == a.size() && "size does not match");
            for (size_t i = 0; i < size_; ++i)
                data[i] = a.data[i] + b.data[i];
        }

        /**
        * @brief Sets the values of the vector to the difference of the values of the given vectors.
        *
        * @param a The first vector.
        * @param b The second vector. 
        */
        void set_sub(const nd_vector<T>& a, const nd_vector<T>& b) const {
            assert(a.size() == b.size() && "size does not match");
            assert(size_ == a.size() && "size does not match");
            for (size_t i = 0; i < size_; ++i)
                data[i] = a.data[i] - b.data[i];
        }
        
        /**
        * @brief Swaps the values of the vector with the given vector.
        *
        * @param vec The vector to swap with.
        */
        void swap(const nd_vector<T>& vec) const {
            assert(this->dim == vec.dim && "dimension does not match");
            assert(isEqual(this->shape, vec.shape) && "shape does not match");
            for (size_t i = 0; i < size_; ++i) {
                T temp = data[i];
                data[i] = vec.data[i];
                vec.data[i] = temp;
            }
        }

        /**
        * @brief Outputs the vector to the given output stream.
        * 
        * @tparam T_ The type of the elements in the vector.
        * @param os The output stream.
        * @param vec The vector to output.
        * @return std::ostream& The output stream.
        */
        template<typename T_>
        friend std::ostream& operator<<(std::ostream& os, const nd_vector<T_>& vec);
        template<typename T_>
        friend std::ostream& print(std::ostream& os, const nd_vector<T_>& vec, const std::string& prefix);
        template<typename T_>
        friend std::istream& operator>>(std::istream& is, nd_vector<T_>& vec);
};

template<typename T>
nd_vector<T>::nd_vector(const size_t& dim, const std::span<size_t>& shape, const std::span<size_t>& unit, T* const data, const size_t& size_) : dim{dim}, shape{shape}, unit{unit}, data{data}, size_{size_} {
    assert(shape.size() == dim && "dimension and shape does not match");
    assert(unit.size() == dim && "dimension and unit does not match");
}

template<typename T>
nd_vector<T> nd_vector<T>::operator[](size_t i) const {
    assert(dim >= 1 && "no inner dimension");
    assert(i < shape[0] && "index out of bounds");
    auto beg_ptr = data + i * unit[0];
    return {dim - 1, 
        shape.subspan(1), 
        unit.subspan(1), 
        beg_ptr,
        unit[0]
    };
}

template<typename T>
inline const T& nd_vector<T>::get(size_t i) const {
    assert((i == 0 || i < size_) && "index out of bounds");
    return data[i];
}

template<typename T>
inline T& nd_vector<T>::get(size_t i) {
    assert((i == 0 || i < size_) && "index out of bounds");
    return data[i];
}

template<typename T>
inline T* const & nd_vector<T>::ptr() const {
    return data;
}

template<typename T>
inline T*& nd_vector<T>::ptr() {
    return data;
}

template<typename T>
inline const std::span<size_t> nd_vector<T>::getShape() const {
    return shape;
}

template<typename T>
inline const std::span<size_t> nd_vector<T>::getUnit() const {
    return unit;
}

template<typename T>
inline size_t nd_vector<T>::size() const {
    return size_;
}

template<typename T>
inline size_t nd_vector<T>::getDim() const {
    return this->dim;
}

template<typename T>
std::ostream& print(std::ostream& os, const nd_vector<T>& vec, const std::string& prefix = "") {
    if (vec.dim == 0) {
        os << vec.get() << "\n";
        return os;
    }
    if (vec.dim == 1) {
        os << prefix <<  "[ ";
        for (size_t i = 0; i < vec.shape[0]; ++i) {
            os << vec.get(i) << " ";
        }
        os << "]\n";
        return os;
    } else {
        os << prefix << "[ \n";
        for (size_t i = 0; i < vec.shape[0]; ++i) {
            print(os, vec[i], "  " + prefix);
        }
        os << prefix << "]\n";
        return os;
    }
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const nd_vector<T>& vec) {
    os << "dim: " << vec.dim << "\n";
    os << "shape: [ ";
    for (const auto& v : vec.shape)
        os << v << " ";
    os << "]\n";
    os << "data: \n";
    print(os, vec);
    return os;
}

template<typename T>
std::istream& operator>>(std::istream& is, nd_vector<T>& vec) {
    for (size_t i = 0; i < vec.size_; ++i)
        is >> vec.data[i];
    return is;
}

template<typename T>
void nd_vector<T>::set_range(const nd_vector<T>& vec) const {
    // set range if latter dimension are the same, or vec[dim] being 1
    if (this->dim == vec.dim) {
        if (this == &vec) return;
        assert(isEqual(this->shape, vec.shape) && "shape does not match");
        for (size_t i = 0; i < size_; ++i) {
            data[i] = vec.data[i];
        }
        return;
    }
    assert(this->dim > vec.dim && "dimension does not match");
    for (size_t i = 0; i < this->shape[0]; ++i) {
        (*this)[i].set_range(vec);
    }
}

#endif // ND_VECTOR_H