#ifndef ND_VECTOR_H
#define ND_VECTOR_H

#include <iostream>
#include <span>
#include <cassert>

// using namespace std;

typedef unsigned int u32;

template<typename T>
class nd_vector {
    private:
        size_t dim;
        std::span<size_t> shape;
        std::span<T> data;
    public:
        /**
        * @brief Constructs an nd_vector object.
        * 
        * @param dim The dimension of the vector.
        * @param shape The shape of the vector.
        * @param data The data of the vector.
        */
        nd_vector(const size_t& dim, const std::span<size_t>& shape, const std::span<T>& data);

        /**
        * @brief Accesses the sub-vector at the given index.
        * 
        * @param i The index to access.
        * @return nd_vector<T> The sub-vector at the given index.
        */
        nd_vector<T> operator[](size_t i);

        /**
        * @brief Sets the value at the given index.
        * 
        * @param val The value to set.
        * @param i The index at which to set the value.
        */
        void set_val(const T& val, size_t i);

        /**
        * @brief Sets a range of values starting at the given index.
        * 
        * @param vec The vector containing the values to set.
        */
        void set_range(const nd_vector<T>& vec);

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
        const std::span<T> span() const;
        std::span<T> span();

        /**
        * @brief Returns the size of the vector.
        * 
        * @return size_t The size of the vector.
        */
        size_t size() const;

        /**
        * @brief Prepends a dimension to the vector.
        * 
        * @param v The dimension to prepend.
        * @return nd_vector<T> The resulting vector with the prepended dimension.
        */
        nd_vector<T> dim_prepend(size_t v);

        /**
        * @brief Appends a dimension to the vector.
        * 
        * @param v The dimension to append.
        * @return nd_vector<T> The resulting vector with the appended dimension.
        */
        nd_vector<T> dim_append(size_t v);
        
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
};

template<typename T>
inline const T& nd_vector<T>::get(size_t i) const {
    assert(i < data.size() && "index out of bounds");
    return data[i];
}

template<typename T>
inline T& nd_vector<T>::get(size_t i) {
    assert(i < data.size() && "index out of bounds");
    return data[i];
}

template<typename T>
inline const std::span<T> nd_vector<T>::span() const {
    return data;
}

template<typename T>
inline std::span<T> nd_vector<T>::span() {
    return data;
}

template<typename T>
inline size_t nd_vector<T>::size() const {
    return data.size();
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const nd_vector<T>& vec) {
    os << "dim: " << vec.dim << "\n";
    os << "shape: [ ";
    for (const size_t& s : vec.shape)
        os << s << " ";
    os << "]\ndata: ";
    os << "[ ";
    for (const T& v : vec.data)
        os << v << " ";
    os << "]\n";
    return os;
}

template<typename T>
bool isEqual(const std::span<T>& a, const std::span<T>& b);

#endif // ND_VECTOR_H