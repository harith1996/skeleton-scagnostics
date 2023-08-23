#ifndef MULTIVAL_HPP
#define MULTIVAL_HPP

#include <array>
#include <type_traits>
#include <algorithm>

template<class T, size_t N>
struct MultiVal {
    std::array<T, N> data;

    template<class S, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
    MultiVal(MultiVal<S, N> other) {
        (*this) = other;
    }

    template<class... A>
    MultiVal(A... args) : data { args... } {
    }

    template<class S, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
    MultiVal(S s) {
        (*this) = s;
    }

    template<class S, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
    MultiVal operator=(S s) {
        for(int i = 0; i < N; i++) {
            data[i] = static_cast<T>(s);
        }
        return *this;
    }

    template<class S, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
    MultiVal operator=(const MultiVal<S, N> other) {
        for(int i = 0; i < N; i++) {
            data[i] = static_cast<T>(other.data[i]);
        }
        return *this;
    }

    MultiVal operator+=(const MultiVal& other) {
        for(int i = 0; i < N; i++) {
            data[i] += other.data[i];
        }
        return *this;
    }

    MultiVal operator-=(const MultiVal& other) {
        for(int i = 0; i < N; i++) {
            data[i] -= other.data[i];
        }
        return *this;
    }

    MultiVal operator*=(const MultiVal& other) {
        for(int i = 0; i < N; i++) {
            data[i] *= other.data[i];
        }
        return *this;
    }

    MultiVal operator/=(const MultiVal& other) {
        for(int i = 0; i < N; i++) {
            data[i] /= other.data[i];
        }
        return *this;
    }

    MultiVal operator+(const MultiVal& other) const {
        MultiVal copy = *this;
        copy += other;
        return copy;
    }

    MultiVal operator-(const MultiVal& other) const {
        MultiVal copy = *this;
        copy -= other;
        return copy;
    }

    MultiVal operator*(const MultiVal& other) const {
        MultiVal copy = *this;
        copy *= other;
        return copy;
    }

    MultiVal operator/(const MultiVal& other) const {
        MultiVal copy = *this;
        copy /= other;
        return copy;
    }

    T& operator[](size_t index) {
        return data[index];
    }

    const T& operator[](size_t index) const {
        return data[index];
    }
};

template<class T, class S, size_t N, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
MultiVal<T,N> operator+(S s, MultiVal<T, N> m) {
    return MultiVal<T, N>(s) + m;
}

template<class T, class S, size_t N, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
MultiVal<T,N> operator-(S s, MultiVal<T, N> m) {
    return MultiVal<T, N>(s) - m;
}

template<class T, class S, size_t N, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
MultiVal<T,N> operator*(S s, MultiVal<T, N> m) {
    return MultiVal<T, N>(s) * m;
}

template<class T, class S, size_t N, class = typename std::enable_if<std::is_convertible<S, T>::value>::type>
MultiVal<T,N> operator/(S s, MultiVal<T, N> m) {
    return MultiVal<T, N>(s) / m;
}

template<class T, size_t N>
bool operator==(MultiVal<T, N> a, MultiVal<T, N> b) {
    for(int i = 0; i < N; i++) {
        if(a[i] != b[i]) return false;
    }
    return true;
}

template<class T, size_t N>
bool operator!=(MultiVal<T, N> a, MultiVal<T, N> b) {
    return !(a == b);
}

template<class T, size_t N>
MultiVal<T, N> min(MultiVal<T, N> a, MultiVal<T, N> b) {
    using std::min;

    MultiVal<T, N> result;
    for(int i = 0; i < N; i++) {
        result[i] = min(a[i], b[i]);
    }
    return result;
}

template<class T, size_t N>
MultiVal<T, N> max(MultiVal<T, N> a, MultiVal<T, N> b) {
    using std::max;

    MultiVal<T, N> result;
    for(int i = 0; i < N; i++) {
        result[i] = max(a[i], b[i]);
    }
    return result;
}

template<class T, size_t N>
T vectorLengthSquared(MultiVal<T, N> v) {
    T result = 0;

    for(int i = 0; i < N; i++) {
        result += v[i] * v[i];
    }

    return result;
}

using Float2 = MultiVal<float, 2>;
using Int2 = MultiVal<int, 2>;

#endif
