#ifndef STAT_HPP
#define STAT_HPP

#include <array>
#include <algorithm>
#include <vector>

template<class T>
struct Stat {
    T average = 0;
    // sigma^2 * count
    T sigmaSquaredN = 0;
    T minValue = 0, maxValue = 0;
    size_t count = 0;
    std::vector<T> values;
    void put(T value) {
        put(1, value);
    }

    void put(size_t n, T value) {
        // http://datagenetics.com/blog/november22017/index.html
        // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
        if(n == 0) return;

        //addes n in the count
        for(int k = 0; k < n; k++) values.push_back(value);

        count += n;

        if(count == n) {
            average = value;
            sigmaSquaredN = 0;
            minValue = value;
            maxValue = value;
            return;
        }

        auto oldAverage = average;
        average += (value - average) * (n * 1.0 / count);
        sigmaSquaredN += n * (value - oldAverage) * (value - average);

        //merge std::min,max into overload set
        using std::min;
        using std::max;

        //Global/Koenig-Lookup should make this work
        minValue = min(minValue, value);
        maxValue = max(maxValue, value);
    }

    void setNewRange(T newMin, T newMax){
        minValue = newMin;
        maxValue = newMax;
    }

    bool empty() const {
        return count == 0;
    }

    T variance() const {
        return sigmaSquaredN / count;
    }
};

template<class T>
Stat<T> operator*(int x, Stat<T> stat) {
    stat.count *= x;
    stat.sigmaSquaredN *= x;
    return stat;
}

template<class T>
Stat<T> operator+(Stat<T> a, Stat<T> b) {
    Stat<T> result;
    if(b.count == 0) return a;
    if(a.count == 0) return b;

    result.count = a.count + b.count;
    result.average = (a.average * a.count + b.average * b.count) / result.count;

    auto delta = a.average - b.average;
    result.sigmaSquaredN = a.sigmaSquaredN + b.sigmaSquaredN + delta * delta * a.count * b.count / result.count;

    using std::min;
    using std::max;

    result.minValue = min(a.minValue, b.minValue);
    result.maxValue = max(a.maxValue, b.maxValue);

    return result;
}

#endif
