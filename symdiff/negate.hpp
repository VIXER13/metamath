#ifndef SYMDIFF_NEGATE_HPP
#define SYMDIFF_NEGATE_HPP

#include "constant.hpp"

namespace metamath::symdiff {

template<class E>
class negate;

template<class U>
using negate_type = std::conditional_t<
    is_integral_constant<U>{},
    integral_constant<integral_constant_t<U>, -integral_constant_v<U>{}>,
    std::conditional_t<is_constant<U>{}, U, negate<U>>
>;

template<class E>
class negate : public expression<negate<E>> {
    const E e;

public:
    template<uintmax_t M>
    using derivative_type = negate_type<typename E::template derivative_type<M>>;

    constexpr negate(const expression<E>& e) :
        e{e()} {}

    template<class T>
    constexpr auto operator()(const T& x) const -> decltype(-e(x)) {
        return -e(x);
    }

    template<uintmax_t M>
    constexpr derivative_type<M> derivative() const {
        return -e.template derivative<M>();
    }
};

template<class T, T N>
constexpr integral_constant<T, -N> operator-(const integral_constant<T, N>&) {
    return integral_constant<T, -N>{};
}

template<class T>
constexpr constant<T> operator-(const constant<T>& e) {
    return constant<T>{-e.value};
}

template<class E>
constexpr negate<E> operator-(const expression<E>& e) {
    return negate<E>{e};
}

}

#endif