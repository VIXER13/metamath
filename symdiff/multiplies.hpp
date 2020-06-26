#ifndef SYMDIFF_MULTIPLIES_HPP
#define SYMDIFF_MULTIPLIES_HPP

#include "plus.hpp"

namespace metamath::symdiff {

template<class E1, class E2>
class multiplies;

template<class U1, class U2>
using multiplies_type = std::conditional_t<
    is_integral_constant<U1>{} && is_integral_constant<U2>{},
    integral_constant<decltype(integral_constant_v<U1>{} * integral_constant_v<U2>{}), integral_constant_v<U1>{} * integral_constant_v<U2>{}>,
    std::conditional_t<
        (is_integral_constant<U1>{} && !integral_constant_v<U1>{}) || (is_integral_constant<U2>{} && integral_constant_v<U2>{} == 1),
        U1,
        std::conditional_t<
            (is_integral_constant<U2>{} && !integral_constant_v<U2>{}) || (is_integral_constant<U1>{} && integral_constant_v<U1>{} == 1),
            U2,
            std::conditional_t<
                is_constant<U1>{} && is_constant<U2>{},
                constant<decltype(constant_t<U1>{} * constant_t<U2>{})>,
                multiplies<U1, U2>
            >
        >
    >
>;

template<class E1, class E2>
class multiplies : public expression<multiplies<E1, E2>> {
    const E1 e1;
    const E2 e2;

public:
    template<uintmax_t M>
    using derivative_type = plus_type<
        multiplies_type<typename E1::template derivative_type<M>, E2>,
        multiplies_type<E1, typename E2::template derivative_type<M>>
    >;

    constexpr multiplies(const expression<E1>& e1, const expression<E2>& e2) :
        e1{e1()}, e2{e2()} {}

    template<class T>
    constexpr auto operator()(const T& x) const -> decltype(e1(x) * e2(x)) {
        return e1(x) * e2(x);
    }

    template<uintmax_t M>
    constexpr derivative_type<M> derivative() const {
        return e1.template derivative<M>() * e2 + e1 * e2.template derivative<M>();
    }
};

template<class T1, T1 N1, class T2, T2 N2>
constexpr integral_constant<decltype(N1 * N2), N1 * N2> operator*(const integral_constant<T1, N1>&, const integral_constant<T2, N2>&) {
    return integral_constant<decltype(N1 * N2), N1 * N2>{};
}

template<class T1, class E2>
constexpr integral_constant<T1, 0> operator*(const integral_constant<T1, 0>&, const expression<E2>&) {
    return integral_constant<T1, 0>{};
}

template<class E1, class T2>
constexpr integral_constant<T2, 0> operator*(const expression<E1>&, const integral_constant<T2, 0>&) {
    return integral_constant<T2, 0>{};
}

template<class T1, class E2>
constexpr const E2& operator*(const integral_constant<T1, 1>&, const expression<E2>& e2) {
    return e2();
}

template<class E1, class T2>
constexpr const E1& operator*(const expression<E1>& e1, const integral_constant<T2, 1>&) {
    return e1();
}

template<class T1, class T2>
constexpr constant<decltype(T1{} * T2{})> operator*(const constant<T1>& e1, const constant<T2>& e2) {
    return constant<decltype(T1{} * T2{})>{e1.value * e2.value};
}

template<class T1, T1 N1, class T2>
constexpr std::enable_if_t<N1 && N1 != 1, constant<decltype(T1{} * T2{})>> operator*(const integral_constant<T1, N1>&, const constant<T2>& e2) {
    return constant<decltype(T1{} * T2{})>{N1 * e2.value};
}

template<class T1, class T2, T2 N2>
constexpr std::enable_if_t<N2 && N2 != 1, constant<decltype(T1{} * T2{})>> operator*(const constant<T1>& e1, const integral_constant<T2, N2>&) {
    return constant<decltype(T1{} * T2{})>{e1.value * N2};
}

template<class E1, class E2>
constexpr multiplies<E1, E2> operator*(const expression<E1>& e1, const expression<E2>& e2) {
    return multiplies<E1, E2>{e1, e2};
}

template<class T1, class E2>
constexpr std::enable_if_t<std::is_arithmetic_v<T1>, multiplies_type<constant<T1>, E2>> operator*(const T1& e1, const expression<E2>& e2) {
    return constant<T1>{e1} * e2();
}

template<class E1, class T2>
constexpr std::enable_if_t<std::is_arithmetic_v<T2>, multiplies_type<E1, constant<T2>>> operator*(const expression<E1>& e1, const T2& e2) {
    return e1() * constant<T2>{e2};
}

}

#endif