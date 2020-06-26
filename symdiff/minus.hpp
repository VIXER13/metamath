#ifndef SYMDIFF_MINUS_HPP
#define SYMDIFF_MINUS_HPP

#include "negate.hpp"

namespace metamath::symdiff {

template<class E1, class E2>
class minus;

template<class U1, class U2>
using minus_type = std::conditional_t<
    is_integral_constant<U1>{} && is_integral_constant<U2>{},
    integral_constant<decltype(integral_constant_v<U1>{} - integral_constant_v<U2>{}), integral_constant_v<U1>{} - integral_constant_v<U2>{}>,
    std::conditional_t<
        is_integral_constant<U1>{} && !integral_constant_v<U1>{},
        negate<U2>,
        std::conditional_t<
            is_integral_constant<U2>{} && !integral_constant_v<U2>{},
            U1,
            std::conditional_t<
                is_constant<U1>{} && is_constant<U2>{},
                constant<decltype(constant_t<U1>{} - constant_t<U2>{})>,
                minus<U1, U2>
            >
        >
    >
>;

template<class E1, class E2>
class minus : public expression<minus<E1, E2>> {
    const E1 e1;
    const E2 e2;

public:
    template<uintmax_t M>
    using derivative_type = minus_type<typename E1::template derivative_type<M>, typename E2::template derivative_type<M>>;

    constexpr minus(const expression<E1>& e1, const expression<E2>& e2) :
        e1{e1()}, e2{e2()} {}

    template<class T>
    constexpr auto operator()(const T& x) const -> decltype(e1(x) - e2(x)) {
        return e1(x) - e2(x);
    }

    template<uintmax_t M>
    constexpr derivative_type<M> derivative() const {
        return e1.template derivative<M>() - e2.template derivative<M>();
    }
};

template<class T1, T1 N1, class T2, T2 N2>
constexpr integral_constant<decltype(N1 - N2), N1 - N2> operator-(const integral_constant<T1, N1>&, const integral_constant<T2, N2>&) {
    return integral_constant<decltype(N1 - N2), N1 - N2>{};
}

template<class T1, class E2>
constexpr negate<E2> operator-(const integral_constant<T1, 0>&, const expression<E2> &e) {
    return -e;
}

template<class E1, class T2>
constexpr const E1& operator-(const expression<E1> &e, const integral_constant<T2, 0>&) {
    return e();
}

template<class T1, class T2>
constexpr constant<decltype(T1{} - T2{})> operator-(const constant<T1>& e1, const constant<T2>& e2) {
    return constant<decltype(T1{} - T2{})>{e1.value - e2.value};
}

template<class T1, T1 N1, class T2>
constexpr std::enable_if_t<N1, constant<decltype(T1{} - T2{})>> operator-(const integral_constant<T1, N1>&, const constant<T2>& e2) {
    return constant<decltype(T1{} - T2{})>{N1 - e2.value};
}

template<class T1, class T2, T2 N2>
constexpr std::enable_if_t<N2, constant<decltype(T1{} - T2{})>> operator-(const constant<T1>& e1, const integral_constant<T2, N2>&) {
    return constant<decltype(T1{} - T2{})>{e1.value - N2};
}

template<class E1, class E2>
constexpr minus<E1, E2> operator-(const expression<E1> &e1, const expression<E2> &e2) {
    return minus<E1, E2>{e1, e2};
}

template<class T1, class E2>
constexpr std::enable_if_t<std::is_arithmetic_v<T1>, minus_type<constant<T1>, E2>> operator-(const T1& e1, const expression<E2>& e2) {
    return constant<T1>{e1} - e2();
}

template<class E1, class T2>
constexpr std::enable_if_t<std::is_arithmetic_v<T2>, minus_type<E1, constant<T2>>> operator-(const expression<E1>& e1, const T2& e2) {
    return e1() - constant<T2>{e2};
}

}

#endif