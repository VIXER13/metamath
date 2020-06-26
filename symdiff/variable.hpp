#ifndef SYMDIFF_VARIABLE_HPP
#define SYMDIFF_VARIABLE_HPP

#include "integral_constant.hpp"

namespace metamath::symdiff {

template<uintmax_t N>
struct variable : expression<variable<N>> {
    template<uintmax_t M>
    using derivative_type = integral_constant<int8_t, N == M>;

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(x[N]) {
        return x[N];
    }

    template<uintmax_t M>
    constexpr derivative_type<M> derivative() const {
        return derivative_type<M>{};
    }

    constexpr operator uintmax_t() const {
        return N;
    }
};

}

#endif