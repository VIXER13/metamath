#ifndef FINITE_ELEMENT_QUBIC_ELEMENT_HPP
#define FINITE_ELEMENT_QUBIC_ELEMENT_HPP

#include <functional>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class Type>
class qubic : protected geometry_1d<Type, standart_segment_geometry> {
protected:
    using geometry_1d<Type, standart_segment_geometry>::xi;

    explicit qubic() noexcept = default;

    // Нумерация узлов на кубическом элементе: 0--1--2--3
    static constexpr std::array<Type, 4> nodes = { -1., -1./3., 1./3., 1. };

    static constexpr auto N1 = -0.5625 * (xi -    1.) * (xi*xi - 1./9.);
    static constexpr auto N2 =  1.6875 * (xi - 1./3.) * (xi*xi -    1.);
    static constexpr auto N3 = -1.6875 * (xi + 1./3.) * (xi*xi -    1.);
    static constexpr auto N4 =  0.5625 * (xi +    1.) * (xi*xi - 1./9.);

    static constexpr auto dN1_dxi = N1.template derivative<xi>();
    static constexpr auto dN2_dxi = N2.template derivative<xi>();
    static constexpr auto dN3_dxi = N3.template derivative<xi>();
    static constexpr auto dN4_dxi = N4.template derivative<xi>();

    static inline const std::array<std::function<Type(const std::array<Type, 1>&)>, 4>
        N   = { [](const std::array<Type, 1>& xi) { return N1(xi); },
                [](const std::array<Type, 1>& xi) { return N2(xi); },
                [](const std::array<Type, 1>& xi) { return N3(xi); },
                [](const std::array<Type, 1>& xi) { return N4(xi); } },

        Nxi = { [](const std::array<Type, 1>& xi) { return dN1_dxi(xi); },
                [](const std::array<Type, 1>& xi) { return dN2_dxi(xi); },
                [](const std::array<Type, 1>& xi) { return dN3_dxi(xi); },
                [](const std::array<Type, 1>& xi) { return dN4_dxi(xi); } };
};

}

#endif