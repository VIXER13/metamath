#ifndef FINITE_ELEMENT_LINEAR_ELEMENT_HPP
#define FINITE_ELEMENT_LINEAR_ELEMENT_HPP

#include <functional>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class Type>
class linear : public geometry_1d<Type, standart_segment_geometry> {
protected:
    using geometry_1d<Type, standart_segment_geometry>::xi;

    explicit linear() noexcept = default;

    // Нумерация узлов на линейном элементе: 0--1
    static constexpr std::array<Type, 2> nodes = { -1., 1. };

    static constexpr auto N1 = 0.5 * (1. - xi);
    static constexpr auto N2 = 0.5 * (1. + xi);

    static constexpr auto dN1_dxi = N1.template derivative<xi>();
    static constexpr auto dN2_dxi = N2.template derivative<xi>();

    static inline const std::array<std::function<Type(const std::array<Type, 1>&)>, 2>
        N   = { [](const std::array<Type, 1>& xi) { return N1(xi); },
                [](const std::array<Type, 1>& xi) { return N2(xi); } },

        Nxi = { [](const std::array<Type, 1>& xi) { return dN1_dxi(xi); },
                [](const std::array<Type, 1>& xi) { return dN2_dxi(xi); } };
};

}

#endif