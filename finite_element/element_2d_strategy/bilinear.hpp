#ifndef FINITE_ELEMENT_BILINEAR_ELEMENT_HPP
#define FINITE_ELEMENT_BILINEAR_ELEMENT_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class Type>
class bilinear : public geometry_2d<Type, rectangle_element_geometry> {
    using geometry_2d<Type, rectangle_element_geometry>::xi;
    using geometry_2d<Type, rectangle_element_geometry>::eta;

protected:
    explicit bilinear() = default;

    // Нумерация узлов на билинейном элементе: 3---2
    //                                         |   |
    //                                         0---1
    static constexpr std::array<std::array<Type, 2>, 4> nodes = { -1., -1.,
                                                                   1., -1.,
                                                                   1.,  1.,
                                                                  -1.,  1. };

    // Базисные функции в локальной системе координат имеют вид: N_i = 0.25 (1 + xi_i x)(1 + eta_i eta), xi_i =+-1, eta_i = +-1, i = 0..3
    static constexpr auto N1 = 0.25 * (1. - xi) * (1. - eta);
    static constexpr auto N2 = 0.25 * (1. + xi) * (1. - eta);
    static constexpr auto N3 = 0.25 * (1. + xi) * (1. + eta);
    static constexpr auto N4 = 0.25 * (1. - xi) * (1. + eta);

    static constexpr auto dN1_dxi = N1.template derivative<xi>();
    static constexpr auto dN2_dxi = N2.template derivative<xi>();
    static constexpr auto dN3_dxi = N3.template derivative<xi>();
    static constexpr auto dN4_dxi = N4.template derivative<xi>();

    static constexpr auto dN1_deta = N1.template derivative<eta>();
    static constexpr auto dN2_deta = N2.template derivative<eta>();
    static constexpr auto dN3_deta = N3.template derivative<eta>();
    static constexpr auto dN4_deta = N4.template derivative<eta>();

    static inline const std::array<std::function<Type(const std::array<Type, 2>&)>, 4>
        N    = { [](const std::array<Type, 2>& xi) { return N1(xi); },
                 [](const std::array<Type, 2>& xi) { return N2(xi); },
                 [](const std::array<Type, 2>& xi) { return N3(xi); },
                 [](const std::array<Type, 2>& xi) { return N4(xi); } },

        Nxi  = { [](const std::array<Type, 2>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_dxi(xi); } },

        Neta = { [](const std::array<Type, 2>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_deta(xi); } };
};

}

#endif