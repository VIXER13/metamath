#ifndef FINITE_ELEMENT_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class Type>
class triangle : public barycentric<Type> {
protected:
    using barycentric<Type>::xi;
    using barycentric<Type>::eta;
    using barycentric<Type>::L1;
    using barycentric<Type>::L2;
    using barycentric<Type>::L3;

    explicit triangle() = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          2--0
    */
    static constexpr std::array<std::array<Type, 2>, 3> nodes = { 1., 0.,
                                                                  0., 1.,
                                                                  0., 0. };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i, i = 1..3
    static constexpr auto N1 = L1;
    static constexpr auto N2 = L2;
    static constexpr auto N3 = L3;

    static constexpr auto dN1_dxi = N1 .template derivative<xi>();
    static constexpr auto dN2_dxi = N2 .template derivative<xi>();
    static constexpr auto dN3_dxi = N3 .template derivative<xi>();

    static constexpr auto dN1_deta = N1 .template derivative<eta>();
    static constexpr auto dN2_deta = N2 .template derivative<eta>();
    static constexpr auto dN3_deta = N3 .template derivative<eta>();

    static inline const std::array<std::function<Type(const std::array<Type, 2>&)>, 3>
        N    = { [](const std::array<Type, 2>& xi) { return N1(xi); },
                 [](const std::array<Type, 2>& xi) { return N2(xi); },
                 [](const std::array<Type, 2>& xi) { return N3(xi); } },

        Nxi  = { [](const std::array<Type, 2>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_dxi(xi); } },

        Neta = { [](const std::array<Type, 2>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_deta(xi); } };
};

}

#endif