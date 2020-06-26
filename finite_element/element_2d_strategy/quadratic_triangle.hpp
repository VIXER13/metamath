#ifndef FINITE_ELEMENT_QUADRATIC_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_QUADRATIC_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class Type>
class quadratic_triangle : public barycentric<Type> {
protected:
    using barycentric<Type>::xi;
    using barycentric<Type>::eta;
    using barycentric<Type>::L1;
    using barycentric<Type>::L2;
    using barycentric<Type>::L3;

    explicit quadratic_triangle() = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          4  3
                                                          |    \
                                                          2--5--0
    */
    static constexpr std::array<std::array<Type, 2>, 6> nodes = { 1.0, 0.0,
                                                                  0.0, 1.0,
                                                                  0.0, 0.0,
                                                                  0.5, 0.5,
                                                                  0.0, 0.5,
                                                                  0.5, 0.0 };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i (2 L_i - 1), i = 1..3,
    //                                                            N_3 = 4 L_1 L_2,
    //                                                            N_4 = 4 L_2 L_3,
    //                                                            N_5 = 4 L_3 L_1.
    static constexpr auto N1 = L1 * (2. * L1 - 1);
    static constexpr auto N2 = L2 * (2. * L2 - 1);
    static constexpr auto N3 = L3 * (2. * L3 - 1);
    static constexpr auto N4 = 4. * L1 * L2;
    static constexpr auto N5 = 4. * L2 * L3;
    static constexpr auto N6 = 4. * L3 * L1;

    static constexpr auto dN1_dxi = N1 .template derivative<xi>();
    static constexpr auto dN2_dxi = N2 .template derivative<xi>();
    static constexpr auto dN3_dxi = N3 .template derivative<xi>();
    static constexpr auto dN4_dxi = N4 .template derivative<xi>();
    static constexpr auto dN5_dxi = N5 .template derivative<xi>();
    static constexpr auto dN6_dxi = N6 .template derivative<xi>();

    static constexpr auto dN1_deta = N1 .template derivative<eta>();
    static constexpr auto dN2_deta = N2 .template derivative<eta>();
    static constexpr auto dN3_deta = N3 .template derivative<eta>();
    static constexpr auto dN4_deta = N4 .template derivative<eta>();
    static constexpr auto dN5_deta = N5 .template derivative<eta>();
    static constexpr auto dN6_deta = N6 .template derivative<eta>();

    static inline const std::array<std::function<Type(const std::array<Type, 2>&)>, 6>
        N    = { [](const std::array<Type, 2>& xi) { return N1(xi); },
                 [](const std::array<Type, 2>& xi) { return N2(xi); },
                 [](const std::array<Type, 2>& xi) { return N3(xi); },
                 [](const std::array<Type, 2>& xi) { return N4(xi); },
                 [](const std::array<Type, 2>& xi) { return N5(xi); },
                 [](const std::array<Type, 2>& xi) { return N6(xi); } },

        Nxi  = { [](const std::array<Type, 2>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN5_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN6_dxi(xi); } },

        Neta = { [](const std::array<Type, 2>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN5_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN6_deta(xi); } };
};

}

#endif