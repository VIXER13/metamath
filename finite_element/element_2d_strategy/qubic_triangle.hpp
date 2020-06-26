#ifndef FINITE_ELEMENT_QUBIC_TRIANGLE_ELEMENT_HPP
#define FINITE_ELEMENT_QUBIC_TRIANGLE_ELEMENT_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class Type>
class qubic_triangle : public barycentric<Type> {
protected:
    using barycentric<Type>::xi;
    using barycentric<Type>::eta;
    using barycentric<Type>::L1;
    using barycentric<Type>::L2;
    using barycentric<Type>::L3;

    explicit qubic_triangle() = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          5 4
                                                          |  \
                                                          6 9 3
                                                          |    \
                                                          2-7-8-0
    */
    static constexpr std::array<std::array<Type, 2>, 10> nodes = {    1.,    0.,
                                                                      0.,    1.,
                                                                      0.,    0.,
                                                                   2./3., 1./3.,
                                                                   1./3., 2./3.,
                                                                      0., 2./3.,
                                                                      0., 1./3.,
                                                                   1./3.,    0.,
                                                                   2./3.,    0.,
                                                                   1./3., 1./3. };

    // Базисные функции в барицентрических координатах имеют вид: N_i = 0.5 L_i (3 L_i - 1)(3 L_i - 2), i = 1..3,
    //                                                            N_3 = 4.5 L_1 L_2 (3 L_1 - 1),
    //                                                            N_4 = 4.5 L_1 L_2 (3 L_2 - 1),
    //                                                            N_5 = 4.5 L_2 L_3 (3 L_2 - 1),
    //                                                            N_6 = 4.5 L_2 L_3 (3 L_3 - 1),
    //                                                            N_7 = 4.5 L_3 L_1 (3 L_3 - 1),
    //                                                            N_8 = 4.5 L_3 L_1 (3 L_1 - 1),
    //                                                            N_9 =  27 L_1 L_2 L_3
    static constexpr auto N1  = 0.5 * L1 * (3.*L1 - 1.) * (3.*L1 - 2.);
    static constexpr auto N2  = 0.5 * L2 * (3.*L2 - 1.) * (3.*L2 - 2.);
    static constexpr auto N3  = 0.5 * L3 * (3.*L3 - 1.) * (3.*L3 - 2.);
    static constexpr auto N4  = 4.5 * L1 * L2 * (3.*L1 - 1.);
    static constexpr auto N5  = 4.5 * L1 * L2 * (3.*L2 - 1.);
    static constexpr auto N6  = 4.5 * L2 * L3 * (3.*L2 - 1.);
    static constexpr auto N7  = 4.5 * L2 * L3 * (3.*L3 - 1.);
    static constexpr auto N8  = 4.5 * L3 * L1 * (3.*L3 - 1.);
    static constexpr auto N9  = 4.5 * L3 * L1 * (3.*L1 - 1.);
    static constexpr auto N10 = 27. * L1 * L2 * L3;

    static constexpr auto dN1_dxi  = N1 .template derivative<xi>();
    static constexpr auto dN2_dxi  = N2 .template derivative<xi>();
    static constexpr auto dN3_dxi  = N3 .template derivative<xi>();
    static constexpr auto dN4_dxi  = N4 .template derivative<xi>();
    static constexpr auto dN5_dxi  = N5 .template derivative<xi>();
    static constexpr auto dN6_dxi  = N6 .template derivative<xi>();
    static constexpr auto dN7_dxi  = N7 .template derivative<xi>();
    static constexpr auto dN8_dxi  = N8 .template derivative<xi>();
    static constexpr auto dN9_dxi  = N9 .template derivative<xi>();
    static constexpr auto dN10_dxi = N10 .template derivative<xi>();

    static constexpr auto dN1_deta  = N1 .template derivative<eta>();
    static constexpr auto dN2_deta  = N2 .template derivative<eta>();
    static constexpr auto dN3_deta  = N3 .template derivative<eta>();
    static constexpr auto dN4_deta  = N4 .template derivative<eta>();
    static constexpr auto dN5_deta  = N5 .template derivative<eta>();
    static constexpr auto dN6_deta  = N6 .template derivative<eta>();
    static constexpr auto dN7_deta  = N7 .template derivative<eta>();
    static constexpr auto dN8_deta  = N8 .template derivative<eta>();
    static constexpr auto dN9_deta  = N9 .template derivative<eta>();
    static constexpr auto dN10_deta = N10 .template derivative<eta>();

    static inline const std::array<std::function<Type(const std::array<Type, 2>&)>, 10>
        N    = { [](const std::array<Type, 2>& xi) { return N1(xi); },
                 [](const std::array<Type, 2>& xi) { return N2(xi); },
                 [](const std::array<Type, 2>& xi) { return N3(xi); },
                 [](const std::array<Type, 2>& xi) { return N4(xi); },
                 [](const std::array<Type, 2>& xi) { return N5(xi); },
                 [](const std::array<Type, 2>& xi) { return N6(xi); },
                 [](const std::array<Type, 2>& xi) { return N7(xi); },
                 [](const std::array<Type, 2>& xi) { return N8(xi); },
                 [](const std::array<Type, 2>& xi) { return N9(xi); },
                 [](const std::array<Type, 2>& xi) { return N10(xi); } },

        Nxi  = { [](const std::array<Type, 2>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN5_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN6_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN7_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN8_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN9_dxi(xi); },
                 [](const std::array<Type, 2>& xi) { return dN10_dxi(xi); } },

        Neta = { [](const std::array<Type, 2>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN3_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN4_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN5_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN6_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN7_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN8_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN9_deta(xi); },
                 [](const std::array<Type, 2>& xi) { return dN10_deta(xi); } };
};

}

#endif