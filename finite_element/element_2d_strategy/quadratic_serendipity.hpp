#ifndef FINITE_ELEMENT_QUADRATIC_SERENDIPITY_ELEMENT_HPP
#define FINITE_ELEMENT_QUADRATIC_SERENDIPITY_ELEMENT_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class Type>
class quadratic_serendipity : public geometry_2d<Type, rectangle_element_geometry> {
    using geometry_2d<Type, rectangle_element_geometry>::xi;
    using geometry_2d<Type, rectangle_element_geometry>::eta;

public:
    Type get_parameter() const noexcept { return _p; }
    void set_parameter(const Type p) noexcept { _p = p; }

protected:
    explicit quadratic_serendipity() = default;

    // В серендиповой аппроксимации высших порядков возникает проблема с негативизмом стандартного базиса в угловых узлах.
    // Для этого вводится специальный параметр p, который позволяет её избежать.
    // В сущности p является значением интеграла по области элемента от угловой функции. Значение интегралов от промежуточных функций есть 1-p.
    Type _p = -1. / 3.; // Значение по умолчанию даёт нам классический вариант квадратичных серендиповых элементов.
    static constexpr symdiff::variable<2> p{}; // В силу особенностей вычисления выражений, дополнительному параметру необходимо завести дополнительную переменную.

    // Нумерация узлов на квадратичном серендиповом элементе: 6---5---4
    //                                                        |       |
    //                                                        7       3
    //                                                        |       |
    //                                                        0---1---2
    static constexpr std::array<std::array<Type, 2>, 8> nodes = { -1., -1.,
                                                                   0., -1.,
                                                                   1., -1.,
                                                                   1.,  0.,
                                                                   1.,  1.,
                                                                   0.,  1.,
                                                                  -1.,  1.,
                                                                  -1.,  0. };

    // Базисные функции в локальной системе координат имеют вид:
    // N_i = 0.0625 (1 + xi_i xi)(1 + eta_i eta)[(36p-1)(1 - xi_i xi - eta_i eta) + (36p+3)xi_i xi eta_i eta], xi_i = +-1, eta_i = +-1, i = 0,2,4,6,
    // N_i = 0.0625 (1 -  xi^2)(1 + eta_i eta)[(5-36p) + (36p+3)eta_i eta], eta_i = +-1, i = 1,5,
    // N_i = 0.0625 (1 - eta^2)(1 +  xi_i  xi)[(5-36p) + (36p+3) xi_i  xi],  xi_i = +-1, i = 3,7.
    static constexpr auto N1 =  (1.-xi)      * (1.-eta) * ((0.5625*p-0.0625)*(1.+xi+eta) + (0.5625*p+0.1875)*xi*eta);
    static constexpr auto N2 = -(1.-xi*xi)   * (1.-eta) * ((0.5625*p-0.3125)             + (0.5625*p+0.1875)*eta);
    static constexpr auto N3 =  (1.+xi)      * (1.-eta) * ((0.5625*p-0.0625)*(1.-xi+eta) - (0.5625*p+0.1875)*xi*eta);
    static constexpr auto N4 = -(1.-eta*eta) * (1.+xi)  * ((0.5625*p-0.3125)             - (0.5625*p+0.1875)*xi);
    static constexpr auto N5 =  (1.+xi)      * (1.+eta) * ((0.5625*p-0.0625)*(1.-xi-eta) + (0.5625*p+0.1875)*xi*eta);
    static constexpr auto N6 = -(1.-xi*xi)   * (1.+eta) * ((0.5625*p-0.3125)             - (0.5625*p+0.1875)*eta);
    static constexpr auto N7 =  (1.-xi)      * (1.+eta) * ((0.5625*p-0.0625)*(1.+xi-eta) - (0.5625*p+0.1875)*xi*eta);
    static constexpr auto N8 = -(1.-eta*eta) * (1.-xi)  * ((0.5625*p-0.3125)             + (0.5625*p+0.1875)*xi);

    static constexpr auto dN1_dxi = N1.template derivative<xi>();
    static constexpr auto dN2_dxi = N2.template derivative<xi>();
    static constexpr auto dN3_dxi = N3.template derivative<xi>();
    static constexpr auto dN4_dxi = N4.template derivative<xi>();
    static constexpr auto dN5_dxi = N1.template derivative<xi>();
    static constexpr auto dN6_dxi = N2.template derivative<xi>();
    static constexpr auto dN7_dxi = N3.template derivative<xi>();
    static constexpr auto dN8_dxi = N4.template derivative<xi>();

    static constexpr auto dN1_deta = N1.template derivative<eta>();
    static constexpr auto dN2_deta = N2.template derivative<eta>();
    static constexpr auto dN3_deta = N3.template derivative<eta>();
    static constexpr auto dN4_deta = N4.template derivative<eta>();
    static constexpr auto dN5_deta = N1.template derivative<eta>();
    static constexpr auto dN6_deta = N2.template derivative<eta>();
    static constexpr auto dN7_deta = N3.template derivative<eta>();
    static constexpr auto dN8_deta = N4.template derivative<eta>();

    static inline const std::array<std::function<Type(const std::array<Type, 3>&)>, 8>
        N    = { [](const std::array<Type, 3>& xi) { return N1(xi); },
                 [](const std::array<Type, 3>& xi) { return N2(xi); },
                 [](const std::array<Type, 3>& xi) { return N3(xi); },
                 [](const std::array<Type, 3>& xi) { return N4(xi); },
                 [](const std::array<Type, 3>& xi) { return N5(xi); },
                 [](const std::array<Type, 3>& xi) { return N6(xi); },
                 [](const std::array<Type, 3>& xi) { return N7(xi); },
                 [](const std::array<Type, 3>& xi) { return N8(xi); } },

        Nxi  = { [](const std::array<Type, 3>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN3_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN4_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN5_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN6_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN7_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN8_dxi(xi); } },

        Neta = { [](const std::array<Type, 3>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN3_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN4_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN5_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN6_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN7_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN8_deta(xi); } };
};

}

#endif