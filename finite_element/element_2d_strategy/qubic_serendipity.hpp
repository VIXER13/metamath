#ifndef FINITE_ELEMENT_QUBIC_SERENDIPITY_ELEMENT_HPP
#define FINITE_ELEMENT_QUBIC_SERENDIPITY_ELEMENT_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class Type>
class qubic_serendipity : public geometry_2d<Type, rectangle_element_geometry> {
    using geometry_2d<Type, rectangle_element_geometry>::xi;
    using geometry_2d<Type, rectangle_element_geometry>::eta;
    
public:
    Type get_parameter() const noexcept { return _p; }
    void set_parameter(const Type p) noexcept { _p = p; }    

protected:
    // В серендиповой аппроксимации высших порядков возникает проблема с негативизмом стандартного базиса в угловых узлах.
    // Для этого вводится специальный параметр p, который позволяет её избежать.
    // В сущности p является значением интеграла по области элемента от угловой функции. Значение интегралов от промежуточных функций есть (1-4p)/2.
    Type _p = -0.5; // Значение по умолчанию даёт нам классический вариант кубических серендиповых элементов.
    static constexpr symdiff::variable<2> p{}; // В силу особенностей вычисления выражений, дополнительному параметру необходимо завести дополнительную переменную.

    explicit qubic_serendipity() = default;

    // Нумерация узлов на кубическом серендиповом элементе: 9---8---7---6
    //                                                      |           |
    //                                                      10          5
    //                                                      |           |
    //                                                      11          4
    //                                                      |           |
    //                                                      0---1---2---3
    static constexpr std::array<std::array<Type, 2>, 12> nodes = {    -1.,    -1.,
                                                                   -1./3.,    -1.,
                                                                    1./3.,    -1.,
                                                                       1.,    -1.,
                                                                       1., -1./3.,
                                                                       1.,  1./3.,
                                                                       1.,     1.,
                                                                    1./3.,     1.,
                                                                   -1./3.,     1.,
                                                                      -1.,     1.,
                                                                      -1.,  1./3.,
                                                                      -1., -1./3. };

    static constexpr auto N1  =  (1.-xi)    * (1.-eta)     * (0.28125*(xi*xi+eta*eta + (2.00000*p+1.000000)*(xi*eta+xi+eta)) + 0.56250*p - 0.031250);
    static constexpr auto N2  = -(1.-xi*xi) * (1.-eta)     * (0.84375*xi             + (0.28125*p+0.140625)*eta              + 0.28125*p - 0.140625);
    static constexpr auto N3  =  (1.-xi*xi) * (1.-eta)     * (0.84375*xi             - (0.28125*p+0.140625)*eta              - 0.28125*p + 0.140625);
    static constexpr auto N4  =  (1.+xi)    * (1.-eta)     * (0.28125*(xi*xi+eta*eta - (2.00000*p+1.000000)*(xi*eta+xi-eta)) + 0.56250*p - 0.031250);
    static constexpr auto N5  = -(1.+xi)    * (1.-eta*eta) * (0.84375*eta            - (0.28125*p+0.140625)*xi               + 0.28125*p - 0.140625);
    static constexpr auto N6  =  (1.+xi)    * (1.-eta*eta) * (0.84375*eta            + (0.28125*p+0.140625)*xi               - 0.28125*p + 0.140625);
    static constexpr auto N7  =  (1.+xi)    * (1.+eta)     * (0.28125*(xi*xi+eta*eta + (2.00000*p+1.000000)*(xi*eta-xi-eta)) + 0.56250*p - 0.031250);
    static constexpr auto N8  =  (1.-xi*xi) * (1.+eta)     * (0.84375*xi             + (0.28125*p+0.140625)*eta              - 0.28125*p + 0.140625);
    static constexpr auto N9  = -(1.-xi*xi) * (1.+eta)     * (0.84375*xi             - (0.28125*p+0.140625)*eta              + 0.28125*p - 0.140625);
    static constexpr auto N10 =  (1.-xi)    * (1.+eta)     * (0.28125*(xi*xi+eta*eta - (2.00000*p+1.000000)*(xi*eta-xi+eta)) + 0.56250*p - 0.031250);
    static constexpr auto N11 =  (1.-xi)    * (1.-eta*eta) * (0.84375*eta            - (0.28125*p+0.140625)*xi               - 0.28125*p + 0.140625);
    static constexpr auto N12 = -(1.-xi)    * (1.-eta*eta) * (0.84375*eta            + (0.28125*p+0.140625)*xi               + 0.28125*p - 0.140625);

    static constexpr auto dN1_dxi  = N1 .template derivative<xi>();
    static constexpr auto dN2_dxi  = N2 .template derivative<xi>();
    static constexpr auto dN3_dxi  = N3 .template derivative<xi>();
    static constexpr auto dN4_dxi  = N4 .template derivative<xi>();
    static constexpr auto dN5_dxi  = N5 .template derivative<xi>();
    static constexpr auto dN6_dxi  = N6 .template derivative<xi>();
    static constexpr auto dN7_dxi  = N7 .template derivative<xi>();
    static constexpr auto dN8_dxi  = N8 .template derivative<xi>();
    static constexpr auto dN9_dxi  = N9 .template derivative<xi>();
    static constexpr auto dN10_dxi = N10.template derivative<xi>();
    static constexpr auto dN11_dxi = N11.template derivative<xi>();
    static constexpr auto dN12_dxi = N12.template derivative<xi>();

    static constexpr auto dN1_deta  = N1 .template derivative<eta>();
    static constexpr auto dN2_deta  = N2 .template derivative<eta>();
    static constexpr auto dN3_deta  = N3 .template derivative<eta>();
    static constexpr auto dN4_deta  = N4 .template derivative<eta>();
    static constexpr auto dN5_deta  = N5 .template derivative<eta>();
    static constexpr auto dN6_deta  = N6 .template derivative<eta>();
    static constexpr auto dN7_deta  = N7 .template derivative<eta>();
    static constexpr auto dN8_deta  = N8 .template derivative<eta>();
    static constexpr auto dN9_deta  = N9 .template derivative<eta>();
    static constexpr auto dN10_deta = N10.template derivative<eta>();
    static constexpr auto dN11_deta = N11.template derivative<eta>();
    static constexpr auto dN12_deta = N12.template derivative<eta>();

    // Базисные функции в локальной системе координат имеют вид:
    // N_i = 1/32 (1 + xi_i xi)(1 + eta_i eta)[9(xi^2 + eta^2) + (18p+9)(xi_i xi eta_i eta - xi_i xi - eta_i eta) + 18p-1], xi_i = +-1, eta_i = +-1, i = 0,3,6,9,
    // N_i = 9/64 (1 -  xi^2)(1 + eta_i eta)[18  xi_i  xi + (2p+1) eta_i eta - 1 + 2p], xi_i = +-1/3, eta_i = +-1  , i = 1,2,7,8,
    // N_i = 9/64 (1 - eta^2)(1 +  xi_i  xi)[18 eta_i eta + (2p+1)  xi_i  xi - 1 + 2p], xi_i = +-1  , eta_i = +-1/3, i = 4,5,10,11.
    static inline const std::array<std::function<Type(const std::array<Type, 3>&)>, 12>
        N    = { [](const std::array<Type, 3>& xi) { return N1(xi); },
                 [](const std::array<Type, 3>& xi) { return N2(xi); },
                 [](const std::array<Type, 3>& xi) { return N3(xi); },
                 [](const std::array<Type, 3>& xi) { return N4(xi); },
                 [](const std::array<Type, 3>& xi) { return N5(xi); },
                 [](const std::array<Type, 3>& xi) { return N6(xi); },
                 [](const std::array<Type, 3>& xi) { return N7(xi); },
                 [](const std::array<Type, 3>& xi) { return N8(xi); },
                 [](const std::array<Type, 3>& xi) { return N9(xi); },
                 [](const std::array<Type, 3>& xi) { return N10(xi); },
                 [](const std::array<Type, 3>& xi) { return N11(xi); },
                 [](const std::array<Type, 3>& xi) { return N12(xi); } },

        Nxi  = { [](const std::array<Type, 3>& xi) { return dN1_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN2_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN3_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN4_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN5_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN6_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN7_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN8_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN9_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN10_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN11_dxi(xi); },
                 [](const std::array<Type, 3>& xi) { return dN12_dxi(xi); } },

        Neta = { [](const std::array<Type, 3>& xi) { return dN1_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN2_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN3_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN4_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN5_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN6_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN7_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN8_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN9_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN10_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN11_deta(xi); },
                 [](const std::array<Type, 3>& xi) { return dN12_deta(xi); } };
};

}

#endif