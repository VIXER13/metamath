#ifndef FINITE_ELEMENT_2D_BASIS_QUADRATIC_TRIANGLE_HPP
#define FINITE_ELEMENT_2D_BASIS_QUADRATIC_TRIANGLE_HPP

#include "barycentric.hpp"

namespace metamath::finite_element {

template<class T>
class quadratic_triangle : public barycentric<T> {
protected:
    using barycentric<T>::xi;
    using barycentric<T>::eta;
    using barycentric<T>::L1;
    using barycentric<T>::L2;
    using barycentric<T>::L3;

    explicit quadratic_triangle() = default;
    ~quadratic_triangle() override = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          4  3
                                                          |    \
                                                          2--5--0
    */
    static constexpr std::array<std::array<T, 2>, 6> nodes = { T{1.0}, T{0.0},
                                                               T{0.0}, T{1.0},
                                                               T{0.0}, T{0.0},
                                                               T{0.5}, T{0.5},
                                                               T{0.0}, T{0.5},
                                                               T{0.5}, T{0.0} };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i (2 L_i - 1), i = 1..3,
    //                                                            N_3 = 4 L_1 L_2,
    //                                                            N_4 = 4 L_2 L_3,
    //                                                            N_5 = 4 L_3 L_1.
    static constexpr auto basis = std::make_tuple(
        L1   * (T{2} * L1 - 1),
        L2   * (T{2} * L2 - 1),
        L3   * (T{2} * L3 - 1),
        T{4} *    L1 * L2,
        T{4} *    L2 * L3,
        T{4} *    L3 * L1
    );
};

}

#endif