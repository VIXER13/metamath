#include "metamath.hpp"
#include <iostream>
#include <array>
#include <vector>

namespace {

// Пример использования библиотеки symdiff
// Важно понимать, что выражения и их производные вычисляются на этапе компиляции
template<class Linear_Container>
void symdiff_example(const Linear_Container& point) {
    // Создаём 3 переменные (важно, они должны быть разных типов, т.е. иметь разные номера)
    static constexpr metamath::symdiff::variable<0> x{};
    static constexpr metamath::symdiff::variable<1> y{};
    static constexpr metamath::symdiff::variable<2> z{};

    // Составляем выражение, оно может использовать многие другие функции по типу sin, cos, tan и т.д.
    static constexpr auto f = x * y * exp(2 * x * z);

    // Берём смешанную производную по переменным x и z
    // Аналитическая производная: 4 * x * y * (1 + x * z) * exp(2 * x * z)
    static constexpr auto df = metamath::symdiff::derivative<x, z>(f);

    // Вычисляем
    std::cout << "x = (" << point[0] << ", " << point[1] << ", " << point[2] << ')' << std::endl;
    std::cout << "f(x)  = " <<  f(point) << std::endl;
    std::cout << "df(x) = " << df(point) << std::endl;
}

// Тест для квадратур
template<class T>
void quadrature_test(const metamath::finite_element::quadrature_1d_base<T>& quadrature) {
    T sum = T{0};
    std::cout << "nodes count: " << quadrature.nodes_count() << std::endl;
    std::cout << "nodes coords: ";
    for(size_t i = 0; i < quadrature.nodes_count(); ++i) {
        std::cout << quadrature.node(i)[0] << ' ';
        sum += quadrature.weight(i);
    }
    std::cout << std::endl;
    std::cout << "weights sum = " << sum << std::endl; // Для гауссовых квадратур сумма должна быть равна 2
}

template<class T>
void element_integrate_test(const metamath::finite_element::element_integrate_base<T>& element) {
    std::cout << "nodes count: " << element.nodes_count() << std::endl;
    std::cout << "quadrature nodes count: " << element.qnodes_count() << std::endl;
}

// Тест для одномерных элементов
template<class T>
void element_1d_integrate_test(const metamath::finite_element::element_1d_integrate_base<T>& element) {
    T len = 0, sum_dxi = 0;
    element_integrate_test(element);
    std::cout << "nodes coords: ";
    for(size_t i = 0; i < element.nodes_count(); ++i) {
        std::cout << element.node(i) << ' ';
        for(size_t q = 0; q < element.qnodes_count(); ++q) {
            len     += element.weight(q) * element.qN  (i, q);
            sum_dxi += element.weight(q) * element.qNxi(i, q);
        }
    }
    std::cout << std::endl;
    std::cout << "len     = " << len    << std::endl   // Длина одномерных элементов равна 2
              << "sum_dxi = " << sum_dxi << std::endl; // Сумма всех производных во всех квадратурных узлах равна 0
}

// Тест для двумерных элементов
template<class T>
void element_2d_integrate_test(const metamath::finite_element::element_2d_integrate_base<T>& element) {
    T area = 0, sum_dxi = 0, sum_deta = 0;
    element_integrate_test(element);
    std::cout << "nodes coords: ";
    for(size_t i = 0; i < element.nodes_count(); ++i) {
        const std::array<T, 2>& node = element.node(i);
        std::cout << "(" << node[0] << ", " << node[1] << ") ";
        for(size_t q = 0; q < element.qnodes_count(); ++q) {
            area     += element.weight(q) * element.qN   (i, q);
            sum_dxi  += element.weight(q) * element.qNxi (i, q);
            sum_deta += element.weight(q) * element.qNeta(i, q);
        }
    }
    std::cout << std::endl;
    std::cout << "area     = " << area    << std::endl   // Площадь элемента, треугольные 0.5, четырёхугольные 4
              << "sum_dxi  = " << sum_dxi  << std::endl  // Сумма всех производных во всех квадратурных узлах равна 0
              << "sum_deta = " << sum_deta << std::endl; // Сумма всех производных во всех квадратурных узлах равна 0
}

}

int main() {
    std::cout << "SYMDIFF_TEST:" << std::endl;
    symdiff_example(std::array{1., 2., 3.});
    std::cout << std::endl;
    symdiff_example(std::vector{0.1, 0.2, 0.3});

    std::cout << std::endl << std::endl;

    using namespace metamath::finite_element;

    // quadrature_1d_base<T> --- абстрактный класс
    // quadrature_1d<T, Quadrature_Type> --- конкретный класс квадратуры
    // T --- тип данных; Quadrature_Type --- класс стратегии, который включает в себя данные о квадратуре
    std::cout << "QUADRATURE_TEST:" << std::endl;
    std::cout << "gauss1:" << std::endl;
    quadrature_test(quadrature_1d<double, gauss1>{});
    std::cout << std::endl;
    std::cout << "gauss2:" << std::endl;
    quadrature_test(quadrature_1d<double, gauss2>{});
    std::cout << std::endl;
    std::cout << "gauss3:" << std::endl;
    quadrature_test(quadrature_1d<double, gauss3>{});
    std::cout << std::endl;
    std::cout << "gauss4:" << std::endl;
    quadrature_test(quadrature_1d<double, gauss4>{});
    std::cout << std::endl;
    std::cout << "gauss5:" << std::endl;
    quadrature_test(quadrature_1d<double, gauss5>{});

    std::cout << std::endl << std::endl;

    // element_integrate_base<T> --- абстрактный класс
    // element_1d_integrate_base<T> --- абстрактный класс, наследуемый от element_integrate_base<T>
    // element_1d_integrate<T, Element_Type> --- конкретный класс, наследуемый от element_1d_integrate_base<T>
    // T --- тип данных; Element_Type --- класс стратегия, который включает в себя данные об элементе
    // Отметим, что производные функций форм вычисляются аналитически при помощи модуля symdiff
    std::cout << "ELEMENT_1D_INTEGRATE_TEST: " << std::endl;
    std::cout << "linear: " << std::endl;
    element_1d_integrate_test(element_1d_integrate<double, linear>{quadrature_1d<double, gauss1>{}});
    std::cout << std::endl;
    std::cout << "Quadratic: " << std::endl;
    element_1d_integrate_test(element_1d_integrate<double, quadratic>{quadrature_1d<double, gauss2>{}});
    std::cout << std::endl;
    std::cout << "Qubic: " << std::endl;
    element_1d_integrate_test(element_1d_integrate<double, qubic>{quadrature_1d<double, gauss2>{}});


    std::cout << std::endl << std::endl;

    // КОМПИЛЯЦИЯ НЕКОТОРЫХ ЭЛЕМЕНТОВ ВЫСШЕГО ПОРЯДКА МОЖЕТ ПРИВЕСТИ К СЛИШКОМ ДОЛГОЙ КОМПИЛЯЦИИ

    // element_integrate_base<T> --- абстрактный класс
    // element_2d_integrate_base<T> --- абстрактный класс, наследуемый от element_integrate_base<T>
    // element_2d_integrate<T, Element_Type> --- конкретный класс, наследуемый от element_2d_integrate_base<T>
    // T --- тип данных; Element_Type --- класс стратегия, который включает в себя данные об элементе
    // Отметим, что производные функций форм вычисляются аналитически при помощи модуля symdiff
    std::cout << "ELEMENT_2D_INTEGRATE_TEST: " << std::endl;
    std::cout << "triangle:" << std::endl;
    element_2d_integrate_test(element_2d_integrate<double, triangle>{quadrature_1d<double, gauss2>{}});
    std::cout << std::endl;
    std::cout << "quadratic_triangle:" << std::endl;
    element_2d_integrate_test(element_2d_integrate<double, quadratic_triangle>{quadrature_1d<double, gauss3>{}});
//    std::cout << std::endl;
//    std::cout << "qubic_triangle:" << std::endl;
//    element_2d_integrate_test(element_2d_integrate<double, qubic_triangle>{quadrature_1d<double, gauss4>{}});
    std::cout << std::endl;
    std::cout << "bilinear:" << std::endl;
    element_2d_integrate_test(element_2d_integrate<double, bilinear>{quadrature_1d<double, gauss2>{}});
    std::cout << std::endl;
    std::cout << "quadratic_serendipity:" << std::endl;
    element_2d_integrate_test(element_2d_integrate<double, quadratic_serendipity>{quadrature_1d<double, gauss3>{}});
    std::cout << std::endl;
    std::cout << "quadratic_lagrange:" << std::endl;
    element_2d_integrate_test(element_2d_integrate<double, quadratic_lagrange>{quadrature_1d<double, gauss3>{}});
//    std::cout << std::endl;
//    std::cout << "qubic_serendipity:" << std::endl;
//    element_2d_integrate_test(element_2d_integrate<double, qubic_serendipity>{quadrature_1d<double, gauss4>{}});
//    std::cout << std::endl;
//    std::cout << "quartic_serendipity:" <<std::endl;
//    element_2d_integrate_test(element_2d_integrate<double,   quartic_serendipity>{quadrature_1d<double, gauss4>{}});
//    std::cout << std::endl;
//    std::cout << "quintic_serendipity:" <<std::endl;
//    element_2d_integrate_test(element_2d_integrate<double,   quintic_serendipity>{quadrature_1d<double, gauss4>{}});

    return EXIT_SUCCESS;
}