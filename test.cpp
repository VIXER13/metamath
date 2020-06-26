#include <iostream>
#include <memory>
#include "finite_element/finite_element.hpp"

using namespace metamath::finite_element;

double check_integral(const std::unique_ptr<element_integrate_base<double>>& e) {
    double sum = 0.;
    for(size_t i = 0; i < e->nodes_count(); ++i)
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            sum += e->weight(q) * e->qN(i, q);
    return sum;
}

int main() {
    std::unique_ptr<element_integrate_base<double>>
        e = std::make_unique<element_1d_integrate<double, linear>>(quadrature<double, gauss1>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_1d_integrate<double, quadratic>>(quadrature<double, gauss2>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_1d_integrate<double, qubic>>(quadrature<double, gauss2>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, triangle>>(quadrature<double, gauss2>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, quadratic_triangle>>(quadrature<double, gauss2>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, qubic_triangle>>(quadrature<double, gauss3>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, bilinear>>(quadrature<double, gauss2>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, quadratic_serendipity>>(quadrature<double, gauss3>{});
    std::cout << check_integral(e) << std::endl;

    e = std::make_unique<element_2d_integrate<double, qubic_serendipity>>(quadrature<double, gauss4>{});
    std::cout << check_integral(e) << std::endl;

    return 0;
}