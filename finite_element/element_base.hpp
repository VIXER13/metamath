#ifndef FINITE_ELEMENT_BASE_HPP
#define FINITE_ELEMENT_BASE_HPP

#include "element_1d_strategy/geometry_1d.hpp"
#include "element_2d_strategy/geometry_2d.hpp"

namespace metamath::finite_element {

class element_base {
public:
    virtual size_t nodes_count() const = 0; // Любой конечный элемент, вне зависимости от его размерности, имеет некоторое количество узлов.
    virtual ~element_base() noexcept = default;
};

template<class Type>
class element_1d_base : public element_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual Type node(const size_t i) const = 0;

    virtual Type N  (const size_t i, const std::array<Type, 1>& xi) const = 0; // Обращение к i-ой функции формы в точке xi.
    virtual Type Nxi(const size_t i, const std::array<Type, 1>& xi) const = 0; // Аналогично для производной.

    virtual Type boundary(const side_1d bound) const = 0; // Геометрия элемента.
};

template<class Type>
class element_2d_base : public element_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

public:
    virtual const std::array<Type, 2>& node(const size_t i) const = 0;

    virtual Type N   (const size_t i, const std::array<Type, 2>& xi) const = 0; // Обращение к i-ой функции формы в точке (xi, eta).
    virtual Type Nxi (const size_t i, const std::array<Type, 2>& xi) const = 0; // Аналогично для производных.
    virtual Type Neta(const size_t i, const std::array<Type, 2>& xi) const = 0;

    virtual Type boundary(const side_2d bound, const Type x) const = 0; // Геометрия элемента.
};

}

#endif