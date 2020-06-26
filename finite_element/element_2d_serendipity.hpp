#ifndef FINITE_ELEMENT_2D_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_SERENDIPITY_HPP

// В силу особенностей серендиповой аппроксимации, были выделены две специализации класса двумерных конечных элементов.

#include "element_base.hpp"
#include "element_2d_strategy/quadratic_serendipity.hpp"
#include "element_2d_strategy/qubic_serendipity.hpp"

namespace metamath::finite_element {

template<class Type, template<class> class Element_Type>
class element_2d_serendipity : public virtual element_2d_base<Type>,
                               public Element_Type<Type> {
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::Nxi.size() &&
                  Element_Type<Type>::N.size() == Element_Type<Type>::Neta.size(),
                  "The number of functions and their derivatives does not match.");

public:
    size_t nodes_count() const override { return Element_Type<Type>::N.size(); }

    const std::array<Type, 2>& node(const size_t i) const override { return Element_Type<Type>::nodes[i]; }

    Type N   (const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::N   [i]({xi[0], xi[1], Element_Type<Type>::_p}); }
    Type Nxi (const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::Nxi [i]({xi[0], xi[1], Element_Type<Type>::_p}); }
    Type Neta(const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::Neta[i]({xi[0], xi[1], Element_Type<Type>::_p}); }

    Type boundary(const side_2d bound, const Type x) const override { return Element_Type<Type>::boundary(bound, x); }
};

// Специализация под квадратичные серендиповы элементы
template<class Type>
class element_2d<Type, quadratic_serendipity> : public element_2d_serendipity<Type, quadratic_serendipity> {};

// Специализация под кубические серендиповы элементы
template<class Type>
class element_2d<Type, qubic_serendipity> : public element_2d_serendipity<Type, qubic_serendipity> {};

}

#endif