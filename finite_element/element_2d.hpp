#ifndef FINITE_ELEMENT_2D_HPP
#define FINITE_ELEMENT_2D_HPP

#include "element_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class Type, template<class> class Element_Type>
class element_2d : public virtual element_2d_base<Type>,
                   public Element_Type<Type> {
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::Nxi.size() &&
                  Element_Type<Type>::N.size() == Element_Type<Type>::Neta.size(),
                  "The number of functions and their derivatives does not match.");

public:
    size_t nodes_count() const override { return Element_Type<Type>::N.size(); }

    const std::array<Type, 2>& node(const size_t i) const override { return Element_Type<Type>::nodes[i]; }

    Type N   (const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::N   [i](xi); }
    Type Nxi (const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::Nxi [i](xi); }
    Type Neta(const size_t i, const std::array<Type, 2>& xi) const override { return Element_Type<Type>::Neta[i](xi); }

    Type boundary(const side_2d bound, const Type x) const override { return Element_Type<Type>::boundary(bound, x); }
};

}

#endif