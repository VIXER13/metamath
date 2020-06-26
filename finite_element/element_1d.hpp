#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_base.hpp"

namespace metamath::finite_element {

// Данная реализация подразумевает, что данные о функциях формы, их производных и геометрия элемента наследуются от класса Element_Type. 
// Таким образом пользователь сможет добавлять свои реализации конечных элементов не прибегая к дублированию интерфейса.
template<class Type, template<class> class Element_Type>
class element_1d : public virtual element_1d_base<Type>,
                   public Element_Type<Type> {
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::nodes.size(),
                  "The number of functions and nodes does not match.");
    static_assert(Element_Type<Type>::N.size() == Element_Type<Type>::Nxi.size(),
                  "The number of functions and their derivatives does not match.");
public:
    size_t nodes_count() const override { return Element_Type<Type>::N.size(); }

    Type node(const size_t i) const override { return Element_Type<Type>::nodes[i]; }
        
    Type N  (const size_t i, const std::array<Type, 1>& xi) const override { return Element_Type<Type>::N  [i](xi); }
    Type Nxi(const size_t i, const std::array<Type, 1>& xi) const override { return Element_Type<Type>::Nxi[i](xi); }

    Type boundary(const side_1d bound) const override { return Element_Type<Type>::boundary(bound); }
};

}

#endif