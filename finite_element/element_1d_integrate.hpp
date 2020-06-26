#ifndef FINITE_ELEMENT_1D_INTEGRATE_HPP
#define FINITE_ELEMENT_1D_INTEGRATE_HPP

#include "element_integrate_base.hpp"
#include "element_1d.hpp"

namespace metamath::finite_element {

template<class Type, template<class> class Element_Type>
class element_1d_integrate : public element_1d_integrate_base<Type>,
                             public element_1d<Type, Element_Type> {
    using element_1d_integrate_base<Type>::_weights;
    using element_1d_integrate_base<Type>::_qN;
    using element_1d_integrate_base<Type>::_qNxi;

public:
    using element_1d_integrate_base<Type>::qnodes_count;
    using element_1d_integrate_base<Type>::nodes_count;
    using element_1d<Type, Element_Type>::boundary;
    using element_1d<Type, Element_Type>::N;
    using element_1d<Type, Element_Type>::Nxi;

    explicit element_1d_integrate(const quadrature_base<Type>& quadrature) { set_quadrature(quadrature); }

    void set_quadrature(const quadrature_base<Type>& quadrature) override {
        std::vector<Type> xi(quadrature.nodes_count());
        _weights.resize(quadrature.nodes_count());
        Type jacobian = (           boundary(side_1d::RIGHT) -            boundary(side_1d::LEFT)) /
                        (quadrature.boundary(side_1d::RIGHT) - quadrature.boundary(side_1d::LEFT));
        for(size_t q = 0; q < quadrature.nodes_count(); ++q) {
            xi[q] = boundary(side_1d::LEFT) + (quadrature.node(q) - quadrature.boundary(side_1d::LEFT)) * jacobian;
            _weights[q] = quadrature.weight(q) * jacobian;
        }

        _qN  .resize(element_1d<Type, Element_Type>::nodes_count() * qnodes_count());
        _qNxi.resize(element_1d<Type, Element_Type>::nodes_count() * qnodes_count());
        for(size_t i = 0; i < nodes_count(); ++i) {
            for(size_t q = 0; q < qnodes_count(); ++q) {
                _qN  [i*qnodes_count() + q] = N  (i, {xi[q]});
                _qNxi[i*qnodes_count() + q] = Nxi(i, {xi[q]});
            }
        }
    }
};

}

#endif