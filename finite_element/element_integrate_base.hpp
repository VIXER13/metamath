#ifndef FINITE_ELEMENT_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_INTEGRATE_BASE_HPP

#include <vector>
#include "element_base.hpp"
#include "quadrature.hpp"

namespace metamath::finite_element {

template<class Type>
class element_integrate_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

protected:
    std::vector<Type> _weights, _qN;
    explicit element_integrate_base() noexcept = default;

public:
    size_t qnodes_count() const noexcept { return _weights.size(); }
    size_t  nodes_count() const noexcept { return _qN.size() / qnodes_count(); }

    Type weight(const size_t q) const noexcept { return _weights[q]; }

    Type qN(const size_t i, const size_t q) const noexcept { return _qN[i*qnodes_count() + q]; }

    virtual ~element_integrate_base() noexcept = default;
};

template<class Type>
class element_1d_integrate_base : public element_integrate_base<Type>,
                                  private virtual element_1d_base<Type> {
protected:
    std::vector<Type> _qNxi;

public:
    using element_integrate_base<Type>::nodes_count;
    using element_integrate_base<Type>::qnodes_count;
    using element_1d_base<Type>::boundary;
    using element_1d_base<Type>::node;
    using element_1d_base<Type>::N;
    using element_1d_base<Type>::Nxi;
    
    virtual void set_quadrature(const quadrature_base<Type>& quadrature) = 0;

    Type qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }
};

template<class Type>
class element_2d_integrate_base : public element_integrate_base<Type>,
                                  public virtual element_2d_base<Type> {
protected:
    std::vector<Type> _qNxi, _qNeta;

public:
    using element_integrate_base<Type>::nodes_count;
    using element_integrate_base<Type>::qnodes_count;
    using element_2d_base<Type>::boundary;
    using element_2d_base<Type>::node;
    using element_2d_base<Type>::N;
    using element_2d_base<Type>::Nxi;
    using element_2d_base<Type>::Neta;
    
    virtual void set_quadrature(const quadrature_base<Type>& quadrature_xi, const quadrature_base<Type>& quadrature_eta) = 0;

    Type qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    Type qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }
};

}

#endif