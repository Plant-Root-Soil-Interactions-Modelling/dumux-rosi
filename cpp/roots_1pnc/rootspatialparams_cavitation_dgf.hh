// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOT_SPATIALPARAMS1PNC_DGF_HH
#define DUMUX_ROOT_SPATIALPARAMS1PNC_DGF_HH

#include <cmath>

#include <dune/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/io/inputfilefunction.hh> // in dumux-rosi

#include "../roots_1p/rootspatialparams_dgf.hh"

namespace Dumux {

/*!
 * \brief Root spatial parameters class for grid files.
 *
 * use initParameters to initialize the class with data from the grid file
 *
 * Permeability takes the effect of xylem caviation into account according to the model of Spery et al.
 *
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParamsCaviationDGF : public RootSpatialParamsDGF<FVGridGeometry, Scalar> {

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using Water = Components::SimpleH2O<Scalar>;

public:

    using PermeabilityType = Scalar; // export permeability type

    using RootSpatialParamsDGF<FVGridGeometry, Scalar>::RootSpatialParamsDGF;

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * Additionally, the model of Spery et al is applied, taking the effect of xylem caviation into account
     *
     * \note Kx has units [m^4/(Pa*s)] so we have to divide by the cross-section area
     *       and multiply with a characteristic viscosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
        const SubControlVolume& scv, const ElementSolution& elemSol) const {
        Scalar mu = Water::liquidViscosity(285.15, 1e5); // temperature, pressure
        auto eIdx = this->gridGeometry().elementMapper().index(element);
        Scalar a = this->radius(eIdx);
        Scalar kx = this->kx(eIdx);
        Scalar p = elemSol[0][0];
        Scalar y = std::abs(p-pRef_)/(b_-pRef_);
        double kappa = std::exp(-std::pow(y, c_));
       // Scalar p_cm = toHead_(p);
       // std::cout << "kappa " << kappa << ", " << p_cm << ", " << kx << "\n";
        return  kappa *kx * mu / (M_PI * a * a);
    }

    /**
     * The two parameters of the Weibull function:
     * b [cm], and
     * c [1]
     *
     * k = exp(-(|psi|/b)^c)
     */
    void setParameter(double b, double c) {
        b_ = toPa_(b);
        c_ = c;
    }

private:

    Scalar b_ = 1.e6; // [Pa]
    Scalar c_ = 100.; // [1]

    static constexpr Scalar g_ = 9.81; // cm / s^2
    static constexpr Scalar rho_ = 1.e3; //1.e3; // kg / m^3
    static constexpr Scalar pRef_ = 1.e5; // Pa

    //! cm pressure head -> Pascal
    Scalar toPa_(Scalar ph) const {
        return pRef_ + ph / 100. * rho_ * g_;
    }

   /* //! Pascal -> cm pressure head
    Scalar toHead_(Scalar p) const {
        return (p - pRef_) * 100. / rho_ / g_;
    }*/

};

} // end namespace Dumux

#endif
