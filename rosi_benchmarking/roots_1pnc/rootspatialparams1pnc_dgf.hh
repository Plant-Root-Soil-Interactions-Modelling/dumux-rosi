// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
#ifndef DUMUX_ROOT_SPATIALPARAMS1PNC_DGF_HH
#define DUMUX_ROOT_SPATIALPARAMS1PNC_DGF_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <cmath>

#include <dumux/io/inputfilefunction.hh> // in dumux-rosi

#include "../roots_1p/rootspatialparams_dgf.hh"

namespace Dumux {

/*!
 * \brief Root spatial parameters class for grid files.
 *
 * use initParameters to initialize the class with data from the grid file
 *
 */
template<class FVGridGeometry, class Scalar>
class RootSpatialParams1pncDGF : public RootSpatialParamsDGF<FVGridGeometry, Scalar> {

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
        auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        Scalar a = this->radius(eIdx);
        Scalar kx = this->kx(eIdx);
        Scalar p = elemSol[0][0];

        std::cout << "pressure? "<< p << ", " ; // for debugging

        Scalar y = -std::abs(p)/b_;
        Scalar k = std::exp(std::exp(c_*std::log(y))); // std::exp(c*std::log(y)) == y^c

         return kx * mu / (M_PI * a * a);
    }

    /**
     * The two parameters of the Weibull function:
     * b [Pa] or [cm] TODO, and
     * c [1]
     *
     * k = exp(-(|psi|/b)^c)
     */
    void setParameter(double b, double c)
    {
        b_ = b;
        c_ = c;
    }

private:

    Scalar b_ = 1.; // [Pa] or shall we convert to [cm] ?
    Scalar c_ = 1.; // [1]

};

} // end namespace Dumux

#endif
