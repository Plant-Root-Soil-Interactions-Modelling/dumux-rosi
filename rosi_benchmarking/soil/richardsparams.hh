// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup RichardsTests
 * \brief spatial parameters for the RichardsLensProblem
 */
#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh> // TODO kick
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsLensProblem
 */
template<class Grid, class FVGridGeometry, class Scalar>
class RichardsParams : public FVSpatialParams<FVGridGeometry, Scalar, RichardsParams<Grid, FVGridGeometry, Scalar>>
{
    using ThisType = RichardsParams<Grid, FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Water = Components::SimpleH2O<Scalar>;

    enum {
        dimWorld = GridView::dimensionworld
    };

public:

    enum GridParameterIndex {
        layerNumber = 0
    };

    using MaterialLaw = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    RichardsParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {

        phi_ = 1; // Richards equation is independent of phi

        /* SimpleH2O is constant in regard to temperature and reference pressure */
        Scalar mu = Water::liquidViscosity(0.,0.); // Dynamic viscosity: 1e-3 Pa·s
        Scalar rho = Water::liquidDensity(0.,0.);  // Density: 1000 kg/m³

        // Get Van Genuchten parameters from the input file
        std::vector<Scalar> Qr = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qr");
        std::vector<Scalar> Qs = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qs");
        std::vector<Scalar> alpha = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Alpha");
        std::vector<Scalar> n = getParam<std::vector<Scalar>>("Soil.VanGenuchten.N");
        Kc_ = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Ks"); // hydraulic conductivity
        homogeneous_ = Qr.size()==1; // more than one set of VG parameters?
        // Qr, Qs, alpha, and n goes to the MaterialLaw VanGenuchten
        for (int i=0; i<Qr.size(); i++) {
            materialParams_.push_back(MaterialLawParams());
            materialParams_.at(i).setSwr(Qr.at(i)/phi_); // Qr
            materialParams_.at(i).setSnr(1.-Qs.at(i)/phi_); // Qs
            Scalar a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g_)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            materialParams_.at(i).setVgn(n.at(i)); // N
            K_.push_back(Kc_.at(i)*mu/(rho*g_)); // Convert to intrinsic permeability
        }

        // check if there is a layer look up table
        try {
            z_ = getParam<std::vector<Scalar>>("Soil.LayerZ");
            layer_ = getParam<std::vector<Scalar>>("Soil.LayerNumber");
            layerTable_ = true;
        } catch(std::exception& e) {
            layerTable_ = false;
        }

    }

    /*!
     * \brief Grid manager is needed for layer number look up
     */
    void setGridManager(GridManager<Grid>* gm)
    {
        gridManager_ = gm;
    }

    /*!
     * \brief \copydoc FVGridGeometry::porosity
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return phi_;
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::permeability
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element, const SubControlVolume& scv, const ElementSolution& elemSol) const
    {
        return K_.at(index_(element));
    }

    /*
     * \brief Hydraulic conductivites [m/s], called by the problem for conversions
     */
    const Scalar hydraulicConductivity(const Element& element) const
    {
        return Kc_.at(index_(element));
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::materialLawParams
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
        const SubControlVolume& scv,
        const ElementSolution& elemSol) const
    {
        return materialParams_.at(index_(element));
    }

    const MaterialLawParams& materialLawParams2(const Element& element) const
    {
        return materialParams_.at(index_(element));
    }

    //! 1d table look up: xx is ascending, returns the index i , so that x>=xx[i] and x<xx[i+1]
    static size_t locate(Scalar x, const std::vector<Scalar>& xx)
    {
        unsigned int jr,jm,jl;
        jl = 0;
        jr = xx.size();
        while (jr-jl > 1) {
            jm=(jr+jl) >> 1; // thats a divided by two
            if (x >= xx[jm])
                jl=jm;
            else
                jr=jm;
        }
        return jl; // left index
    }

    //! returns linearly interpolated values of a 1-D function at specific query point x. Vector xx contains the sample points, and vv contains the corresponding values
    static Scalar interp1(Scalar x, const std::vector<Scalar>& vv, const std::vector<Scalar>& xx)
    {
        size_t i = locate(x, xx);
        Scalar t = (x - xx[i])/(xx[i+1] - xx[i]);
        t = std::min(std::max(t,0.),1.);
        Scalar v = vv[i]*(1.-t) + vv[i+1]*t;
        return v;
    }


private:

    //! returns the index of the soil layer
    size_t index_(const Element& element) const {
        if (homogeneous_) {
            return 0;
        } else {
            if (layerTable_) { // obtain from look up table
                double z = element.geometry().center()[dimWorld-1];
                return (size_t) (interp1(z,layer_, z_) + 0.5); // round
            } else { // obtain from grid
                return (size_t)( gridManager_->getGridData()->parameters(element).at(layerNumber)+0.5 );
            }
        }
    }

    bool homogeneous_; // soil is homogeneous
    bool layerTable_; // layer data in look up table
    std::vector<Scalar> z_;
    std::vector<Scalar> layer_;
    Scalar phi_; // porosity

    GridManager<Grid>* gridManager_ = nullptr;

    std::vector<Scalar> K_; // permeability [m²]
    std::vector<Scalar> Kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)

};

} // end namespace Dumux

#endif
