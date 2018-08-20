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
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsLensProblem
 */
template<class FVGridGeometry, class Scalar>
class RichardsParams : public FVSpatialParams<FVGridGeometry, Scalar, RichardsParams<FVGridGeometry, Scalar>>
{
    using ThisType = RichardsParams<FVGridGeometry, Scalar>;
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
        for (int i=0; i<Qr.size(); i++)
        {
            materialParams_.push_back(MaterialLawParams());
            materialParams_.at(i).setSwr(Qr.at(i)/phi_); // Qr
            materialParams_.at(i).setSnr(1.-Qs.at(i)/phi_); // Qs
            Scalar a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g_)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            materialParams_.at(i).setVgn(n.at(i)); // N
            K_.push_back(Kc_.at(i)*mu/(rho*g_)); // Convert to intrinsic permeability
        }
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::permeability
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element, const SubControlVolume& scv, const ElementSolution& elemSol) const
    {
        if (homogeneous_)
        {
            return K_.at(0);
        }
        else // TODO we need the grid manager to look up stuff
        {
            GlobalPosition pos = scv.center();
            if (pos[dimWorld-1]>1.5)
            { // hard coded for specific example
                return K_.at(0);
            } else
            {
                return K_.at(1);
            }
        }
    }

    /*
     * \brief Hydraulic conductivites [m/s], called by the problem for conversions
     */
    const Scalar hydraulicConductivity(const GlobalPosition& pos) const
    {
        if (homogeneous_)
        {
            return Kc_.at(0);
        }
        else // TODO we need the grid manager to look up stuff
        {
            if (pos[dimWorld-1]>1.5)
            { // hard coded for specific example
                return Kc_.at(0);
            } else
            {
                return Kc_.at(1);
            }
        }
    }


    /*!
     * \brief \copydoc FVGridGeometry::porosity
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return phi_;
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::materialLawParams
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        return materialLawParamsAtPos(globalPos);
    }

    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (homogeneous_)
        {
            return materialParams_.at(0);
        }
        else // TODO we need the grid manager to look up stuff
        {
            if (globalPos[dimWorld-1]>1.5) { // hard coded for specific example
                return materialParams_.at(0);
            } else {
                return materialParams_.at(1);
            }
        }
    }

private:

    bool homogeneous_; // soil is homogeneous

    Scalar phi_; // porosity

    std::vector<Scalar> K_; // permeability [m²]
    std::vector<Scalar> Kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

    static constexpr Scalar g_ = 9.81; // cm / s^2 (for type conversions)

};

} // end namespace Dumux

#endif
