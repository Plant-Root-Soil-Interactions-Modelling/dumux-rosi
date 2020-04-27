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
#ifndef RICHARDS_PARAMETERS_HH
#define RICHARDS_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
// #include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/io/inputfilefunction.hh>

namespace Dumux {

/*!
 * The SpatialParams class of RichardsProblem
 *
 * supports multiple soil layers (in z-direction),
 * with different VG parameters sets
 */
template<class FVGridGeometry, class Scalar>
class RichardsParams : public FVSpatialParams<FVGridGeometry, Scalar, RichardsParams<FVGridGeometry, Scalar>>
{
public:

    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Water = Components::SimpleH2O<Scalar>;
    using MaterialLaw = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    enum { dimWorld = GridView::dimensionworld };

    RichardsParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : FVSpatialParams<FVGridGeometry, Scalar, RichardsParams<FVGridGeometry, Scalar>>(fvGridGeometry)
    {

        /* SimpleH2O is constant in regard to temperature and reference pressure */
        Scalar mu = Water::liquidViscosity(0.,0.); // Dynamic viscosity: 1e-3 [Pa s]
        Scalar rho = Water::liquidDensity(0.,0.);  // Density: 1000 [kg/m³]

        /* Get Van Genuchten parameters from the input file */
        std::vector<Scalar> qr = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qr"); // [1]
        std::vector<Scalar> qs = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Qs"); // [1]
        std::vector<Scalar> alpha = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Alpha");  // [1/cm]
        std::vector<Scalar> n = getParam<std::vector<Scalar>>("Soil.VanGenuchten.N"); // [1]
        kc_ = getParam<std::vector<Scalar>>("Soil.VanGenuchten.Ks"); // hydraulic conductivity [cm/day]
        std::transform(kc_.begin (), kc_.end (), kc_.begin (), std::bind1st(std::multiplies<Scalar>(), 1./100./24./3600.)); // convert from [cm/day] to [m/s]
        homogeneous_ = qr.size()==1; // more than one set of VG parameters?

        phi_.resize(qr.size());
        // Qr, Qs, alpha, and n goes to the MaterialLaw VanGenuchten
        for (int i=0; i<qr.size(); i++) {

            phi_[i] = qs.at(i); // Richards equation is independent of phi [1]
            materialParams_.push_back(MaterialLawParams());
            materialParams_.at(i).setSwr(qr.at(i)/phi_[i]); // Qr
            materialParams_.at(i).setSnr(1.-qs.at(i)/phi_[i]); // Qs
            Scalar a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g_)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            materialParams_.at(i).setVgn(n.at(i)); // N
            k_.push_back(kc_.at(i)*mu/(rho*g_)); // Convert to intrinsic permeability
            // Regularisation parameters
            double eps = 1.e-6; // with 1.e-9 benchmark 3 does not work anymore (and everything becomes slower)
            materialParams_.at(i).setPcLowSw(eps);
            materialParams_.at(i).setPcHighSw(1. - eps);
            materialParams_.at(i).setKrnLowSw(eps);
            materialParams_.at(i).setKrwHighSw(1 - eps);
        }
        layerIdx_ = Dumux::getParam<int>("Soil.Grid.layerIdx", 1);
        layer_ = InputFileFunction("Soil.Layer", "Number", "Z", layerIdx_, 0); // [1]([m])

        // std::cout << "RichardsParams created: homogeneous " << homogeneous_ << " " << "\n" << std::endl;
    }

    /*!
     * \brief \copydoc FVGridGeometry::porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const {
        return phi_.at(index_(element));
    }

    /*!
     * \brief \copydoc FVGridGeometry::porosity
     * simper interface
     */
    Scalar porosity(const Element& element) const {
        return phi_.at(index_(element));
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::permeability
     * [m^2]\
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
        const SubControlVolume& scv, const ElementSolution& elemSol) const {
        return permeability(element);
    }

    //! simpler interface
    PermeabilityType permeability(const Element& element) const {
        return k_.at(index_(element));
    }

    /*
     * \brief Hydraulic conductivities [m/s], called by the problem for conversions
     */
    const Scalar hydraulicConductivity(const Element& element) const {
        return kc_.at(index_(element));
    }

    //! set of VG parameters for the element
    const MaterialLawParams& materialLawParams(const Element& element) const {
        return materialParams_.at(index_(element));
    }

    /*!
     * \brief \copydoc FVSpatialParamsOneP::materialLawParams
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
        const SubControlVolume& scv,
        const ElementSolution& elemSol) const {
        return materialParams_.at(index_(element));
    }

    //! pointer to the soils layer input file function
    InputFileFunction* layerIFF() {
        return &layer_;
    }

private:

    //! returns the index of the soil layer
    size_t index_(const Element& element) const {
        if (homogeneous_) {
            return 0;
        } else {
            auto eIdx = this->fvGridGeometry().elementMapper().index(element);
            Scalar z = element.geometry().center()[dimWorld - 1];
            //std::cout << z << "\n";
            return size_t(layer_.f(z, eIdx)-1); // layer number starts with 1 in the input file
        }
    }

    std::vector<Scalar> phi_; // porosity

    bool homogeneous_; // soil is homogeneous
    InputFileFunction layer_;
    int layerIdx_; // index of layer data within the grid

    std::vector<Scalar> k_; // permeability [m²]
    std::vector<Scalar> kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

    static constexpr Scalar g_ = 9.81; // cm / s^2

};

} // end namespace Dumux

#endif
