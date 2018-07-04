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
#ifndef RICHARDS_PARAMS
#define RICHARDS_PARAMS

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>



namespace Dumux //
{
/*!
 * \ingroup RichardsTests
 * \brief spatial parameters for the RichardsLensProblem
 */
// forward declaration
template<class TypeTag>
class RichardsParams;



namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsParams);

// Set the spatial parameters
SET_TYPE_PROP(RichardsParams, SpatialParams, RichardsParams<TypeTag>);

//// Set the material law
//SET_PROP(RichardsParams, MaterialLaw)
//{
//private:
//    // define the material law which is parameterized by effective
//    // saturations
//    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//public:
//    // define the material law parameterized by absolute saturations
//    using type = EffToAbsLaw<VanGenuchten<Scalar>>;
//};
}


/*!
 * \ingroup Richards
 * \brief The spatial parameters class for the Richards problem
 */
template<class TypeTag>
class RichardsParams : public FVSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), typename GET_PROP_TYPE(TypeTag, Scalar), RichardsParams<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, RichardsParams<TypeTag>>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexSet = typename GridView::IndexSet;
    using Water = Components::SimpleH2O<Scalar>;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    using EffectiveLaw = EffToAbsLaw<VanGenuchten<Scalar>>;

public:

    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    /*!
     * \brief Constructor
     */
    RichardsParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry)  {
    	phi_ = 1; // Richards equation is independent of phi
        Scalar g =  9.81; // std::abs(problem.gravity()[dimWorld-1]);
        Scalar temp =  0.;
        Scalar pnRef = 0.;
        Scalar mu = Water::liquidViscosity(temp,pnRef); // h2o: 1e-3 Pa·s Dynamic viscosity of water(independent of temp and p)
        Scalar rho = Water::liquidDensity(temp,pnRef);  // h2o: 1000 kg/m³ Density of water(independent of temp and p)

        // get van genuchten parameters from the input file
        std::vector<Scalar> Qr = getParam<std::vector<Scalar>>("VanGenuchten.Qr");
        std::vector<Scalar> Qs = getParam<std::vector<Scalar>>("VanGenuchten.Qs");
        std::vector<Scalar> alpha = getParam<std::vector<Scalar>>("VanGenuchten.alpha");
        std::vector<Scalar> n = getParam<std::vector<Scalar>>("VanGenuchten.n");
        Kc_ = getParam<std::vector<Scalar>>("VanGenuchten.Ks"); // hydraulic conductivity
        more_ = Qr.size()>1; // more than one set of VG parameters?

        // Qr, Qs, alpha, and n goes to the MaterialLaw VanGenuchten
        for (int i=0; i<Qr.size(); i++) {
            materialParams_.push_back(MaterialLawParams());
            materialParams_.at(i).setSwr(Qr.at(i)/phi_); // Qr
            materialParams_.at(i).setSnr(1.-Qs.at(i)/phi_); // Qs
            Scalar a = alpha.at(i) * 100.; // alpha, from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            materialParams_.at(i).setVgn(n.at(i)); // n
            K_.push_back(Kc_.at(i)*mu/(rho*g)); // Kc, convert to intrinsic permeability (from the hydraulic conductivity)
        }
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$. override from FVSpatialParamsOneP
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const {
    	 return K_.at(getDI(element));
    }

    /*!
     * The (absolute) hydraulic conductivity. [m/s] called by the problem class
     */
    const Scalar hydraulicConductivity(const Element &element) const
    {
        return Kc_.at(getDI(element));
    }

   /*!
    * \brief Define the porosity in [-]. ovveride from FVSpatialParamsOneP
    */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const {
    	return phi_; // porosity*saturation = watercontent
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const    {
    	return materialParams_.at(getDI(element));
    }

    // why
    const MaterialLawParams& materialLawParams(const Element& element) const    {
    	return materialParams_.at(getDI(element));
    }


private:

    /**
     * returns the domain index
     */
    int getDI(const Element &element) const {
        if (more_) {
//            int i = (GridCreator::parameters(element)).at(1)-1; // starting from 1
//            return i; // TODO update
            return 0;
        } else {
            return 0;
        }
    }

    bool more_;
    Scalar phi_;
    std::vector<Scalar> K_; // permeability [m²]
    std::vector<Scalar> Kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

};

} // end namespace Dumux

#endif

