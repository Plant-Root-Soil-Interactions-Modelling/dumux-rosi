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

#include <dumux/material/components/simpleh2o.hh>

#include <dumux/porousmediumflow/richards/model.hh>

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

// Set the material law
SET_PROP(RichardsParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<VanGenuchten<Scalar>>;
};
}



/*!
 * \ingroup Richards
 * \brief The spatial parameters class for the Richards problem
 */
template<class TypeTag>
class RichardsParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexSet = typename GridView::IndexSet;
    using Water = SimpleH2O<Scalar>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);  //  todo make sure thats the right grid (set somewhere else)
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief Constructor
     */
    RichardsParams(const Problem& problem) : ParentType(problem)  {

    	phi_ = 1; // Richards equation is independent of phi

        Scalar g =  9.81; // std::abs(problem.gravity()[dimWorld-1]);
        Scalar temp =  problem.temperature();
        Scalar pnRef = problem.nonWettingReferencePressure();
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
            // QR
            materialParams_.at(i).setSwr(Qr.at(i)/phi_);
            // materialParams_.at(i).setSwr(0.05);
            // QS
            materialParams_.at(i).setSnr(1.-Qs.at(i)/phi_);
            // materialParams_.at(i).setSnr(0.);
            // ALPHA
            Scalar a = alpha.at(i) * 100.; // from [1/cm] to [1/m]
            materialParams_.at(i).setVgAlpha(a/(rho*g)); //  psi*(rho*g) = p  (from [1/m] to [1/Pa])
            // materialParams_.at(i).setVgAlpha(0.0037);
            // N
            materialParams_.at(i).setVgn(n.at(i));
            // materialParams_.at(i).setVgn(4.7);
            // Kc
            K_.push_back(Kc_.at(i)*mu/(rho*g)); // convert to intrinsic permeability (from the hydraulic conductivity)
            // K_.push_back(5.e-12);
            // Debug
            std::cout << "\nVan Genuchten Parameters are " << Qr.at(i) << ", " << 1.-Qs.at(i)/phi_ <<
                    ", "<< a << ", " << a/(rho*g)<< ", "<< n.at(i) << ", "<< Kc_.at(i)*mu/(rho*g) << ", " << rho << ", "<< g <<"\n";
        }
        // while(std::cin.get()!='\n');
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$. (todo Who uses that?)
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    decltype(auto) permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const {
    	GlobalPosition pos = scv.center();
    	if (more_) {
    		if (pos[2]>1.5) { // hard coded for specific example
    			return K_.at(0);
    		} else {
    			return K_.at(1);
    		}
    	} else {
    		return K_.at(0);
    	}
    }

    /*
     * \brief my own function for defining the (absolute) hydraulic conductivity. [m/s], called by the problem
     */
    const Scalar hydraulicConductivity(const GlobalPosition& pos) const
    {
    	if (more_) {
    		if (pos[2]>1.5) { // hard coded for specific example
    			return Kc_.at(0);
    		} else {
    			return Kc_.at(1);
    		}
    	} else {
    		return Kc_.at(0);
    	}
    }

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate (todo Who uses that?)
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const {
    	return phi_; // porosity*saturation = watercontent
    }


    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \return the material parameters object
     * \param globalPos The position of the center of the element
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
    	if (more_) {
    		if (globalPos[2]>1.5) { // hard coded for specific example
    			return materialParams_.at(0);
    		} else {
    			return materialParams_.at(1);
    		}
    	} else {
    		return materialParams_.at(0);
    	}
    }


private:

    bool more_;
    Scalar phi_;
    std::vector<Scalar> K_; // permeability [m²]
    std::vector<Scalar> Kc_; // hydraulic conductivity [m/s]
    std::vector<MaterialLawParams> materialParams_;

};

} // end namespace Dumux

#endif

