// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */

#ifndef DUMUX_COMPOSITIONAL_CYL_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_CYL_LOCAL_RESIDUAL_HH

#include <vector>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Element-wise calculation of the local residual for problems
 *        using compositional fully implicit model.
 */
template<class TypeTag>
class CompositionalLocalResidual: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static constexpr int replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

public:

    using ParentType::ParentType;

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };
        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
            	auto eqIdx = conti0EqIdx + compIdx;
            	Scalar b = problem.bufferPower(scv, volVars, compIdx);
                if (eqIdx != replaceCompEqIdx) {
                    storage[eqIdx] += (volVars.porosity()*volVars.saturation(phaseIdx)+b)
                                      * massOrMoleDensity(volVars, phaseIdx)
                                      * massOrMoleFraction(volVars, phaseIdx, compIdx)*scv.center()[0];
					// std::cout<<"computeStorage, storage[eqIdx]: phaseIdx "<<phaseIdx<<" compIdx "<<compIdx<<" eqIdx "<<eqIdx
					// <<" numComponents "<<numComponents<<" numPhases "<<numPhases <<" storage "<<storage[eqIdx]
					// <<" poro "<<volVars.porosity()<<" saturation "<<volVars.saturation(phaseIdx)<<std::endl;
                }
            }

            // in case one balance is substituted by the total mole balance
			// std::cout<<"replaceCompEqIdx "<<replaceCompEqIdx <<" "<< useTotalMoleOrMassBalance 
			 // <<" "<< numComponents<<std::endl;
            if (useTotalMoleOrMassBalance) {
            	auto compIdx = replaceCompEqIdx- conti0EqIdx;
            	Scalar b = problem.bufferPower(scv, volVars, compIdx);
                storage[replaceCompEqIdx] += massOrMoleDensity(volVars, phaseIdx)
                                             *(volVars.porosity()*volVars.saturation(phaseIdx)+b)*scv.center()[0];
				//std::cout<<"computeStorage, storage[replaceCompEqIdx] "<<compIdx<<" "<<replaceCompEqIdx<<" "<<storage[replaceCompEqIdx]<<std::endl;
            }

            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

		// std::cout<<"computeStorage "<<scv.center()<<" ";
		// for(int i = 0; i< storage.size(); i++)
		// {
			// std::cout<<storage[i]<<" ";
			
		// }std::cout<<std::endl;
		
        return storage;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        NumEqVector flux(0.0);

		//default: use mass fraction => density == pure water density [kg/m^3]
        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
			// see @dumux/porousmediumflow/fluxvariables.hh (advectiveFlux())
			//see @dumux/flux/box/fickslaw.h (flux())
            auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx);
			
			// std::cout<<"diffusiveFluxes ";
			// for(int i = 0; i < diffusiveFluxes.size();i++)
			// {std::cout<<diffusiveFluxes[i]<<" ";}std::cout<<std::endl;
			
			diffusiveFluxes *= scvf.center()[0];//factor according to the position along the 1d model
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding * relative_permability / viscosity
				// relative_permability for the wetting phase of the medium implied by van Genuchten's parameterization.
				// [kg_h2o/m^3_h2o] * kg_a/kg_h2o * [-] / [Pa*s] = [kg_a/m^3_h2o] / [Pa*s]
				
                const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { 
					
					auto d = massOrMoleDensity(volVars, phaseIdx);
					auto fr = massOrMoleFraction(volVars, phaseIdx, compIdx);
					auto m = volVars.mobility(phaseIdx);
					//std::cout<<"upwindterm "<<d<<" "<<fr<<" "<<m<<std::endl;
					return d * fr * m ;
					
				};

                if (eqIdx != replaceCompEqIdx)
				{
					// see @dumux/porousmediumflow/fluxvariables.hh (advectiveFlux())
					//see @dumux/flux/box/darcyslaw.h (flux())
					// gives upwindTerm() * Ks [m/s] * dp [Pa] ==> kg_a / s
					auto advF = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
					//std::cout<<"advF "<<advF<<std::endl;
                    flux[eqIdx] += (advF*scvf.center()[0]); //
				}

                // diffusive fluxes (only for the component balances)
                if(eqIdx != replaceCompEqIdx)
				{
                    flux[eqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                            : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
				}
            }

            // in case one balance is substituted by the total mole balance
            //std::cout<<"replaceCompEqIdx "<<replaceCompEqIdx <<" "<< useTotalMoleOrMassBalance 
			 //<<" "<< numComponents<<std::endl;
            if (useTotalMoleOrMassBalance)
            {
            //DUNE_THROW(Dune::NotImplemented, "localresidual::computeFlux: using useTotalMoleOrMassBalance");
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, phaseIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*volVars.mobility(phaseIdx); };

                flux[replaceCompEqIdx] += (fluxVars.advectiveFlux(phaseIdx, upwindTerm)*scvf.center()[0]); // ;

                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                    flux[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                                       : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            }

            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        EnergyLocalResidual::heatConductionFlux(flux, fluxVars);

		// //std::cout<<"computeFlux "<<scvf.center()<<" ";
		// for(int i = 0; i< flux.size(); i++)
		// {
			// std::cout<<flux[i]<<" ";
			
		// }std::cout<<std::endl;
        return flux;
    }

protected:

    Implementation *asImp_()
    { return static_cast<Implementation *> (this); }

    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

} // end namespace Dumux

#endif
