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

#ifndef DUMUX_COMPOSITIONAL_10C_LOCAL_RESIDUAL_HH
#define DUMUX_COMPOSITIONAL_10C_LOCAL_RESIDUAL_HH

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
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
	using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
	using FVGridGeometry = typename Problem::GridGeometry;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int numFluidPhases = ModelTraits::numFluidPhases();
    static constexpr int numFluidComponents = ModelTraits::numFluidComponents();
    static constexpr int numSolidComps =  ModelTraits::numSolidComponents();
    static constexpr int numInertSolidComps =  ModelTraits::numInertSolidComps();
    static constexpr int soilIdx =  SolidSystem::mainCompIdx;
    static constexpr bool useMoles = ModelTraits::useMoles();
	static constexpr int dimWorld = GridView::dimensionworld;

    enum { conti0EqIdx = Indices::conti0EqIdx };

    //! The index of the component balance equation that gets replaced with the total mass balance
    static constexpr int replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
    static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numFluidComponents;

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
		double pos0 = 1.;
        if(dimWorld==1)//cylindrical model
		{
			pos0 = scv.center()[0];
		}
		
		auto storage = computeSolidStorage(problem, scv, volVars);//compute storage for solid phase

        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };
        // compute storage term of all components within all phases
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numFluidComponents; ++compIdx)
            {
            	auto eqIdx = conti0EqIdx + compIdx;
                if (eqIdx != replaceCompEqIdx) {
						int dofIndex = scv.dofIndex();
						Scalar b = problem.bufferPower(dofIndex, volVars, compIdx, scv);
						
					// mol solute / m3 space
                    storage[eqIdx] += std::max(0.,//manually force it to remaine above 0.
                                        volVars.porosity()*volVars.saturation(phaseIdx) //m3 liquide / m3 space
									  * b // [-]
                                      * massOrMoleDensity(volVars, phaseIdx) // mol liquid / m3 liquid
                                      * massOrMoleFraction(volVars, phaseIdx, compIdx) 
									  * pos0);// mol solute / mol liquid 
									  
                }
            }

            if (useTotalMoleOrMassBalance) {//just for water phase
                storage[replaceCompEqIdx] += massOrMoleDensity(volVars, phaseIdx)
											* (volVars.porosity()*volVars.saturation(phaseIdx) 
											* pos0);
            }

            //! The energy storage in the fluid phase with index phaseIdx
            EnergyLocalResidual::fluidPhaseStorage(storage, scv, volVars, phaseIdx);
        }

        //! The energy storage in the solid matrix
        EnergyLocalResidual::solidPhaseStorage(storage, scv, volVars);

        return storage;
    }
	
	
    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * We consider the volume-average here (e.g. phase mass inside a
     * sub control volume divided by the volume). The volume is multiplied
     * onto it afterwards in the local residual of the respective spatial
     * discretization scheme.
     *
     * \param problem The problem (Initial/Boundary conditions...) to be solved
     * \param scv The sub-control volume of the finite volume grid
     * \param volVars The volume variables (primary/secondary variables) in the scv
     * \return Amount per volume of the conserved quantities
     */
    NumEqVector computeSolidStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
		double pos0 = 1.;
        if(dimWorld==1)//cylindrical model
		{
			pos0 = scv.center()[0];
		}
		
        NumEqVector storage(0.0);
		if(numSolidComps > numInertSolidComps)
		{
                if(problem.verbose_local_residual)
                {
                     std::cout<<"solid storage ";
                }
			//soil mass or mol density (mol soil C / m3 soil)
		double massOrMoleDensity =  useMoles ? volVars.solidComponentMolarDensity(soilIdx) : volVars.solidComponentDensity(soilIdx);
			//component mole fraction (mol comp C/mol soil C)
        const auto massOrMoleFraction= [](const auto& volVars, const int compIdx)
        { double mOMF = useMoles ? volVars.solidMoleFraction(compIdx) : volVars.solidMassFraction(compIdx); return mOMF; };

			// compute storage term of all components within all fluid phases
			for (int sCompIdx = 0; sCompIdx < numSolidComps - numInertSolidComps; ++sCompIdx)
			{
				auto eqIdx = Indices::conti0EqIdx + numFluidComponents + sCompIdx ;
				// mol comp / m3 space
				storage[eqIdx] += std::max(0.,//manually force it to remaine above 0.
                                (1 - volVars.porosity())// m3 solide / m3 space
								* massOrMoleDensity	//mol solid / m3 solide
								* massOrMoleFraction(volVars, sCompIdx) // mol comp / mol solid
								* pos0); // mol comp / mol solid
								
                if(problem.verbose_local_residual)
                {
                     std::cout<<", ("<<storage[eqIdx]<<", "<<  massOrMoleFraction(volVars, sCompIdx)<<")";
                }
			}
            if(problem.verbose_local_residual)
                {
                     std::cout<<std::endl;
                }
		}
		
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
		double pos0 = 1.;
		double scvf_area = 1.;
        if(dimWorld==1)//cylindrical model
		{
			pos0 = scvf.center()[0];
			scvf_area = 2 * M_PI * scvf.center()[0] * problem.segLength;//m2
		}
		
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        // get upwind weights into local scope
        NumEqVector flux = computeSolidFlux();//(0.0);

		//default: use mass fraction => density == pure water density [kg/m^3]
        const auto massOrMoleDensity = [](const auto& volVars, const int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        const auto massOrMoleFraction= [](const auto& volVars, const int phaseIdx, const int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // advective fluxes
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
        {
			// see @dumux/porousmediumflow/fluxvariables.hh (advectiveFlux())
			//see @dumux/flux/box/fickslaw.h (flux())
            auto diffusiveFluxes = fluxVars.molecularDiffusionFlux(phaseIdx) * pos0;
			
			
			
            for (int compIdx = 0; compIdx < numFluidComponents; ++compIdx)
            {
                // get equation index
                const auto eqIdx = conti0EqIdx + compIdx;

                // the physical quantities for which we perform upwinding * relative_permability / viscosity
				// relative_permability for the wetting phase of the medium implied by van Genuchten's parameterization.
				// [kg_h2o/m^3_h2o] * kg_a/kg_h2o * [-] / [Pa*s] = [kg_a/m^3_h2o] / [Pa*s]
				
                const auto upwindTerm = [&massOrMoleDensity, &massOrMoleFraction, phaseIdx, compIdx] (const auto& volVars)
                { 
					return massOrMoleDensity(volVars, phaseIdx) * massOrMoleFraction(volVars, phaseIdx, compIdx) * volVars.mobility(phaseIdx) ; 
				};
				
				Scalar b = 1;

                if ((eqIdx != replaceCompEqIdx)&&(compIdx != problem.mucilIdx))//&&(problem.doSoluteFlow))
				{
					
					// see @dumux/porousmediumflow/fluxvariables.hh (advectiveFlux())
					//see @dumux/flux/box/darcyslaw.h (flux())
					// gives upwindTerm() * Ks [m/s] * dp [Pa] ==> kg_a / s
                    flux[eqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm)/b*pos0; //
				}

                // diffusive fluxes (only for the component balances)
                if(eqIdx != replaceCompEqIdx)
				{
                    flux[eqIdx] += (( useMoles ? diffusiveFluxes[compIdx]
                                            : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx))/b);
				}
				
            }

            // // in case one balance is substituted by the total mole balance
            if (useTotalMoleOrMassBalance)//water phase
            {
                // the physical quantities for which we perform upwinding
                const auto upwindTerm = [&massOrMoleDensity, phaseIdx] (const auto& volVars)
                { return massOrMoleDensity(volVars, phaseIdx)*volVars.mobility(phaseIdx); };

                flux[replaceCompEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm) * pos0; // ;

                for(int compIdx = 0; compIdx < numFluidComponents; ++compIdx)
                    flux[replaceCompEqIdx] += useMoles ? diffusiveFluxes[compIdx]
                                                       : diffusiveFluxes[compIdx]*FluidSystem::molarMass(compIdx);
            }

            //! Add advective phase energy fluxes. For isothermal model the contribution is zero.
            // EnergyLocalResidual::heatConvectionFlux(flux, fluxVars, phaseIdx);
        }

        //! Add diffusive energy fluxes. For isothermal model the contribution is zero.
        // EnergyLocalResidual::heatConductionFlux(flux, fluxVars);
		
		if(scvf.numOutsideScvs() > 1) //can be 0 (boundary) or 1
		{
			DUNE_THROW(Dune::InvalidStateException, "wrong value for numOutsideScvs");
		}else{
			problem.setFaceFlux(flux/pos0*scvf_area, scvf.index());
		}
			
		return flux;
    }


    NumEqVector computeSolidFlux() const
    {
        NumEqVector flux(0.0);
		if(numSolidComps > numInertSolidComps)
		{
			// compute storage term of all components within all fluid phases
			for (int sCompIdx = 0; sCompIdx < numSolidComps - numInertSolidComps; ++sCompIdx)
			{
				auto eqIdx = Indices::conti0EqIdx + numFluidComponents + sCompIdx ;
				flux[eqIdx] += 0.;
			}
		}
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
