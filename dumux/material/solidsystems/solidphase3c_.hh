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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::LiquidPhaseThreeC
 */
#ifndef DUMUX_SOLID_THREEC_PHASE_HH
#define DUMUX_SOLID_THREEC_PHASE_HH

#include <cassert>
#include <limits>

#include <dune/common/exceptions.hh>
//#include <dumux/material/solidsystems/base.hh>
//#include <dumux/material/binarycoefficients/h2o_constant.hh>
#include <dumux/io/name.hh>

namespace Dumux {
namespace SolidSystems {


/*!
 * \ingroup Solidsystems
 * \brief A Solid phase consisting of a two components,
 *        a main component and a conservative tracer component
 */
template <class Scalar, class MainComponent, class SecondComponent, class ThirdComponent>
class SolidPhaseThreeC
: public Base<Scalar, SolidPhaseThreeC<Scalar, MainComponent, SecondComponent, ThirdComponent> >
{
    using ThisType = SolidPhaseThreeC<Scalar, MainComponent, SecondComponent, ThirdComponent>;
    using Base = Dumux::SolidSystems::Base<Scalar, ThisType>;
    //file only for h2o + one componant. better re-write own liquidDiffCoeff() function
	//using BinaryCoefficients = BinaryCoeff::H2O_Component<Scalar, SecondComponent>;//??

public:
    using ParameterCache = NullParameterCache;

    static constexpr int numPhases = 1; //!< Number of phases in the Solid system
    static constexpr int numComponents = 3; //!< Number of components in the Solid system

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the only phase

    static constexpr int comp0Idx = 0; //!< index of the frist component
    static constexpr int comp1Idx = 1; //!< index of the second component
    static constexpr int comp2Idx = 2; //!< index of the 3rd component
    static constexpr int mainCompIdx = comp0Idx; //!< index of the main component
    static constexpr int secondCompIdx = comp1Idx; //!< index of the secondary component
    static constexpr int thirdCompIdx = comp2Idx; //!< index of the secondary component

    /*!
    * \brief Initialize the Solid system's static parameters generically
    */
    static void init() {}

    /****************************************
     * Solid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a Solid phase
     *
     * \param phaseIdx The index of the Solid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return IOName::solidPhase(); }

    /*!
     * \brief Returns whether the component is inert (doesn't react)
     */
    static constexpr bool isInert()
    {  return true; }

    /*!
     * \brief Returns whether the Solids are miscible
     * \note There is only one phase, so miscibility makes no sense
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    { 
		std::string componentName_;
		switch(compIdx)
		{
			case comp0Idx:
			{
				componentName_ = MainComponent::name();
				break;
			}
			case comp1Idx:
			{
				componentName_ = SecondComponent::name();
				break;
			}
			case comp2Idx:
			{
				componentName_ = ThirdComponent::name();
				break;
			}
			default:
			{
				DUNE_THROW(Dune::InvalidStateException, "liquidphase3c::componentName: compIdx "<<compIdx<<" not recognised");
			}				
		}		
		return componentName_; 
	}

    /*!
     * \brief A human readable name for the Solid system.
     */
    static std::string name()
    { return "SolidPhaseThreeC"; }

    /*!
     * \brief Returns whether the Solid is gaseous
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns true if and only if a Solid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a Solid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the Solid composition. This assumption is true
     * if only a single component is involved. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the Solid phase to consider
     */
    static bool isIdealMixture(int phaseIdx = 0)
    { return true; }

    /*!
     * \brief Returns true if the Solid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return false; }

    /*!
     * \brief Returns true if the Solid is assumed to be an ideal gas
     */
    static bool isIdealGas(int phaseIdx = 0)
    { return false; /* we're a liquid! */ }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx)
    {  
		Scalar molarMass_;
		switch(compIdx)
		{
			case comp0Idx:
			{
				molarMass_ = MainComponent::molarMass();
				break;
			}
			default:
			{
				DUNE_THROW(Dune::InvalidStateException, "liquidphase3c::molarMass: compIdx "<<compIdx<<" not recognised");
			}
		}		
		//DUNE_THROW(Dune::InvalidStateException, "liquidphase3c::molarMass: do not use molar mass");
		return molarMass_; 
	}

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the phase at a given pressure and temperature.
     */
    static Scalar density(Scalar temperature)
    {  return MainComponent::solidDensity(temperature); }

    using Base::density;
    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the phase at a given pressure and temperature.
     */
    template <class SolidState>
    static Scalar density(const SolidState &solidState)
    {
        return density(solidState.temperature());		  
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        return MainComponent::solidMolarDensity(T, p);
    }


   /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    static Scalar thermalConductivity(Scalar temperature)
    { return MainComponent::solidThermalConductivity(temperature); }

    /*!
     * \brief Thermal conductivity of the solid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class SolidState>
    static Scalar thermalConductivity(const SolidState &solidState)
    {
        return thermalConductivity(solidState.temperature());
    }

    /*!
     * \brief Specific isobaric heat capacity of the solid \f$\mathrm{[J/(kg K)]}\f$.
     */
    static Scalar heatCapacity(Scalar temperature)
    { return MainComponent::solidHeatCapacity(temperature); }

    /*!
     * \brief Specific isobaric heat capacity of the solid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class SolidState>
    static Scalar heatCapacity(const SolidState &solidState)
    {
        return heatCapacity(solidState.temperature());
    }
    
};

} // namespace FluidSystems

} // namespace Dumux

#endif
